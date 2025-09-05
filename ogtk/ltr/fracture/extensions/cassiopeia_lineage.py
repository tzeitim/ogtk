import polars as pl
from pathlib import Path
from typing import Set,Dict, Optional, Type, List, Tuple
from enum import Enum
from .base import PostProcessorExtension
from .registry import extension_registry
from .config import ExtensionConfig
from ..pipeline.types import StepResults,FractureXp
from dataclasses import dataclass, field
from ogtk.utils.log import CustomLogger
from ogtk.utils.general import fuzzy_match_str

def plug_cassiopeia(
        ldf: pl.LazyFrame,
        ann_intbc_mod: pl.DataFrame,
        workdir: Path|str ='.',
        logger: None|CustomLogger=  None,
        barcode_interval: List|Tuple = (0, 7),
        cutsite_locations: List =  [40, 67, 94, 121, 148, 175, 202, 229, 256, 283],
        cutsite_width: int = 12,
        context: bool = True,
        context_size: int = 50,
        ) -> StepResults:
    """ 
    readName - A unique identifier for each row/sequence
    cellBC - The cell barcode
    UMI - The UMI (Unique Molecular Identifier)
    readCount - The number of reads for this sequence
    seq - The actual sequence to be aligned

    """
    import cassiopeia as cas

    allele_params = {
        'barcode_interval': barcode_interval,
        'cutsite_locations': cutsite_locations,
        'cutsite_width': cutsite_width, 
        'context': context,
        'context_size': context_size,
    }

    cass_ldf = (
            ldf.with_columns(
            readName=pl.col('umi'),
            cellBC=pl.col('umi'), 
            UMI=pl.col('umi'),
            readCount=pl.col('reads'),
            seq=pl.col('intBC')+pl.col('contig'),
            )
            .select('readName', 'cellBC', 'UMI', 'readCount', 'seq', 'intBC', 'sbc')
            #TODO: change `how`?
            .join(ann_intbc_mod.lazy(), left_on='intBC', right_on='intBC', how='inner')
    )

    res = []
    umi_tables = []
    nones = 0


    for (intBC, mod, sbc), queries in cass_ldf.collect().partition_by('intBC', 'mod', 'sbc', as_dict=True).items():
        if mod is None:
            nones+=1
        else:
            if logger:
                logger.debug(f"{mod=} {intBC=} {queries.shape=}")

            umi_table = cas.pp.align_sequences(
                 queries=queries.to_pandas(),
                 ref_filepath=f'{workdir}/{mod}.fasta',
                 n_threads=1,
                 method='global'
             )
    
            allele_table =  cas.pp.call_alleles(
                            umi_table,
                            ref_filepath = f'{workdir}/{mod}.fasta',
                            **allele_params,
                        )

            #replace intBC for real intBC since Cassiopeia does something I don't understand yet
            # include colums dropped by cass?
            allele_table['intBC'] = intBC
            allele_table['mod'] = mod
            allele_table['sbc'] = sbc

            umi_tables.append(umi_table.copy())
            #self.xp.logger.info(f"found {umi_table['intBC'].n_unique()} integrations for {mod}")

            res.append(
                pl.DataFrame(allele_table).with_columns(mod=pl.lit(mod), mols=allele_table.shape[0])
                )

    return StepResults(
            results={"alleles_pl" : pl.concat(res), "alleles_pd" : umi_tables}
    )

def save_ref_to_fasta(refs, out_dir='.', field='mod') -> None : 
    """ write fastas in current dir"""
    for i in refs.get_column(field):
        with open(f'{out_dir}/{i}.fasta', 'wt') as out:
            out.write(
            refs.filter(pl.col(field)==i).dna.to_fasta(read_id_col=field, read_col='seq').get_column('seq_fasta')[0]
            )

def kmer_classify_cassettes(ldf, refs: pl.DataFrame, K: int = 25) -> StepResults:
    """ Annotates contigs based on the top occurring mod based on kmer matches """

    k_ref = refs.kmer.explode_kmers(k=K, seq_field='seq')

    cont_k = (
        ldf
        .kmer.explode_kmers(k=K, seq_field='contig', only_unique=False)
        .filter(pl.col('kmer').is_in(set(k_ref.get_column('kmer').to_list())))
        .join(k_ref.drop('seq').lazy(), left_on='kmer', right_on='kmer', how='inner')
        .group_by('intBC', 'mod', 'sbc')
                .agg(fps=pl.col('umi').n_unique())
        .group_by('intBC', 'mod')
                .agg(fps=pl.col('fps').sum())
        .sort('fps', descending=True)
        .group_by('intBC', maintain_order=True)
        .first()
        .collect(engine="streaming")
        )
    # TODO return how many intBCs didn't get a mod
    # TODO return the number of ties 
    return StepResults(results={'ann_intbc_mod':cont_k},
                       metrics={'n_ann_intbc':cont_k.height})

def parse_contigs(
        ldf: pl.LazyFrame,
        int_anchor1: str,
        int_anchor2 : str,
        sbc_dict: Optional[Dict] = None,
        annotation:str = 'sbc',
        ) ->StepResults:
    """Extracts integration barcodes (intBC) and sample barcodes (sbc)."""
    
    # 1. Extract/use sbcs
    if sbc_dict is not None:
        ldf = ldf.with_columns(pl.col('sbc').replace(sbc_dict).alias(annotation))
    
    # 2. Extract intBC
    ldf = (
        ldf   
        .with_columns(
            pl.col('contig').str
            .extract(f'{int_anchor1}(.+?){int_anchor2}', 1).alias('intBC')
        )
    )
    
    return StepResults(results={'ldf':ldf},
                       metrics={'n_parsed_contigs':ldf.select(pl.len()).collect().item()})


def generate_refs_from_fasta(refs_fasta_path: str|Path, anchor1: str, anchor2: str) -> pl.DataFrame:
    """ Given a FASTA file with references it generates a trimmed version of the cassettes for aligment"""

    refs =  (
            pl.read_csv(refs_fasta_path,
                        has_header=False,
                        )
            .unstack(step=2, how='horizontal')
            .rename({"column_1_0":"mod", "column_1_1":"seq"})
            .with_columns(
                pl.col('mod').str.extract(">(v.+?)_",1).str.replace("v","mod_"),
                pl.col('seq').str.to_uppercase())
            # trim sequence
            .with_columns(
                pl.col('seq').str.replace(f".+?({anchor1})", anchor1).str.replace(f"({anchor2}).+?$", anchor2)
                )
            )

    return refs

@dataclass
class CassiopeiaConfig(ExtensionConfig):
    # Required fields
    int_anchor1: str
    int_anchor2: str
    
    # Optional fields with defaults
    sbc_dict: Optional[Dict] = None
    annotation: str = 'sbc'
    refs_fasta_path: Optional[str] = None
    anchor1: Optional[str] = None
    anchor2: Optional[str] = None
    
    # Runtime parameters (set during processing)
    ann_intbc_mod: Optional[pl.DataFrame] = None

    barcode_interval: Tuple[int, int] = (0, 7)
    cutsite_locations: List[int] = field(default_factory=lambda: [40, 67, 94, 121, 148, 175, 202, 229, 256, 283])
    cutsite_width: int = 12
    context: bool = True
    context_size: int = 50
    

class CassiopeiaStep(Enum):
    """Steps within the Cassiopeia lineage extension"""
    PARSE_CONTIGS = 'parse_contigs'
    CLASSIFY_CASSETTES = "classify_cassettes" 
    PLUG_CASSIOPEIA = "plug_cassiopeia"
    EXTRACT_BARCODES = "extract_barcodes"
    GENERATE_MATRIX = "generate_matrix"
    GENERATE_METADATA = "generate_metadata"

class CassiopeiaLineageExtension(PostProcessorExtension):
    xp: FractureXp
    config: ExtensionConfig
    def __init__(self, xp: FractureXp):
        super().__init__(xp)
        self.temp_data = {}
    
    def get_config_class(self) -> Type[ExtensionConfig]:
        return CassiopeiaConfig

    @property
    def required_params(self) -> Set[str]:
        return self.config.get_required_fields()
    
    @property
    def name(self) -> str:
        return "cassiopeia_lineage"
    
    def process(self, contigs_path: Path) -> StepResults:
        """Main entry point - orchestrates all sub-steps"""
        
        # Initialize
        #self.temp_data['contigs'] = pl.read_parquet(contigs_path)
        config = self.config
        force = getattr(config, 'force', False)  
        
        # Run steps based on configuration
        self.workdir = contigs_path.parent
        parsed_path = contigs_path.with_stem(contigs_path.stem + '_parsed')
        cass_mols_path = contigs_path.with_stem(contigs_path.stem + '_cass_mols')
        cass_allele_path = contigs_path.with_stem(contigs_path.stem + '_cass_allele')
        refs_path = contigs_path.parent / "refs.parquet"
        ann_intbc_mod_path = contigs_path.parent / "ann_intbc_mod_path.parquet"

        self.ldf = pl.scan_parquet(contigs_path).filter(pl.col('contig').str.len_chars() > 0)

        if self.config.refs_fasta_path is not None:
            self.refs = generate_refs_from_fasta(
                    **self.config.get_function_config(generate_refs_from_fasta)
            )
            self.refs.write_parquet(refs_path)
            save_ref_to_fasta(self.refs, out_dir=self.workdir, field='mod')

        #self.temp_data['xp'] = xp
        #self.temp_data['outputs_dir'] = contigs_path.parent / "cassiopeia_outputs"
        #self.temp_data['outputs_dir'].mkdir(exist_ok=True)
        #xp.logger.info(f"{self.temp_data['outputs_dir']=}")

        final_results = {}
        final_metrics = {}
        

        if not parsed_path.exists() or force:
            self.xp.logger.info("Running parse_contigs step")
            result = self._parse_contigs()
            self.ldf = result.results['ldf']
            self.ldf.sink_parquet(parsed_path)
            self.xp.logger.info(f"Writing parsed contigs to {parsed_path}")
            final_metrics.update(result.metrics)
        else:
            self.xp.logger.info(f"Loading parsed contigs from {parsed_path}")
            self.ldf = pl.scan_parquet(parsed_path)
        
        if self.should_run_step(CassiopeiaStep.CLASSIFY_CASSETTES.value):
            self.xp.logger.info("Running classify_cassettes step")

            #if annotation_path is None:
            if True:
                result = self._kmer_classify_cassettes()
            else:
                pass
                #result = use a given intbc_mod_map 
            self.xp.logger.io(f"Saving intBC mod annotations to {ann_intbc_mod_path}")
            # result.results['ann_intbc_mod'] # maps intBC to mod
            self.ann_intbc_mod = result.results['ann_intbc_mod']
            self.ann_intbc_mod.write_parquet(ann_intbc_mod_path) 

            final_results.update(result.results)
            final_metrics.update(result.metrics)

        if self.should_run_step(CassiopeiaStep.PLUG_CASSIOPEIA.value):
            self.xp.logger.info("Running plug_cassiopeia step")
            result = self._plug_cassiopeia()
            self.alleles_pl = result.results['alleles_pl']
            self.alleles_pl.write_parquet(f"{self.workdir}/alleles_pl.parquet")
            self.xp.logger.io(f"Saving cassiopeia allele tables to {cass_allele_path}")
            self.ldf.sink_parquet(cass_allele_path)

            
        if self.should_run_step(CassiopeiaStep.GENERATE_MATRIX.value):
            self.xp.logger.info("Running generate_matrix step")
            result = self._generate_matrix()
            final_results.update(result.results)
            final_metrics.update(result.metrics)
        
        if self.should_run_step(CassiopeiaStep.GENERATE_METADATA.value):
            self.xp.logger.info("Running generate_metadata step")
            result = self._generate_metadata()
            final_results.update(result.results)
            final_metrics.update(result.metrics)
        
        return StepResults(results=final_results, metrics=final_metrics)

    def _template(self) -> StepResults:
        """ """
        return StepResults(
            results={"":[]},
            metrics={"":[]},
        )

    def _parse_contigs(self) -> StepResults:
        return parse_contigs(ldf=self.ldf, 
                             **self.config.get_function_config(parse_contigs))
    
    def _kmer_classify_cassettes(self) -> StepResults:
        """ Guesses what reference a given intBC corresponds to based on kmer composition analysis"""
        return kmer_classify_cassettes(
                ldf=self.ldf, 
                refs=self.refs, 
                **self.config.get_function_config(kmer_classify_cassettes)
                )


    def _convert_to_allele_table(self) -> StepResults:
        """ """
        return StepResults(
            results={"":[]},
            metrics={"":[]},
        )
    
    def _plug_cassiopeia(self) -> StepResults:
        """ 
        readName - A unique identifier for each row/sequence
        cellBC - The cell barcode
        UMI - The UMI (Unique Molecular Identifier)
        readCount - The number of reads for this sequence
        seq - The actual sequence to be aligned

        """
        # Update config with runtime parameter
        self.config.ann_intbc_mod = self.ann_intbc_mod
        config = self.config.get_function_config(plug_cassiopeia)

        return plug_cassiopeia(ldf=self.ldf,
                               workdir=self.workdir,
                               logger=self.xp.logger,
                               **config)
    
    def _pycea_explore(self):
        
       """ """ 
       solvers = ['vanilla', 'mcgs', 'nj']
       solvers = ['vanilla']
       #for index, n_cells in [(i, ii.shape) for i,ii in enumerate(umi_tables) if ii.shape[0]>200]:
       for (index, intbc, mod, n_ori, well), allele_table in alg_df.partition_by('xid', 'intBC', 'mod', 'n_ori', 'well', as_dict=True).items():
           allele_table = allele_table.to_pandas()
           tdata = None
           character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(
               allele_table,
               allele_rep_thresh = 1.0)
           
           cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix)
           
           print(f"{index = } {intbc=}\tnumber of cells {cas_tree.n_cell}, number of characters {cas_tree.n_character} ")
           meta = '-'.join(pl.DataFrame(allele_table).select('xid', 'well', 'mod', 'intBC', 'n_ori').unique()[0].to_dummies(separator="_").columns)
           collapse = False
           collapse = True
           
           allele_colors_hex = dict(map(lambda x: (x, ColorHash(x).hex), pl.DataFrame(allele_table).unpivot(on=rcols, value_name="allele").get_column('allele').unique()))
           
           allele_matrix = pl.DataFrame(allele_table).select('UMI', '^r\d+$').to_pandas().set_index('UMI')
           allele_matrix
           for solver in solvers :
               match solver:
                   case "shared":
                       solve = cas.solver.SharedMutationJoiningSolver()
                   case "vanilla":
                       solve = cas.solver.VanillaGreedySolver()
                   case "UPGMA":
                       solve = cas.solver.SharedMutationJoiningSolver()
                   case "mcgs":
                       solve = cas.solver.MaxCutGreedySolver()
                   case "mc": # this one is very slow
                       solve = cas.solver.MaxCutSolver()
                   case "nj":
                       solve = cas.solver.NeighborJoiningSolver(add_root=True)

               print(solver)
               solve.solve(cas_tree, collapse_mutationless_edges=collapse)
               if tdata is None:
                   tdata = td.TreeData(
                       X=None,  
                       allow_overlap=True,
                       obs=allele_matrix.loc[cas_tree.leaves], 
                      # obsm={"alleles": character_matrix.loc[cas_tree.leaves].values}, # can't use this since must be encoded, e.g. character_matrix
                       obst={solver: cas_tree.get_tree_topology()}
                   )
               
               name = f'pctree_{cas_tree.n_cell}_{meta}_{solver}_{collapse}'

               tdata.obst[solver] = cas_tree.get_tree_topology()
               pycea.pp.add_depth(tdata, tree=solver)

               if True:
                   fig, ax = plt.subplots(1,1, figsize=(10, 100))
                   pc.pl.tree(tdata, 
                          tree=solver,
                     keys=rcols,
                     polar=False, 
                     extend_branches=True,
                     palette=allele_colors_hex,
                              branch_linewidth=0.15,
                     ax=ax)
                   fig.savefig(f'{name}.png')
                   
                   if False:
                       fig, ax = plt.subplots(1,1, figsize=(10, 10), dpi=1200, subplot_kw={"projection": "polar"})
                       pc.pl.tree(tdata, 
                              tree=solver,
                         keys=rcols,
                         polar=True, 
                         extend_branches=True,
                                  branch_linewidth=0.15,
                         palette=allele_colors_hex,
                         ax=ax)
                       fig.savefig(f'{name}_circ.png')

                   plt.close('all')
                
    def _extract_barcodes(self) -> StepResults:
        """Extract integration and sample barcodes"""
        xp = self.temp_data['xp']
        df_contigs = self.temp_data['contigs']
        
        # Your barcode extraction logic here
        df_annotated = (
            df_contigs
            .with_columns([
                pl.col('contig')
                .str.extract(f'({xp.intbc_5prime}[ATCG]{{10,20}})')
                .alias('integration_barcode'),
                
                pl.col('contig')
                .str.slice(0, xp.sbc_len)
                .alias('sample_barcode')
            ])
            .filter(pl.col('integration_barcode').is_not_null())
        )
        
        output_path = self.temp_data['outputs_dir'] / "barcodes_extracted.parquet"

        df_annotated.write_parquet(output_path)
        self.temp_data['annotated_contigs'] = df_annotated
        
        return StepResults(
            results={"barcodes_extracted": str(output_path)},
            metrics={
                "contigs_with_integration_bc": df_annotated.height,
                "unique_integration_bcs": df_annotated.select('integration_barcode').n_unique()
            }
        )
    
    
    def _generate_matrix(self) -> StepResults:
        # Implementation here  
        pass
        
    def _generate_metadata(self) -> StepResults:
        # Implementation here
        pass


def generate_mask_PEtracer_expression(features_csv,
                                      column_name="seq",
                                      fuzzy_pattern:bool = True,
                                      fuzzy_kwargs:dict = None,
                                      logger: None|CustomLogger = None
                                      ) -> pl.Expr:
    """
    Generate a Polars expression to replace TARGET sequences with scrambled versions
    based on preceding META sequences.

    Example of usage:
        expr = generate_mask_PEtracer_expression(
            'features.csv',
            fuzzy_pattern=True,
            fuzzy_kwargs={
                'wildcard': '.{0,2}',      # Allow up to 2 character variations
                'include_original': False,  # Don't include exact match
                'sep': '|',                # Use pipe separator
                'max_length': 150          # Only fuzzify sequences up to 150 chars
            }
        )
    
    Expected csv format:
        feature,seq,kind
        META01,AGAAGCCGTGTGCCGGTCTA,META
        META02,ATCGTGCGGACGAGACAGCA,META
        RNF2,TGGCAGTCATCTTAGTCATTACGACAGGTGTTCGTTGTAACTCATATA,TARGET
        HEK3,CTTGGGGCCCAGACTGAGCACGACTTGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCA,TARGET
        EMX1,GGCCTGAGTCCGAGCAGAAGAACTTGGGCTCCCATCACATCAACCGGTGG,TARGET
    
    Args:
        features_csv: Path to CSV file containing feature definitions
        column_name: Name of the column to apply replacements to (default: "seq")
        fuzzy_pattern: Whether to perform fuzzy string matching to account for sequencing errors
        fuzzy_kwargs: Dictionary of keyword arguments to pass to fuzzy_match_str function
                     Supported keys: wildcard, include_original, sep, max_length
                     
    Returns:
        pl.Expr: Polars expression with chained replacements
       
    Usage:
        # Basic usage
        expr = generate_mask_PEtracer_expression('features.csv')
        
        # With custom fuzzy matching parameters
        expr = generate_mask_PEtracer_expression(
            'features.csv',
            fuzzy_pattern=True,
            fuzzy_kwargs={
                'wildcard': '.{0,2}',  # Allow up to 2 character variations
                'include_original': False,  # Don't include exact match
                'max_length': 100  # Only fuzzify sequences up to 100 chars
            }
        )
        
        result = (
            df.lazy()
            .with_columns(expr.alias("result"))
            .collect()
        )
    """
    import polars as pl
    import random
    
    if fuzzy_kwargs is None:
        fuzzy_kwargs = {}
    
    def fuzzy_match_str(string, wildcard=".{0,1}", include_original=True, sep='|', max_length=200):
        """Generate fuzzy regex pattern allowing single character variations.
        
        Args:
            string: Input string to fuzzify
            wildcard: Regex pattern for character substitution
            include_original: Whether to include exact match
            sep: Separator for alternatives
            max_length: Skip fuzzy variants for strings longer than this
            
        Returns:
            Regex pattern string for fuzzy matching
        """
        if not string:
            return string
            
        fuzz = []
        if include_original:
            fuzz.append(string)
            
        # Skip fuzzy generation for very long strings to avoid performance issues
        if len(string) <= max_length:
            for i in range(len(string)):
                # Use list comprehension + join instead of character-by-character building
                variant = ''.join(wildcard if j == i else c for j, c in enumerate(string))
                fuzz.append(variant)
                
        return sep.join(fuzz)
    
    def scramble_dna_sequence(seed_string, sequence, fuzzy_pattern):
        """
        Scramble a DNA sequence using a deterministic seed.
        """
        seed = hash(seed_string) % (2**32)  # Ensure positive 32-bit integer
        random.seed(seed)
        seq_list = list(sequence.upper())
        random.shuffle(seq_list)
        scrambled_seq = ''.join(seq_list)
        if fuzzy_pattern:
            return fuzzy_match_str(scrambled_seq, **fuzzy_kwargs)
        return scrambled_seq 
    
    meta = pl.read_csv(features_csv)
    
    patterns_df = (
        meta
        .filter(pl.col('kind') == 'META')
        .join(
            meta.filter(pl.col('kind') == 'TARGET'),
            suffix="_target",
            how='cross'
        )
        .with_columns(
            # Generate scrambled sequence for each META-TARGET combination
            # Original version (META-specific): 
            # lambda x: scramble_dna_sequence(x['feature'], x['seq_target'])
            pl.struct(['feature', 'feature_target']).map_elements(
                lambda x: scramble_dna_sequence(
                    seed_string = f"{x['feature']}_{x['feature_target']}", 
                    sequence  = meta.filter(pl.col('feature') == x['feature_target']).get_column('seq')[0],
                    fuzzy_pattern = fuzzy_pattern),
                return_dtype=pl.Utf8
                ).alias("scrambled_seq")
            )
        .with_columns(
            # pattern = rf"("+pl.col('seq')+")(.*?)("+pl.col('seq_target')+")(.*$)",
            # replacement = rf"$1${{2}}"+pl.col('scrambled_seq')+"$4"
            pattern=pl.concat_str([
                pl.lit("("),
                pl.col('seq'),  # META sequence
                pl.lit(")(.*?)("),
                pl.col('seq_target'),  # TARGET sequence
                pl.lit(")(.*$)")
            ]),
            replacement=pl.concat_str([
                pl.lit("$1${2}"),
                pl.col('scrambled_seq'),
                pl.lit("$4")
            ])
        )
    )
    
    if logger is not None:
        logger.debug(f"Generated {len(patterns_df)} replacement patterns")
        if len(patterns_df) > 0:
            logger.debug("Sample patterns:")
            sample = patterns_df.head(3).select(['feature', 'feature_target', 'pattern', 'replacement'])
            for row in sample.iter_rows(named=True):
                logger.debug(f"  {row['feature']} + {row['feature_target']}: {row['pattern'][:60]}...")
                logger.debug(f"    -> {row['replacement'][:60]}...")
    
    # Build the chained expression
    expr = pl.col(column_name)
    
    for row in patterns_df.iter_rows(named=True):
        pattern = row['pattern']
        replacement = row['replacement']
        expr = expr.str.replace(pattern, replacement)
    
    print(f"Built expression with {len(patterns_df)} chained replacements")
    return expr
# Register the extension
extension_registry.register(CassiopeiaLineageExtension)

