import polars as pl
from pathlib import Path
from typing import Set,Dict, Optional, Type
from enum import Enum
from .base import PostProcessorExtension
from .registry import extension_registry
from .config import ExtensionConfig
from ..pipeline.types import StepResults,FractureXp
from dataclasses import dataclass

def save_ref_to_fasta(refs, out_dir='.', field='mod') -> None : 
    """ write fastas in current dir"""
    for i in refs.get_column(field):
        with open(f'{out_dir}/{i}.fasta', 'wt') as out:
            out.write(
            refs.filter(pl.col(field)==i).dna.to_fasta(read_id_col=field, read_col='seq').get_column('seq_fasta')[0]
            )

def kmer_classify_cassettes(ldf, refs, **config) -> StepResults:
    """ Annotates contigs based on the top occurring mod based on kmer matches """
    K = config.get('K', 25)

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
    return StepResults(results={'ann_intbc_mod':cont_k}, metrics={'n_ann_intbc':cont_k.height})

def parse_contigs(ldf, **config):
    """Extracts integration barcodes (intBC) and sample barcodes (sbc)."""
    int_anchor1 = config['int_anchor1']
    int_anchor2 = config['int_anchor2'] 
    sbc_dict = config.get('sbc_dict')
    annotation = config.get('annotation', 'sbc')
    
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
    
    return StepResults(results={'ldf':ldf}, metrics={'n_parsed_contigs':df.height})


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
    
    # Function-specific parameter mapping
    _FUNCTION_PARAMS = {
        'parse_contigs': {'int_anchor1', 'int_anchor2', 'sbc_dict', 'annotation'},
        'generate_refs_from_fasta': {'refs_fasta_path', 'anchor1', 'anchor2'},
        'classify_cassettes': {'int_anchor1', 'refs_fasta_path', 'annotation_path'},
        'plug_cassiopeia': {},
    }

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
                self.config.refs_fasta_path, 
                anchor1=self.config.anchor1, 
                anchor2=self.config.anchor2
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
        return parse_contigs(df=self.ldf, **self.config.get_function_config('parse_contigs'))
    
    def _kmer_classify_cassettes(self) -> StepResults:
        """ Guesses what reference a given intBC corresponds to based on kmer composition analysis"""
        return kmer_classify_cassettes(ldf=self.ldf, 
                                       refs=self.refs, 
                                       **self.config.get_function_config('classify_cassettes'))


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
        import cassiopeia as cas
        ann_intbc_mod = self.ann_intbc_mod

        cass_ldf = (
                self.ldf.with_columns(
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

        for (intBC, mod), queries in cass_ldf.collect().partition_by('intBC', 'mod', as_dict=True).items():
            if mod is None:
                nones+=1
            else:
                self.xp.logger.debug(f"{mod=} {intBC=} {queries.shape=}")

                umi_table = cas.pp.align_sequences(
                     queries=queries.to_pandas(),
                     ref_filepath=f'{self.workdir}/{mod}.fasta',
                     n_threads=1,
                     method='global'
                 )

                allele_table =  cas.pp.call_alleles(
                                umi_table,
                                ref_filepath=f'{self.workdir}/{mod}.fasta',
                                barcode_interval=(0, 7),
                                cutsite_locations=[40, 67, 94, 121, 148, 175, 202, 229, 256, 283],
                                cutsite_width=12,
                                context=True,
                                context_size=5,
                            )

                #replace intBC for real intBC since Cassiopeia does something I don't understand yet
                allele_table['intBC'] = intBC
                allele_table['mod'] = mod

                umi_tables.append(umi_table.copy())
                #self.xp.logger.info(f"found {umi_table['intBC'].n_unique()} integrations for {mod}")

                res.append(
                    pl.DataFrame(allele_table).with_columns(mod=pl.lit(mod), mols=allele_table.shape[0])
                    )

        return StepResults(
                results={"alleles_pl" : pl.concat(res), "alleles_pd" : umi_tables}
        )

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

# Register the extension
extension_registry.register(CassiopeiaLineageExtension)
