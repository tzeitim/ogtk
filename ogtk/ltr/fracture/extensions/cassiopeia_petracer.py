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
        gap_open_penalty: float = 20.0,
        gap_extend_penalty: float = 1.0,
        alignment_method: str = 'global',
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
                 method=alignment_method,
                 gap_open_penalty=gap_open_penalty,
                 gap_extend_penalty=gap_extend_penalty,
             )
    
            allele_table =  cas.pp.call_alleles(
                            umi_table,
                            ref_filepath = f'{workdir}/{mod}.fasta',
                            **allele_params,
                        )

            # Enrich allele columns with actual insertion sequences from CIGAR
            # allele_table already contains Seq and CIGAR columns from umi_table
            pl_allele = pl.DataFrame(allele_table)
            rcols = pl_allele.select(pl.col('^r\\d+$')).columns
            if rcols and 'Seq' in pl_allele.columns and 'CIGAR' in pl_allele.columns:
                pl_allele = pl_allele.with_columns([
                    pl.col(col).cigar.enrich_insertions(pl.col('Seq'), pl.col('CIGAR'))
                    for col in rcols
                ])
            allele_table = pl_allele.to_pandas()

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

def save_ref_to_fasta(refs: pl.DataFrame, out_dir: str|Path = '.', field: str = 'mod') -> None:
    """Write reference sequences to individual FASTA files.

    Args:
        refs: DataFrame with 'mod' and 'seq' columns
        out_dir: Output directory for FASTA files
        field: Column to use for file naming and filtering
    """
    if refs.height == 0:
        raise ValueError("refs DataFrame is empty - no references to save")

    unique_values = refs.get_column(field).unique()
    if unique_values.null_count() > 0:
        raise ValueError(f"refs DataFrame has null values in '{field}' column")

    for i in unique_values:
        filtered = refs.filter(pl.col(field) == i)
        if filtered.height == 0:
            raise ValueError(f"No rows found for {field}={i}")

        fasta_content = filtered.dna.to_fasta(read_id_col=field, read_col='seq').get_column('seq_fasta')[0]
        out_path = Path(out_dir) / f"{i}.fasta"
        with open(out_path, 'wt') as out:
            out.write(fasta_content)

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
    
    # 2. Extract intBC and trim contig to remove the intBC region
    ldf = (
        ldf
        .with_columns(
            pl.col('contig').str
            .extract(f'{int_anchor1}(.+?){int_anchor2}', 1).alias('intBC'),
            pl.col('contig').str.replace(f'.*{int_anchor2}', '').alias('contig')
        )
    )
    
    # Get count lazily - will be computed when needed
    n_contigs = ldf.select(pl.len()).collect().item()
    return StepResults(results={'ldf':ldf},
                       metrics={'n_parsed_contigs': n_contigs})


def generate_refs_from_fasta(refs_fasta_path: str|Path, anchor1: str, anchor2: str) -> pl.DataFrame:
    """ Given a FASTA file with references it generates a trimmed version of the cassettes for aligment

    Supports FASTA headers in multiple formats:
    - ">5mer" → mod_5mer
    - ">v5mer_xxx" → mod_5mer
    - ">10mer" → mod_10mer
    """

    refs =  (
            pl.read_csv(refs_fasta_path,
                        has_header=False,
                        )
            .unstack(step=2, how='horizontal')
            .rename({"column_1_0":"mod", "column_1_1":"seq"})
            .with_columns(
                # Handle multiple header formats:
                # ">5mer" -> "5mer" -> "mod_5mer"
                # ">v5mer_xxx" -> "5mer" -> "mod_5mer"
                pl.concat_str([
                    pl.lit("mod_"),
                    pl.col('mod')
                        .str.strip_chars(">")  # Remove leading >
                        .str.extract(r"^v?(\d+mer)", 1)  # Extract Nmer pattern (optional v prefix)
                ]).alias('mod'),
                pl.col('seq').str.to_uppercase())
            # trim sequence
            .with_columns(
                pl.col('seq').str.replace(f".+?({anchor1})", anchor1).str.replace(f"({anchor2}).+?$", anchor2)
                )
            )

    # Validate we got valid mod values
    null_mods = refs.filter(pl.col('mod').is_null()).height
    if null_mods > 0:
        raise ValueError(f"Failed to parse {null_mods} FASTA headers. Expected format: >5mer, >10mer, or >v5mer_xxx")

    return refs


# ============================================================================
# Segment-based allele table generation
# ============================================================================

CASSETTE_START_MARKER = "_CASSETTE_START_"
CASSETTE_END_MARKER = "_CASSETTE_END_"
TARGET_ORDER = ["RNF2", "HEK3", "EMX1"]  # Order within each triplet

CASSETTE_CONFIGS = {
    "5mer": {"n_metas": 4, "n_targets": 12},
    "10mer": {"n_metas": 9, "n_targets": 27},
    "20mer": {"n_metas": 19, "n_targets": 57},
}


def filter_segment_types_for_cassette(
    segment_types: pl.DataFrame,
    metas: pl.DataFrame,
    cassette_type: str,
    min_freq_threshold: float = 0.7,
    total_molecules: int | None = None,
    logger: Optional[CustomLogger] = None,
) -> pl.DataFrame:
    """
    Filter segment_types to only include expected metas for the cassette type.

    For a 5mer cassette, only META01-META04 are expected. Any other metas
    (META05+) are filtered out UNLESS they appear at high frequency (>min_freq_threshold),
    which may indicate they are real (e.g., a different cassette version).

    Args:
        segment_types: DataFrame with start_meta, end_meta, count columns
        metas: DataFrame with feature column (META01, META02, etc.)
        cassette_type: One of "5mer", "10mer", "20mer"
        min_freq_threshold: Minimum frequency (0-1) for an unexpected meta to be kept.
                           Default 0.7 (70%). Set to 1.0 to strictly filter.
        total_molecules: Total number of molecules for frequency calculation.
                        If None, uses sum of counts in segment_types.
        logger: Optional logger for debugging

    Returns:
        Filtered segment_types DataFrame
    """
    if cassette_type not in CASSETTE_CONFIGS:
        if logger:
            logger.warning(f"Unknown cassette type '{cassette_type}', skipping meta filtering")
        return segment_types

    n_metas = CASSETTE_CONFIGS[cassette_type]["n_metas"]

    # Get ordered meta names from metas DataFrame
    meta_only = metas.filter(pl.col('kind') == 'META') if 'kind' in metas.columns else metas
    all_meta_names = meta_only['feature'].to_list()

    # Expected metas for this cassette type (first n_metas)
    expected_metas = set(all_meta_names[:n_metas])
    # Add edge markers as always valid
    expected_metas.add(CASSETTE_START_MARKER)
    expected_metas.add(CASSETTE_END_MARKER)

    if logger:
        logger.debug(f"Cassette {cassette_type}: expecting {n_metas} metas: {sorted(expected_metas)}")

    # Calculate frequency of each meta across all transitions
    if total_molecules is None:
        total_molecules = segment_types['count'].sum()

    if total_molecules == 0:
        return segment_types

    # Find all unique metas in segment_types
    all_starts = set(segment_types['start_meta'].unique().to_list())
    all_ends = set(segment_types['end_meta'].unique().to_list())
    all_observed_metas = all_starts | all_ends

    # Check for unexpected metas with high frequency
    unexpected_metas = all_observed_metas - expected_metas
    high_freq_unexpected = set()

    for meta in unexpected_metas:
        # Count molecules with this meta (in either start or end)
        meta_count = segment_types.filter(
            (pl.col('start_meta') == meta) | (pl.col('end_meta') == meta)
        )['count'].sum()
        freq = meta_count / total_molecules

        if freq >= min_freq_threshold:
            high_freq_unexpected.add(meta)
            if logger:
                logger.info(f"Keeping unexpected meta '{meta}' - appears in {freq:.1%} of molecules (>= {min_freq_threshold:.0%} threshold)")
        elif logger:
            logger.debug(f"Filtering out meta '{meta}' - appears in {freq:.1%} of molecules (< {min_freq_threshold:.0%} threshold)")

    # Final whitelist: expected + high-frequency unexpected
    whitelist = expected_metas | high_freq_unexpected

    # Filter segment_types to only include transitions with whitelisted metas
    filtered = segment_types.filter(
        pl.col('start_meta').is_in(list(whitelist)) &
        pl.col('end_meta').is_in(list(whitelist))
    )

    n_removed = segment_types.height - filtered.height
    if n_removed > 0 and logger:
        logger.info(f"Filtered {n_removed} transitions with unexpected metas (kept {filtered.height})")

    return filtered


def build_target_position_map(
    metas_df: pl.DataFrame,
    cassette_type: str = "5mer",
) -> Dict[Tuple[str, str], List[Tuple[str, int]]]:
    """
    Build mapping from (start_meta, end_meta) segment boundaries to TARGET positions.

    Logic:
    - For 5mer: 4 METAs (META01-META04), creating 4 segments each with 3 TARGETs = 12 positions
    - Segment (_CASSETTE_START_, META01) → positions r1, r2, r3 (RNF2, HEK3, EMX1)
    - Segment (META01, META02) → positions r4, r5, r6
    - etc.

    Returns:
        Dict mapping (start_meta, end_meta) -> [(target_name, position_index), ...]
    """
    raise NotImplementedError("TODO: implement with polars")


def segments_to_allele_table(
    segments_df: pl.DataFrame,
    metas_df: pl.DataFrame,
    int_anchor1: str,
    int_anchor2: str,
    cassette_type: str = "5mer",
    min_consensus_support: float = 0.5,
    logger: Optional[CustomLogger] = None,
) -> pl.DataFrame:
    """
    Transform segmentation table to cassiopeia allele table format.

    Input columns: [sbc, umi, start_meta, end_meta, segment_seq]

    Output columns:
    - cellBC, UMI, readCount, intBC, sbc, mod
    - r1, r2, r3, ... rN (one per TARGET position)

    Logic:
    1. Extract intBC using: pl.col('segment_seq').str.extract(f'{int_anchor1}(.+?){int_anchor2}', 1)
    2. Map (start_meta, end_meta) to TARGET positions using build_target_position_map()
    3. For each TARGET in segment, extract insertion using flanks from metas_df
       - left_flank and right_flank columns identify TARGET boundaries
       - Insertion = sequence between left_flank and right_flank
       - Empty string = wild-type, None = flanks not found
    4. Group by (sbc, umi) and consensus across multiple observations
    5. Pivot to wide format: r1, r2, r3, ...

    Values in rN columns:
    - "" (empty): Wild-type, no insertion
    - "ACTGT": Insertion sequence (lineage mark)
    - None: Missing data
    """
    raise NotImplementedError("TODO: implement with polars")


@dataclass
class CassiopeiaConfig(ExtensionConfig):
    # Required fields
    int_anchor1: str
    int_anchor2: str
    force: bool = False

    # Optional fields with defaults
    sbc_dict: Optional[Dict] = None
    annotation: str = 'sbc'
    refs_fasta_path: Optional[str] = None
    anchor1: Optional[str] = None
    anchor2: Optional[str] = None

    barcode_interval: Tuple[int, int] = (0, 7)
    cutsite_locations: List[int] = field(default_factory=lambda: [40, 67, 94, 121, 148, 175, 202, 229, 256, 283])
    cutsite_width: int = 12
    context: bool = True
    context_size: int = 50

    # Alignment parameters for cas.pp.align_sequences
    gap_open_penalty: float = 20.0  # Penalty for opening a gap (default: 20)
    gap_extend_penalty: float = 1.0  # Penalty for extending a gap (default: 1)
    alignment_method: str = 'global'  # 'local' (Smith-Waterman) or 'global' (Needleman-Wunsch)

    # Segmented allele extraction fields
    metas_flanks_csv: Optional[str] = None  # Path to PEtracer_metas_flanks.csv
    cassette_type: str = "5mer"  # Options: 5mer, 10mer, 20mer
    min_consensus_support: float = 0.5  # Minimum support for consensus calls

    # Runtime fields (populated during processing, not from config file)
    ann_intbc_mod: Optional[pl.DataFrame] = None  # intBC → mod mapping from classify_cassettes
    

class CassiopeiaStep(Enum):
    """Steps within the Cassiopeia lineage extension"""
    PARSE_CONTIGS = 'parse_contigs'
    CLASSIFY_CASSETTES = "classify_cassettes"
    REGENERATE_FILTERED_QC = "regenerate_filtered_qc"  # Regenerate QC with cassette-type filtering
    PLUG_CASSIOPEIA = "plug_cassiopeia"
    EXTRACT_BARCODES = "extract_barcodes"
    GENERATE_MATRIX = "generate_matrix"
    GENERATE_METADATA = "generate_metadata"
    SEGMENTED_ALLELE = "segmented_allele"  # Direct allele table from segments

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
        return "cassiopeia_petracer"
    
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

        if self.should_run_step(CassiopeiaStep.REGENERATE_FILTERED_QC.value):
            self.xp.logger.info("Running regenerate_filtered_qc step")
            result = self._regenerate_filtered_qc()
            final_results.update(result.results)
            final_metrics.update(result.metrics)

        if self.should_run_step(CassiopeiaStep.PLUG_CASSIOPEIA.value):
            self.xp.logger.info("Running plug_cassiopeia step")
            result = self._plug_cassiopeia()
            self.alleles_pl = result.results['alleles_pl']
            self.alleles_pl.write_parquet(f"{self.workdir}/alleles_pl.parquet")
            self.xp.logger.io(f"Saving cassiopeia allele tables to {cass_allele_path}")
            self.ldf.sink_parquet(cass_allele_path)

        if self.should_run_step(CassiopeiaStep.SEGMENTED_ALLELE.value):
            self.xp.logger.info("Running segmented_allele step")
            result = self._segmented_allele()
            self.alleles_segmented = result.results['allele_table_segmented']
            final_results.update(result.results)
            final_metrics.update(result.metrics)

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

    def _regenerate_filtered_qc(self) -> StepResults:
        """
        Regenerate segmentation QC plots with cassette-type filtering.

        After cassette classification, we know the cassette type (5mer, 10mer, etc.)
        and can filter out spurious meta matches in the QC plots.

        This step:
        1. Loads segmentation data from intermediate files
        2. Determines cassette_type from refs (e.g., mod_5mer -> "5mer")
        3. Regenerates QC plots with filtered segment_types
        """
        from ..pipeline.segmentation import (
            generate_segmentation_report,
            plot_segmentation_qc,
        )

        # Determine cassette_type from refs
        # refs has 'mod' column like "mod_5mer", "mod_10mer"
        if not hasattr(self, 'refs') or self.refs is None:
            self.xp.logger.warning("No refs available, skipping filtered QC regeneration")
            return StepResults(results={}, metrics={})

        # Get the dominant cassette type from refs (or use config)
        mods = self.refs['mod'].unique().to_list()
        # Extract cassette type: "mod_5mer" -> "5mer"
        cassette_types = [m.replace("mod_", "") for m in mods if m.startswith("mod_")]

        if not cassette_types:
            self.xp.logger.warning("Could not determine cassette type from refs")
            return StepResults(results={}, metrics={})

        # Use the first (or only) cassette type
        cassette_type = cassette_types[0]
        self.xp.logger.info(f"Regenerating QC with cassette_type={cassette_type} filtering")

        # Load intermediate files
        intermediate_dir = self.workdir / "intermediate"
        segments_path = intermediate_dir / "segments.parquet"
        assembled_path = intermediate_dir / "assembled.parquet"

        # Also try legacy path
        if not segments_path.exists():
            segments_path = self.workdir / "segments_debug.parquet"

        if not segments_path.exists():
            self.xp.logger.warning(f"Segments file not found at {segments_path}, skipping filtered QC")
            return StepResults(results={}, metrics={})

        # Load metas
        metas_csv = self.xp.fracture.get('metas_csv') or getattr(self.xp, 'features_csv', None)
        if not metas_csv:
            self.xp.logger.warning("No metas_csv configured, skipping filtered QC")
            return StepResults(results={}, metrics={})

        metas = pl.read_csv(metas_csv)
        segments_df = pl.read_parquet(segments_path)

        # Load assembled if available
        if assembled_path.exists():
            assembled_df = pl.read_parquet(assembled_path)
        else:
            # Create minimal assembled_df from segments
            assembled_df = segments_df.select(['umi', 'start_meta', 'end_meta']).unique()

        # Load contigs
        contigs_path = list(self.workdir.glob("contigs_segmented_*.parquet"))
        if contigs_path:
            contigs_df = pl.read_parquet(contigs_path[0])
            if 'stitched_seq' in contigs_df.columns and 'contig' not in contigs_df.columns:
                contigs_df = contigs_df.rename({'stitched_seq': 'contig'})
        else:
            # Minimal contigs from ldf
            contigs_df = self.ldf.select(['umi', 'contig']).collect()

        # Get anchors from config
        cassette_start_anchor = getattr(self.xp, 'start_anchor', None)
        cassette_end_anchor = getattr(self.xp, 'end_anchor', None)

        # Get filtering threshold from config (default 70%)
        min_meta_freq = self.xp.fracture.get('min_meta_freq_threshold', 0.7)

        # Get heterogeneity threshold from config (default 0.20 = 20% of dominant)
        heterogeneity_threshold = self.xp.fracture.get('heterogeneity_threshold', 0.20)

        # Generate unfiltered report first
        report = generate_segmentation_report(
            segments_df=segments_df,
            assembled_df=assembled_df,
            contigs_df=contigs_df,
            metas=metas,
            cassette_start_anchor=cassette_start_anchor,
            cassette_end_anchor=cassette_end_anchor,
            heterogeneity_threshold=heterogeneity_threshold,
        )

        # Now filter the segment_types in the report
        segment_types_df = pl.DataFrame(report['segments']['segment_types'])
        unique_umis = report['segments']['unique_umis']

        filtered_segment_types = filter_segment_types_for_cassette(
            segment_types=segment_types_df,
            metas=metas,
            cassette_type=cassette_type,
            min_freq_threshold=min_meta_freq,
            total_molecules=unique_umis,
            logger=self.xp.logger,
        )

        # Update the report with filtered segment_types
        report['segments']['segment_types'] = filtered_segment_types.to_dicts()

        # Generate filtered QC plots - use central figures directory per Xp
        figures_dir = Path(self.xp.sample_figs)
        figures_dir.mkdir(parents=True, exist_ok=True)

        sample_name = getattr(self.xp, 'target_sample', 'sample')
        plot_segmentation_qc(
            report=report,
            contigs_df=contigs_df,
            metas=metas,
            output_dir=figures_dir,
            sample_name=f"{sample_name}_filtered",
            contig_col='contig',
            logger=self.xp.logger,
        )

        self.xp.logger.info(f"Saved filtered QC plots to {figures_dir}")

        return StepResults(
            results={"filtered_qc_dir": str(figures_dir)},
            metrics={"cassette_type": cassette_type},
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

    def _segmented_allele(self) -> StepResults:
        """
        Generate allele table directly from segments (no assembly needed).

        Reads segments_debug.parquet from workdir and generates allele_table_segmented.parquet.
        """
        segments_path = self.workdir / "segments_debug.parquet"

        if not segments_path.exists():
            raise FileNotFoundError(
                f"Segments file not found: {segments_path}. "
                "Run fracture with use_segmentation=True first."
            )

        segments_df = pl.read_parquet(segments_path)

        if self.config.metas_flanks_csv is None:
            raise ValueError("metas_flanks_csv must be set in config for segmented_allele step")

        metas_df = pl.read_csv(self.config.metas_flanks_csv)

        allele_table = segments_to_allele_table(
            segments_df=segments_df,
            metas_df=metas_df,
            int_anchor1=self.config.int_anchor1,
            int_anchor2=self.config.int_anchor2,
            cassette_type=self.config.cassette_type,
            min_consensus_support=self.config.min_consensus_support,
            logger=self.xp.logger,
        )

        output_path = self.workdir / "allele_table_segmented.parquet"
        allele_table.write_parquet(output_path)
        self.xp.logger.io(f"Saved segmented allele table to {output_path}")

        return StepResults(
            results={'allele_table_segmented': allele_table},
            metrics={
                'total_molecules': allele_table.height,
                'unique_intbcs': allele_table['intBC'].n_unique() if allele_table.height > 0 else 0,
            }
        )

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

# Register the extension
extension_registry.register(CassiopeiaLineageExtension)
