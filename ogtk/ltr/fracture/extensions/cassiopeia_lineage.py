import polars as pl
from pathlib import Path
from typing import Set,Dict, Optional, Type
from enum import Enum
from .base import PostProcessorExtension
from .registry import extension_registry
from .config import ExtensionConfig
from ..pipeline.types import StepResults,FractureXp
from dataclasses import dataclass

def parse_contigs(df, **config):
    """Extracts integration barcodes (intBC) and sample barcodes (sbc)."""
    int_anchor1 = config['int_anchor1']
    int_anchor2 = config['int_anchor2'] 
    sbc_dict = config.get('sbc_dict')
    annotation = config.get('annotation', 'sbc')
    
    # 1. Extract/use sbcs
    if sbc_dict is not None:
        df = df.with_columns(pl.col('sbc').replace(sbc_dict).alias(annotation))
    
    # 2. Extract intBC
    df = (
        df   
        .with_columns(
            pl.col('contig').str
            .extract(f'{int_anchor1}(.+?){int_anchor2}', 1).alias('intBC')
        )
    )
    
    return df


def _generate_refs_from_fasta(refs_fasta_path: str|Path, anchor1: str, anchor2: str) -> pl.DataFrame:
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
        'classify_cassettes': {'int_anchor1', 'refs_fasta_path'},
    }

class CassiopeiaStep(Enum):
    """Steps within the Cassiopeia lineage extension"""
    PARSE_CONTIGS = 'parse_contigs'
    CLASSIFY_CASSETTES = "classify_cassettes" 
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
        self.df = pl.read_parquet(contigs_path)

        if self.config.refs_fasta_path is not None:
            refs = _generate_refs_from_fasta(
                self.config.refs_fasta_path, 
                anchor1=self.config.anchor1, 
                anchor2=self.config.anchor2
            )

        #self.temp_data['xp'] = xp
        #self.temp_data['outputs_dir'] = contigs_path.parent / "cassiopeia_outputs"
        #self.temp_data['outputs_dir'].mkdir(exist_ok=True)
        #xp.logger.info(f"{self.temp_data['outputs_dir']=}")

        final_results = {}
        final_metrics = {}
        
        # Run steps based on configuration

        if self.should_run_step(CassiopeiaStep.PARSE_CONTIGS.value):
            self.xp.logger.info("Running parse_contigs step")
            result = self._parse_contigs()
            final_results.update(result.results)
            final_metrics.update(result.metrics)

        if self.should_run_step(CassiopeiaStep.EXTRACT_BARCODES.value):
            self.xp.logger.info("Running extract_barcodes step")
            result = self._extract_barcodes()
            final_results.update(result.results)
            final_metrics.update(result.metrics)
        
        if self.should_run_step(CassiopeiaStep.CLASSIFY_CASSETTES.value):
            self.xp.logger.info("Running classify_cassettes step")
            result = self._classify_cassettes()
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
        return parse_contigs(df=self.df, **self.config.get_function_config('parse_contigs'))

    def _convert_to_allele_table(self) -> StepResults:
        """ """
        return StepResults(
            results={"":[]},
            metrics={"":[]},
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
    
    # Implement other step methods...
    def _classify_cassettes(self) -> StepResults:
        """ Guesses what reference a given intBC corresponds to based on kmer compositin analysis"""

        
        return StepResults(
                results={
                'parsed_fn': "Tururu",
                },
                metrics = {}
            )
    
    def _generate_matrix(self) -> StepResults:
        # Implementation here  
        pass
        
    def _generate_metadata(self) -> StepResults:
        # Implementation here
        pass

# Register the extension
extension_registry.register(CassiopeiaLineageExtension)
