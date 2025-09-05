from ogtk.utils.db import Xp
from pathlib import Path
from typing import Any, List
from ogtk.utils.log import CustomLogger, Rlogger, call
from .plotting import PlotDB
from typing import NamedTuple,Dict

class FractureXp(Xp):
    """Extended Xp class with pipeline-specific functionality"""
    steps: Any
    dry: bool
    make_test: bool
    logger: CustomLogger
    pro_datain: str
    modality: str
    umi_len: int
    sbc_len: int
    cbc_len: int
    rev_comp: bool
    anchor_ont: str
    samples: List[str]
    pro_workdir: str
    plotdb: PlotDB
    do_plot: bool    
    fracture: dict
    start_anchor: str
    end_anchor: str
    force_tab: bool
    allow_wildcards: bool
    intbc_5prime: str
    parse_read1: bool
    extensions: List[str]
    extension_steps: Dict[str, List[str]]
    extension_config: Dict[str, Dict[str, Any]]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.steps = getattr(self, 'steps', None)
        self.dry = getattr(self, 'dry', False)
        self.make_test = getattr(self, 'make_test', False)
        self.parse_read1 = getattr(self, 'parse_read1', False)
        self.sbc_len = getattr(self, 'sbc_len', 6)
        
        # Set defaults based on modality
        modality = getattr(self, 'modality', 'single-molecule')
        if modality == 'single-cell':
            # For single-cell data: no sample barcode, set cell barcode length
            self.sbc_len = getattr(self, 'sbc_len', 0)  # Override default for single-cell
            self.cbc_len = getattr(self, 'cbc_len', 16)  # Default cell barcode length
        else:
            # For single-molecule data: standard settings
            self.cbc_len = getattr(self, 'cbc_len', 0)  # No cell barcode for single-molecule

        if not hasattr(self, 'anchor_ont'):
            raise ValueError("Missing required parameter 'anchor_ont' for BAM processing")

        # extensions
        self.extensions = getattr(self, 'extensions', [])
        self.extension_steps = getattr(self, 'extension_steps', {})
        self.extension_config = getattr(self, 'extension_config', {})


    @call
    def organize_files_by_sample(self, files, samples, max_files=None):
        # Create a dictionary to store files for each sample
        sample_files = {sample['id']: [] for sample in samples}
        
        # Sort files into appropriate sample groups
        empty_keys = []
        for f in files:
            for sample_id in sample_files:
                if Path(f).name.startswith(sample_id):
                    sample_files[sample_id].append(f)
                    break
                else:
                    empty_keys.append(sample_id)

        if max_files is not None and isinstance(max_files, int):
            if any([len(i)>max_files for i in sample_files.values()]):
                self.logger.error(f"Some sample(s) matched more than {max_files}\n{sample_files}", with_traceback=True)
                raise ValueError("Verify that there are files that coorespond to all samples.")
            
        return sample_files

class StepResults(NamedTuple):
    #TODO is it necessary so many levels? (results, metrics)
    """Container for data to be passed to plotting or QCs
        - results: is meant to contain elements or objects related to the pipeline progression
        - metrics: keeps track of values that are used to generate the final summary and compute stats 
    """
    results: dict  
    metrics: dict = {}

    def has_metrics(self):
        """Check if a step returned any metrics """
        return len(self.metrics)>0
