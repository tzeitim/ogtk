from ogtk.utils.db import Xp
from pathlib import Path
from typing import Any, List
from ogtk.utils.log import CustomLogger, Rlogger, call
from .plotting import PlotDB


class FractureXp(Xp):
    """Extended Xp class with pipeline-specific functionality"""
    steps: Any
    dry: bool
    make_test: bool
    logger: CustomLogger
    pro_datain: str
    modality: str
    umi_len: int
    rev_comp: bool
    anchor_ont: str
    samples: List[str]
    pro_workdir: str
    plotdb: PlotDB
    do_plot: bool    
    fracture: dict
    start_anchor: str
    end_anchor: str

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.steps = getattr(self, 'steps', None)
        self.dry = getattr(self, 'dry', False)
        self.make_test = getattr(self, 'make_test', False)


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
                self.logger.error(f"Some sample(s) matched more than {max_files}\n{sample_files}")
                raise ValueError("Verify that there are files that coorespond to all samples.")
            
        return sample_files

