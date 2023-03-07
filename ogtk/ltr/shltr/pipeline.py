from ogtk.utils import db
import ogtk.utils as ut
import ogtk.ltr.shltr as shltr
import scanpy as sc

class Xp(db.Xp):
    def load_10_mtx(self, cache = True, batch = None):
        ''' Loads the corresponding 10x matrix for a given experiment
        '''
        matrix_dir = f'{self.wd_scrna}/{self.sample_id}/outs/filtered_feature_bc_matrix/'
        self.print(f"importing matrix from:\n {matrix_dir}")

        adata = sc.read_10x_mtx(matrix_dir,
            var_names = 'gene_symbols',
            cache = cache
            )

        if batch is None:
            batch = self.sample_id

        adata.obs['batch'] = batch
        return(adata)

    def list_guide_tables(self,  suffix = 'shrna'):
        ''' Returns the file names of the available guide reads
        '''
        path = f'{self.wd_samplewd}/{suffix}'
        files = ut.sfind(path=path, pattern='*.parquet')
        return(files)

    def load_guide_molecules(
            self,
            clone,
            sample_id = None,
            suffix = 'shrna',
            index = 0,
            valid_ibars = None,
            down_sample = None,
            use_cache = True,
            corr_dir_fn=None):
        ''' Returns data frame at the molecule level
        '''
        if sample_id is None:
            sample_id = self.sample_id

        files = self.list_guide_tables(suffix)
        parquet_ifn = files[index]

        shltr.sc.reads_to_molecules(
            sample_id=sample_id,
            parquet_ifn=parquet_ifn,
            valid_ibars = valid_ibars, 
            use_cache=use_cache,
            corr_dict_fn=corr_dir_fn,
            down_sample=down_sample,
            #min_cells_per_ibar=int(len(obs27['cbc']) * 0.2),
            min_reads=2, 
            max_reads=1e6, 
            clone=clone)

    def init_mols(self):
        ''' Creates the .mols atribute
        '''
        self.mols = self.load_guide_molecules()
