import pickle
import ogtk
import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    
    parser.add_argument('csv', type=str, help="path to the CSV samplesheet. currently using 'bc', 'alias', 'mother', 'father', 'injected'")
    parser.add_argument('figdir', type=str, help="path to dump all figures for all experiments")
    parser.add_argument('-e', '--end', type=int, help="how many reads to load from the original fastqs. end = 0 means all reads", default=0) 
     
    
    args=parser.parse_args()
    
    #sdf = pd.read_csv('/home/polivar/src/projects/mltracer/conf/fbli_experiments/fbli_bulk_experiments.csv')
    sdf = pd.read_csv(args.csv)
    #'/home/polivar/src/projects/mltracer/conf/fbli_experiments/fbli_bulk_mstlt10x.csv')
    sdf['fullpath'] = sdf[['wd', 'r1']].apply(lambda x: '/'.join(x), axis=1)
    sdf['name'] = sdf[['bc', 'alias', 'mother', 'father', 'injected']].apply(lambda x: '_'.join(x).replace('series', ''), axis=1)
    xps = sdf.groupby('fullpath')
    end = args.end
    if end == 0:
        end = None
    fig_dir = args.figdir #
    #"/local/users/polivar/src/projects/mltracer/figures/"
    all_counts = {}
    
    for i in sorted(xps.groups):
        print(i)
        df = xps.get_group(i)
        samples_dict = dict(zip(df['bc'], df['name']))
        # load the readset
        rs = ogtk.UM.pfastq_collapse_UMI(
                    fastq_ifn1 = i,
                    fastq_ifn2 = i.replace("R1", "R2"),
                    umi_len = 28, end = end)
    
        #determine the counts
        umi_counts = rs[1].umi_counts()
        #generate umi indexes
        umis_per_sample = {}
        for bc in df['bc']:
            umis_per_sample[samples_dict[bc]] = [umi for umi in umi_counts.index if umi.endswith(bc)]
    
        # export full readset
        set_name = os.path.basename(i.replace("R1", "").replace(".fastq.gz", ""))
        pickle_ofn = i.replace("R1", "pairs").replace(".fastq.gz", ".readset.pickle")
        print(f'saving full readset {os.path.basename(pickle_ofn)}')
    
        # save individual readsets for each index
        for k,umi_ls in umis_per_sample.items():
            print(f'saving subrs {k}')
            sub_pickle_ofn = pickle_ofn.replace("pairs", k)
            with open(sub_pickle_ofn, 'wb') as sub_pickle:
                pickle.dump(ogtk.UM.subrs(name = k, umilist = umi_ls, rs =rs[1], keep_intact = False), sub_pickle)
        
        # keep the "undetermined"
        with open(pickle_ofn, 'wb') as full_pickle:
            pickle.dump(rs, full_pickle)
    
        # plot umi counts
        xp_name = df['xp'].to_list()[0] + "_" + set_name
        png_out = f"{fig_dir}/{xp_name}.png"
        plt.cla()
        fig = plt.figure(figsize= (15,5))
        ax  = fig.add_subplot(111)
    
    
        plt.title(xp_name)
        plt.yscale('log')
        plt.xscale('log')
        for bc,umis in umis_per_sample.items():
            all_counts[xp_name + bc] = umi_counts[umis]
            plt.plot(all_counts[xp_name + bc].to_list())
    
        box = ax.get_position()
        ax.set_position([0.3, 0.4, box.width*0.3, box.height])
    
        plt.legend(umis_per_sample.keys(), bbox_to_anchor = (1.0, 0.5))
        plt.grid(True)
        plt.savefig(png_out, bbox_inches='tight', dpi=100) #, loc=(1.1, 0.5))
    
