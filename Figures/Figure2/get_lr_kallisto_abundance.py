from scipy.io import mmread
import pandas as pd
import numpy as np
import os

platforms = ['PacBio','ONT']
def LRGASP_comp(sim_name, species, k):
    RNA_seq_data = pd.read_excel('Extended Data Table 1: RNA-Seq Data Matrix.xlsx')
    for i,row in RNA_seq_data.iterrows():
        if row['species'] == species and row['platform'] == sim_name and "simulation" not in row['sample'] and row['file_contents'] == 'reads' and row['library_prep'] != 'R2C2':
            platform = row['platform']
            if os.path.exists(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/matrix.abundance.mtx'):
                count = mmread(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/matrix.abundance.mtx')

                labels = pd.read_csv(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/transcripts.txt', header=None, sep='\t')

                count_bus = pd.DataFrame(count.todense().T, columns=['bus_counts'])
                count_bus['transcript_id'] = [labels.values[i][0] for i in range(np.shape(labels.values)[0])]
                count_bus = count_bus[count_bus['bus_counts'] > 0]
                count_bus.to_csv('er_quant_tcc/'+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_bus_quant_tcc.tsv', sep="\t", columns=['transcript_id','bus_counts'], header=1, index=0)
                
                if row['replicate'] == 1:
                    count_all_reps = count_bus[['transcript_id', 'bus_counts']]
                else: 
                    count_all_reps = count_all_reps.merge(count_bus, how='outer', on='transcript_id')

            if os.path.exists(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/lr_quant_tcc/matrix.abundance.mtx'):
                count = mmread(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/lr_quant_tcc/matrix.abundance.mtx')

                labels = pd.read_csv(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/transcripts.txt', header=None, sep='\t')

                count_bus = pd.DataFrame(count.todense().T, columns=['bus_counts'])
                count_bus['transcript_id'] = [labels.values[i][0] for i in range(np.shape(labels.values)[0])]

                count_bus = count_bus[count_bus['bus_counts'] > 0]
                count_bus.to_csv('er_quant_tcc/'+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_bus_lr_quant_tcc.tsv', sep="\t", columns=['transcript_id','bus_counts'], header=1, index=0)
                if row['replicate'] == 1:
                    count_all_reps_lr = count_bus[['transcript_id', 'bus_counts']]
                else: 
                    count_all_reps_lr = count_all_reps_lr.merge(count_bus, how='outer', on='transcript_id')
            
            if os.path.exists(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/lr_quant_tcc/matrix.abundance.mtx') and row['replicate'] == 3:
                count_all_reps.fillna(0, inplace=True)
                count_all_reps_lr.fillna(0, inplace=True)
                count_all_reps.to_csv('er_quant_tcc/'+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_k-'+k+'_bus_quant_tcc.tsv', sep="\t", header=1, index=0)
                count_all_reps_lr.to_csv('er_quant_tcc/'+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_k-'+k+'_bus_lr_quant_tcc.tsv', sep="\t", header=1, index=0)    


LRGASP_comp("PacBio",'mouse',str(63))

LRGASP_comp("ONT",'mouse',str(63))

LRGASP_comp("PacBio",'human',str(63))

LRGASP_comp("ONT",'human',str(63))


