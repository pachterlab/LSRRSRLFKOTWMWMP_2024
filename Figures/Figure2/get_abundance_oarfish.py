import pandas as pd
import numpy as np
import os

platforms = ['PacBio','ONT']
dir='v0.5.1LRGASP_data/'
def LRGASP_comp(sim_name, species):
    RNA_seq_data = pd.read_excel('Extended Data Table 1: RNA-Seq Data Matrix.xlsx')
    for i,row in RNA_seq_data.iterrows():
        platform=row['platform']
        if row['species'] == species and row['platform'] == sim_name and "simulation" not in row['sample'] and row['file_contents'] == 'reads' and row['library_prep'] != 'R2C2':
            location = species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_oarfish'
            if row['replicate'] == 3 and os.path.exists(dir+location+'.quant'):
                location = species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(1)+'_oarfish'
                filename = dir+location+'.quant'
                print(filename)
                count1 = pd.read_csv(filename, sep='\t')
                count1 = count1[count1['num_reads'] > 0]
                location = species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(2)+'_oarfish'
                count2 = pd.read_csv(dir+location+'.quant', sep='\t')
                count2 = count2[count2['num_reads'] > 0]
                location = species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(3)+'_oarfish'
                count3 = pd.read_csv(dir+location+'.quant', sep='\t')
                count3 = count3[count3['num_reads'] > 0]

                location = species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']
                count_all_reps_lr = count1[['tname', 'num_reads']].merge(count2[['tname','num_reads']], how='outer', on='tname')
                count_all_reps_lr = count_all_reps_lr.merge(count3[['tname','num_reads']], how='outer', on='tname')
                count_all_reps_lr.fillna(0, inplace=True)
                count_all_reps_lr.to_csv(location+'_oarfish.tsv', sep="\t", header=1, index=0)    


LRGASP_comp("PacBio",'mouse')

LRGASP_comp("ONT",'mouse')

LRGASP_comp("PacBio",'human')

LRGASP_comp("ONT",'human')


