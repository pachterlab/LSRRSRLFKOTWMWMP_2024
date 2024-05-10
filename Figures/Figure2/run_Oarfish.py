import pandas as pd
import os
import time 

def run(species, platform, sample):
    f = open('/home/rloving/LRGASP_data/'+species+'_'+platform+'_'+sample+'_time.txt', 'a')
    RNA_seq_data = pd.read_excel('Extended Data Table 1: RNA-Seq Data Matrix.xlsx')
    if species == 'human':
        ref = 'LRGASP.human.cdna.fasta'
    else:
        ref = 'LRGASP.cdna.fasta'
    for i,row in RNA_seq_data.iterrows():
        if row['sample'] == sample and row['species'] == species and row['platform'] == platform and row['file_contents'] == 'reads':
            if row['library_prep'] == 'dRNA' and not os.path.exists(row['file_acc']+'.perlm.fastq.gz'):
                os.system('gunzip '+row['file_acc']+'.fastq.gz')
                os.system('perl -pe \'tr/uU/tT/ unless(/@+/)\''+' < '+row['file_acc']+'.fastq'+' > '+row['file_acc']+'.perlm.fastq')
                os.system('gzip '+row['file_acc']+'.perlm.fastq')
            if True: 
                location=species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_oarfish'
                tic = time.perf_counter()
                print(location, tic)
                if row['library_prep'] == 'dRNA':
                    os.system("./minimap2-2.28_x64-linux/minimap2 -t 32 -ax map-ont -N 100 "+ref+" /home/rloving/LRGASP_data/"+row['file_acc']+".perlm.fastq.gz > /home/rloving/LRGASP_data/"+row['file_acc']+".sam")
                    os.system("samtools view -bS -o /home/rloving/LRGASP_data/"+row['file_acc']+".bam /home/rloving/LRGASP_data/"+row['file_acc']+".sam")
                    os.system("oarfish --alignments /home/rloving/LRGASP_data/"+row['file_acc']+".bam --threads 32 --output "+location+" --model-coverage --filter-group no-filters")
                elif row['platform'] == 'ONT':
                    os.system("./minimap2-2.28_x64-linux/minimap2 -t 32 -ax map-ont -N 100 "+ref+" /home/rloving/LRGASP_data/"+row['file_acc']+".fastq.gz > /home/rloving/LRGASP_data/"+row['file_acc']+".sam")
                    os.system("samtools view -bS -o /home/rloving/LRGASP_data/"+row['file_acc']+".bam /home/rloving/LRGASP_data/"+row['file_acc']+".sam")
                    os.system("oarfish --alignments /home/rloving/LRGASP_data/"+row['file_acc']+".bam --threads 32 --output "+location+" --model-coverage --filter-group no-filters")
                else:
                    os.system("./minimap2-2.28_x64-linux/minimap2 -t 32 -ax map-pb -N 100 "+ref+" /home/rloving/LRGASP_data/"+row['file_acc']+".fastq.gz > /home/rloving/LRGASP_data/"+row['file_acc']+".sam")
                    os.system("samtools view -bS -o /home/rloving/LRGASP_data/"+row['file_acc']+".bam /home/rloving/LRGASP_data/"+row['file_acc']+".sam")
                    os.system("oarfish --alignments /home/rloving/LRGASP_data/"+row['file_acc']+".bam --threads 32 --output "+location+" --model-coverage --filter-group no-filters")
                toc = time.perf_counter()
                f.write(location+'\t'+str(toc-tic)+'\n')
    f.close()


run('mouse', 'PacBio', 'ES')

run('mouse', 'ONT', 'ES')

run('human', 'PacBio', 'WTC11')

run('human', 'ONT', 'WTC11')
