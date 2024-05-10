import pandas as pd
import os

kallisto='../kallisto_synced/kallisto/build/src/kallisto'
bustools='../bustools/build/src/bustools'
def run_human(species, platform, k):
	RNA_seq_data = pd.read_excel('Extended Data Table 1: RNA-Seq Data Matrix.xlsx')
	for i,row in RNA_seq_data.iterrows():
		print(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er')
		if row['sample'] == 'WTC11' and row['species'] == species and row['platform'] == platform and row['file_contents'] == 'reads':
			if row['library_prep'] == 'dRNA' and not os.path.exists(row['file_acc']+'.perlm.fastq.gz'):
				os.system('gunzip '+row['file_acc']+'.fastq.gz')
				os.system('perl -pe \'tr/uU/tT/ unless(/@+/)\''+' < '+row['file_acc']+'.fastq'+' > '+row['file_acc']+'.perlm.fastq')
				os.system('gzip '+row['file_acc']+'.perlm.fastq')
			if True: 
        if row['library_prep'] == 'dRNA':
					os.system(kallisto+' bus -x bulk --error-rate .003 --threshold 0.8 --unmapped -t 30 --long -i LRGASP.human.k-'+k+'.dl.idx '+row['file_acc']+'.perlm.fastq.gz -o '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er')
				elif row['platform'] == 'ONT': 
					os.system(kallisto+' bus -x bulk --error-rate .004 --threshold 0.8 --unmapped -t 30 --long -i LRGASP.human.k-'+k+'.dl.idx '+row['file_acc']+'.fastq.gz -o '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er')
				else:
					os.system(kallisto+' bus -x bulk --error-rate .0009 --threshold 0.8 --unmapped -t 30 --long -i LRGASP.human.k-'+k+'.dl.idx '+row['file_acc']+'.fastq.gz -o '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er')
				os.system(bustools+' sort -t 30 '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/output.bus '\
					'-o '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/sorted.bus')
				os.system(bustools+' count '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/sorted.bus '\
					'-t '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/transcripts.txt '\
					'-e '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/matrix.ec '\
					'-o '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/count --cm -m '\
					'-g LRGASP.human.t2g')
				os.system(kallisto+' quant-tcc -t 30 '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/count.mtx '\
					'-i LRGASP.human.k-'+k+'.dl.idx '\
					'-f '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/flens.txt '\
					'-e '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/count.ec.txt '\
					'-o '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/lr_quant_tcc')
				os.system(kallisto+' quant-tcc -t 30 '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/count.mtx '\
					'-i LRGASP.human.k-'+k+'.dl.idx '\
					'-e '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/count.ec.txt '\
					'-o '+species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er')

def run(species, platform, k):
	f = open('/home/rebekah/LRGASP_data/'+species+'_'+platform+'_ES_time_lr-kallisto.txt', 'w')
	RNA_seq_data = pd.read_excel('Extended Data Table 1: RNA-Seq Data Matrix.xlsx')
	for i,row in RNA_seq_data.iterrows():
		if row['species'] == species and row['platform'] == platform and row['file_contents'] == 'reads':
			print(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er')
			tic = time.perf_counter()
			if row['library_prep'] == 'dRNA' and not os.path.exists(row['file_acc']+'.perlm.fastq.gz'):
				os.system('gunzip '+row['file_acc']+'.fastq.gz')
				os.system('perl -pe \'tr/uU/tT/ unless(/@+/)\''+' < '+row['file_acc']+'.fastq'+' > '+row['file_acc']+'.perlm.fastq')
				os.system('gzip '+row['file_acc']+'.perlm.fastq')
				#os.system('cp '+row['file_acc']+'.perlm.fastq.gz '+row['file_acc']+'.fastq.gz')
			if True:#not os.path.exists(species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_k-'+k+'_er/output.bus'):
				location=species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])
				if row['library_prep'] == 'dRNA':
					os.system(kallisto+' bus -x bulk --threshold 0.8 -t 32 --long -i ../LRGASP.k-'+k+'.dl.idx '+row['file_acc']+'.perlm.fastq.gz -o '+location+'_k-'+k+'_er')
				elif row['platform'] == 'ONT': 
					os.system(kallisto+' bus -x bulk --threshold 0.8 -t 32 --long -i ../LRGASP.k-'+k+'.dl.idx '+row['file_acc']+'.fastq.gz -o '+location+'_k-'+k+'_er')
				else:
					os.system(kallisto+' bus -x bulk --threshold 0.8 -t 32 --long -i ../LRGASP.k-'+k+'.dl.idx '+row['file_acc']+'.fastq.gz -o '+location+'_k-'+k+'_er')
				os.system(bustools+' sort -t 30 '+location+'_k-'+k+'_er/output.bus '\
					'-o '+location+'_k-'+k+'_er/sorted.bus')
				os.system(bustools+' count '+location+'_k-'+k+'_er/sorted.bus '\
					'-t '+location+'_k-'+k+'_er/transcripts.txt '\
					'-e '+location+'_k-'+k+'_er/matrix.ec '\
					'-o '+location+'_k-'+k+'_er/count --cm -m '\
					'-g ../LRGASP.t2g')
				os.system(kallisto+' quant-tcc -t 32 --long -P '+platform+' '+location+'_k-'+k+'_er/count.mtx '\
					'-i ../LRGASP.k-'+k+'.dl.idx '\
					'-e '+location+'_k-'+k+'_er/count.ec.txt '\
					'-o '+location+'_k-'+k+'_er')
			toc = time.perf_counter()
			f.write(location+'\t'+str(toc-tic)+'\n')
	f.close()

run('mouse', 'PacBio', str(63))			

run('mouse', 'ONT', str(63))

run_human('human', 'PacBio', str(63))			

run_human('human', 'ONT', str(63))
