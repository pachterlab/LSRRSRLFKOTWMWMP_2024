import pandas as pd
import os

def run(species, platform):
	RNA_seq_data = pd.read_excel('Extended Data Table 1: RNA-Seq Data Matrix.xlsx')
	for i,row in RNA_seq_data.iterrows():
		if row['sample'] == 'WTC11' and row['species'] == species and row['platform'] == platform and row['file_contents'] == 'reads':
			if row['library_prep'] == 'dRNA' and not os.path.exists(row['file_acc']+'.perlm.fastq.gz'):
				os.system('gunzip '+row['file_acc']+'.fastq.gz')
				os.system('perl -pe \'tr/uU/tT/ unless(/@+/)\''+' < '+row['file_acc']+'.fastq'+' > '+row['file_acc']+'.perlm.fastq')
				os.system('gzip '+row['file_acc']+'.perlm.fastq')
			if True: 
        location=species+'_'+row['sample']+'_'+platform+'_'+row['library_prep']+'_'+str(row['replicate'])+'_IQ'
				if row['library_prep'] == 'dRNA':
					os.system("../IsoQuant/isoquant.py -t 32 -d nanopore --stranded forward --fastq /home/rebekah/LRGASP_data/"+row['file_acc']+".perlm.fastq.gz --reference ../lrgasp_grch38_sirvs.fasta --genedb ../lrgasp_gencode_v38_sirvs.gtf --complete_genedb --output "+location+" --prefix "+location)
				elif row['platform'] == 'ONT':
					os.system("../IsoQuant/isoquant.py -t 32 -d nanopore --stranded forward --fastq /home/rebekah/LRGASP_data/"+row['file_acc']+".fastq.gz --reference ../lrgasp_grch38_sirvs.fasta --genedb ../lrgasp_gencode_v38_sirvs.gtf --complete_genedb --output "+location+" --prefix "+location)
				else:
					os.system("../IsoQuant/isoquant.py -t 32 -d pacbio --stranded forward --fastq /home/rebekah/LRGASP_data/"+row['file_acc']+".fastq.gz --reference ../lrgasp_grch38_sirvs.fasta --genedb ../lrgasp_gencode_v38_sirvs.gtf --complete_genedb --output "+location+" --prefix "+location)


run('human', 'PacBio')			
run('human', 'ONT')
