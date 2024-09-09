import numpy as np
import pandas as pd
import scipy.stats as st

def compute_ccc(file, rep1, rep2):
	abundance = pd.read_csv(file, sep='\t')
	y_true = abundance.iloc[:,rep1]
	y_pred = abundance.iloc[:,rep2]
	y_true, y_pred = y_true/np.sum(y_true)*1000000,y_pred/np.sum(y_pred)*1000000
	ccc = c_ccc(y_true, y_pred)
	return ccc

def c_ccc(y_true, y_pred):
	cor = np.corrcoef(y_true, y_pred)[0][1]
	# Means
	mean_true = np.mean(y_true)
	mean_pred = np.mean(y_pred)
	# Population variances
	var_true = np.var(y_true)
	var_pred = np.var(y_pred)
	# Population standard deviations
	sd_true = np.std(y_true)
	sd_pred = np.std(y_pred)
	# Calculate CCC
	numerator = 2 * cor * sd_true * sd_pred
	denominator = var_true + var_pred + (mean_true - mean_pred)**2
	if denominator == 0:
		denominator = 1
	ccc = numerator / denominator
	
	return ccc

def compute_mean_ccc_COV2(file, rep1, rep2, rep3):
	abundance = pd.read_csv(file, sep='\t') 
	abundance_col0 = abundance.columns[0]
	if 'mouse' in file:
		gene_id = pd.read_csv('bambu_IQ_quants/mouse_ES_ONT_CapTrap_bambu/counts_transcript.txt', sep='\t')
	else:
		gene_id = pd.read_csv('bambu_IQ_quants/human_WTC11_ONT_CapTrap_bambu/counts_transcript.txt', sep='\t') 
	abundance = abundance.merge(gene_id[['TXNAME','GENEID']],how='inner', left_on=abundance.columns[0], right_on=gene_id.columns[0])
	abundance = abundance.drop(columns=['TXNAME', abundance_col0], axis=1)
	cols = [abundance.columns[rep1-1], abundance.columns[rep2-1], abundance.columns[rep3-1]]
	abundance[cols] = 1000000*abundance[cols] / abundance[cols].sum()
	std = abundance.groupby(['GENEID']).std()
	mean = abundance.groupby(['GENEID']).mean()
	std = std.fillna(0)
	#mean = mean.fillna(1)
	mean[mean == 0] = 1
	cov = pd.DataFrame(mean)
	cov['1'] = std.iloc[:,[rep1-1]].values/mean.iloc[:,[rep1-1]].values
	cov['2'] = std.iloc[:,[rep2-1]].values/mean.iloc[:,[rep2-1]].values
	cov['3'] = std.iloc[:,[rep3-1]].values/mean.iloc[:,[rep3-1]].values
	cov.dropna()
	cov = cov.loc[~(cov==0).all(axis=1)]
	ccc1 = c_ccc(np.square(cov['1']), np.square(cov['2']))
	#print('CCC(1, 2) = ', ccc1)
	ccc2 = c_ccc(np.square(cov['2']), np.square(cov['3']))
	#print('CCC(2, 3) = ', ccc2)
	ccc3 = c_ccc(np.square(cov['1']), np.square(cov['3']))
	#print('CCC(1, 3) = ', ccc3)
	reps = [ccc1, ccc2, ccc3]
	interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
	return str(round(np.mean(reps), 4))+'\t'+str((interval[1]-interval[0])/2.)
	#return 0
begin='human_WTC11'
end='_end_abundance_filtered_reduced.tsv'
rep1 = compute_ccc('cDNA_PacBio'+end, 1, 2)
#print('CCC(1, 2) = ', rep1)
rep2 = compute_ccc('cDNA_PacBio'+end, 2, 3)
#print('CCC(2, 3) = ', rep2)
rep3 = compute_ccc('cDNA_PacBio'+end, 1, 3)
#print('CCC(1, 3) = ', rep3)
print('tool\tsample\tProtocol_Platform\tmeanCCC\tmeanCCC_error')
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1, 
              loc=np.mean(reps), 
              scale=st.sem(reps))
print('TALON\t'+begin+'\tcDNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc('CapTrap_PacBio'+end, 1, 2)
rep2 = compute_ccc('CapTrap_PacBio'+end, 2, 3)
rep3 = compute_ccc('CapTrap_PacBio'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1, 
              loc=np.mean(reps),
              scale=st.sem(reps))
print('TALON\t'+begin+'\tCapTrap_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc('cDNA_ONT'+end, 1, 2)
rep2 = compute_ccc('cDNA_ONT'+end, 2, 3)
rep3 = compute_ccc('cDNA_ONT'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1, 
              loc=np.mean(reps),
              scale=st.sem(reps))
print('TALON\t'+begin+'\tcDNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc('dRNA_ONT'+end, 1, 2)
#print('CCC(1, 2) = ', rep1)
rep2 = compute_ccc('dRNA_ONT'+end, 2, 3)
#print('CCC(2, 3) = ', rep2)
rep3 = compute_ccc('dRNA_ONT'+end, 1, 3)
#print('CCC(1, 3) = ', rep1)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1, 
              loc=np.mean(reps),
              scale=st.sem(reps))
print('TALON\t'+begin+'\tdRNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

begin='mouse_ES_'
end='_oarfish.tsv'
rep1 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 2)
#print('CCC(1, 2) = ', rep1)
rep2 = compute_ccc(begin+'PacBio_cDNA'+end, 2, 3)
#print('CCC(2, 3) = ', rep2)
rep3 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 3)
#print('CCC(1, 3) = ', rep3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1, 
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Oarfish\tmouse_ES\tcDNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_CapTrap'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Oarfish\tmouse_ES\tCapTrap_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Oarfish\tmouse_ES\tcDNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_dRNA'+end, 1, 2)
#print('CCC(1, 2) = ', rep1)
rep2 = compute_ccc(begin+'ONT_dRNA'+end, 2, 3)
#print('CCC(2, 3) = ', rep2)
rep3 = compute_ccc(begin+'ONT_dRNA'+end, 1, 3)
#print('CCC(1, 3) = ', rep3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Oarfish\tmouse_ES\tdRNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)


begin='human_WTC11_'
end='_oarfish.tsv'
rep1 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Oarfish\thuman_WTC11\tcDNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_CapTrap'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Oarfish\thuman_WTC11\tCapTrap_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Oarfish\thuman_WTC11\tcDNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_dRNA'+end, 1, 2)
#print('CCC(1, 2) = ', rep1)
rep2 = compute_ccc(begin+'ONT_dRNA'+end, 2, 3)
#print('CCC(2, 3) = ', rep2)
rep3 = compute_ccc(begin+'ONT_dRNA'+end, 1, 3)
#print('CCC(1, 3) = ', rep3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Oarfish\thuman_WTC11\tdRNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

begin='mouse_ES_'
end='_k-63_bus_quant_tcc.tsv'
rep1 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 2)
#print('CCC(1, 2) = ', rep1)
rep2 = compute_ccc(begin+'PacBio_cDNA'+end, 2, 3)
#print('CCC(2, 3) = ', rep2)
rep3 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 3)
#print('CCC(1, 3) = ', rep3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('lr-kallisto\tmouse_ES\tcDNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_CapTrap'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('lr-kallisto\tmouse_ES\tCapTrap_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('lr-kallisto\tmouse_ES\tcDNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_dRNA'+end, 1, 2)
#print('CCC(1, 2) = ', rep1)
rep2 = compute_ccc(begin+'ONT_dRNA'+end, 2, 3)
#print('CCC(2, 3) = ', rep2)
rep3 = compute_ccc(begin+'ONT_dRNA'+end, 1, 3)
#print('CCC(1, 3) = ', rep3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('lr-kallisto\tmouse_ES\tdRNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

begin='human_WTC11_'
end='_k-63_bus_quant_tcc.tsv'
rep1 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('lr-kallisto\thuman_WTC11\tcDNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_CapTrap'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('lr-kallisto\thuman_WTC11\tCapTrap_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('lr-kallisto\thuman_WTC11\tcDNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_dRNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_dRNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_dRNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('lr-kallisto\thuman_WTC11\tdRNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

begin='bambu_IQ_quants/mouse_ES_'
end='_bambu/counts_transcript_lrgasp.tsv'
rep1 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Bambu\tmouse_ES\tcDNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_CapTrap'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Bambu\tmouse_ES\tCapTrap_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Bambu\tmouse_ES\tcDNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_dRNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_dRNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_dRNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Bambu\tmouse_ES\tdRNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

begin='bambu_IQ_quants/mouse_ES_'
end='_IQ.tsv'
rep1 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('IsoQuant\tmouse_ES\tcDNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_CapTrap'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('IsoQuant\tmouse_ES\tCapTrap_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('IsoQuant\tmouse_ES\tcDNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_dRNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_dRNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_dRNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('IsoQuant\tmouse_ES\tdRNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

begin='bambu_IQ_quants/human_WTC11_'
end='_bambu/counts_transcript_lrgasp.tsv'
rep1 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Bambu\thuman_WTC11\tcDNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_CapTrap'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Bambu\thuman_WTC11\tCapTrap_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Bambu\thuman_WTC11\tcDNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_dRNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_dRNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_dRNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('Bambu\thuman_WTC11\tdRNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

begin='bambu_IQ_quants/human_WTC11_'
end='_IQ.tsv'
rep1 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('IsoQuant\thuman_WTC11\tcDNA_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 2)
rep2 = compute_ccc(begin+'PacBio_CapTrap'+end, 2, 3)
rep3 = compute_ccc(begin+'PacBio_CapTrap'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('IsoQuant\thuman_WTC11\tCapTrap_PacBio\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_cDNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_cDNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_cDNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('IsoQuant\thuman_WTC11\tcDNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

rep1 = compute_ccc(begin+'ONT_dRNA'+end, 1, 2)
rep2 = compute_ccc(begin+'ONT_dRNA'+end, 2, 3)
rep3 = compute_ccc(begin+'ONT_dRNA'+end, 1, 3)
reps = [rep1, rep2, rep3]
interval = st.t.interval(confidence=0.90, df=len(reps)-1,
              loc=np.mean(reps),
              scale=st.sem(reps))
print('IsoQuant\thuman_WTC11\tdRNA_ONT\t', round(np.mean(reps), 4),"\t", (interval[1]-interval[0])/2.)

begin='human_WTC11'
end='_end_abundance_filtered_reduced.tsv'
print('tool\tsample\tProtocol_Platform\tmeanCCCofCOV2\tmeanCCCofCOV2_error')

print('TALON\t'+begin+'\tcDNA_PacBio\t', compute_mean_ccc_COV2('cDNA_PacBio'+end, 1, 2, 3))

print('TALON\t'+begin+'\tCapTrap_PacBio\t', compute_mean_ccc_COV2('CapTrap_PacBio'+end, 1, 2, 3))

print('TALON\t'+begin+'\tcDNA_ONT\t', compute_mean_ccc_COV2('cDNA_ONT'+end, 1, 2, 3)) 

print('TALON\t'+begin+'\tdRNA_ONT\t', compute_mean_ccc_COV2('dRNA_ONT'+end, 1, 2, 3))

begin='mouse_ES_'
end='_oarfish.tsv'
print('Oarfish\tmouse_ES\tcDNA_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_cDNA'+end, 1, 2, 3))

print('Oarfish\tmouse_ES\tCapTrap_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_CapTrap'+end, 1, 2, 3))

print('Oarfish\tmouse_ES\tcDNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_cDNA'+end, 1, 2, 3))

print('Oarfish\tmouse_ES\tdRNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_dRNA'+end, 1, 2, 3))


begin='human_WTC11_'
end='_oarfish.tsv'
print('Oarfish\thuman_WTC11\tcDNA_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_cDNA'+end, 1, 2, 3))

print('Oarfish\thuman_WTC11\tCapTrap_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_CapTrap'+end, 1, 2, 3))

print('Oarfish\thuman_WTC11\tcDNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_cDNA'+end, 1, 2, 3))

print('Oarfish\thuman_WTC11\tdRNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_dRNA'+end, 1, 2, 3))

begin='mouse_ES_'
end='_k-63_bus_quant_tcc.tsv'
print('lr-kallisto\tmouse_ES\tcDNA_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_cDNA'+end, 1, 2, 3))

print('lr-kallisto\tmouse_ES\tCapTrap_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_CapTrap'+end, 1, 2, 3))

print('lr-kallisto\tmouse_ES\tcDNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_cDNA'+end, 1, 2, 3))

print('lr-kallisto\tmouse_ES\tdRNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_dRNA'+end, 1, 2, 3))

begin='human_WTC11_'
end='_k-63_bus_quant_tcc.tsv'
print('lr-kallisto\thuman_WTC11\tcDNA_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_cDNA'+end, 1, 2, 3))

print('lr-kallisto\thuman_WTC11\tCapTrap_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_CapTrap'+end, 1, 2, 3))

print('lr-kallisto\thuman_WTC11\tcDNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_cDNA'+end, 1, 2, 3))

print('lr-kallisto\thuman_WTC11\tdRNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_dRNA'+end, 1, 2, 3))

begin='bambu_IQ_quants/mouse_ES_'
end='_bambu/counts_transcript_lrgasp.tsv'
print('Bambu\tmouse_ES\tcDNA_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_cDNA'+end, 1, 2, 3))

print('Bambu\tmouse_ES\tCapTrap_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_CapTrap'+end, 1, 2, 3))

print('Bambu\tmouse_ES\tcDNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_cDNA'+end, 1, 2, 3))

print('Bambu\tmouse_ES\tdRNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_dRNA'+end, 1, 2, 3))

begin='bambu_IQ_quants/mouse_ES_'
end='_IQ.tsv'
print('IsoQuant\tmouse_ES\tcDNA_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_cDNA'+end, 1, 2, 3))

print('IsoQuant\tmouse_ES\tCapTrap_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_CapTrap'+end, 1, 2, 3))

print('IsoQuant\tmouse_ES\tcDNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_cDNA'+end, 1, 2, 3))

print('IsoQuant\tmouse_ES\tdRNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_dRNA'+end, 1, 2, 3))

h='human_WTC11'
begin='bambu_IQ_quants/'+h+'_'
end='_bambu/counts_transcript_lrgasp.tsv'
print('Bambu\t'+h+'\tcDNA_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_cDNA'+end, 1, 2, 3))

print('Bambu\t'+h+'\tCapTrap_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_CapTrap'+end, 1, 2, 3))

print('Bambu\t'+h+'\tcDNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_cDNA'+end, 1, 2, 3))

print('Bambu\t'+h+'\tdRNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_dRNA'+end, 1, 2, 3))

begin='bambu_IQ_quants/'+h+'_'
end='_IQ.tsv'
print('IsoQuant\t'+h+'\tcDNA_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_cDNA'+end, 1, 2, 3))

print('IsoQuant\t'+h+'\tCapTrap_PacBio\t', compute_mean_ccc_COV2(begin+'PacBio_CapTrap'+end, 1, 2, 3))

print('IsoQuant\t'+h+'\tcDNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_cDNA'+end, 1, 2, 3))

print('IsoQuant\t'+h+'\tdRNA_ONT\t', compute_mean_ccc_COV2(begin+'ONT_dRNA'+end, 1, 2, 3))
