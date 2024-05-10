import os
import pandas as pd
from scipy.io import mmread
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

def c_ccc(y_pred, y_true):
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
        ccc = numerator / denominator
        return ccc

def comparison(long, short, ref):
        print(long, short, ref)
        count = mmread(short+'/matrix.abundance.tpm.mtx')
        labels = pd.read_csv(short+'/transcripts.txt', header=None, sep='\t')

        count_bus_sr = pd.DataFrame(count.todense().T, columns=['sr_bus_counts'])
        count_bus_sr['transcript_id'] = [labels.values[i][0] for i in range(np.shape(labels.values)[0])]
        count_bus_sr = count_bus_sr[count_bus_sr['sr_bus_counts'] > 0]
        count_bus_sr.to_csv(short+'_bus_sr_quant_tcc.tsv', sep="\t", columns=['transcript_id','sr_bus_counts'], header=1, index=0)
        
        count = mmread(long+'/matrix.abundance.tpm.mtx')
        labels = pd.read_csv(long+'/transcripts.txt', header=None, sep='\t')

        count_bus_lr = pd.DataFrame(count.todense().T, columns=['lr_bus_counts'])
        count_bus_lr['transcript_id'] = [labels.values[i][0] for i in range(np.shape(labels.values)[0])]

        count_bus_lr = count_bus_lr[count_bus_lr['lr_bus_counts'] > 0]
        count_bus_lr.to_csv(long+'_bus_lr_quant_tcc.tsv', sep="\t", columns=['transcript_id','lr_bus_counts'], header=1, index=0)

        count = count_bus_sr.merge(count_bus_lr, how='outer', on='transcript_id')
        count = count.fillna(0)

        x = (count['lr_bus_counts'])
        y = (count['sr_bus_counts'])
        x = np.array(np.log2(x+1)).flatten()
        y = np.array(np.log2(y+1)).flatten()
        print("lr-kallisto compared with short: "+short)
        mrd = np.abs(y - x) / y
        print('MRD', np.median(mrd))
        print('nrmse', np.sqrt(np.mean(np.square(x-y),axis=0)) / np.std(y, axis=0))

        r, p = scipy.stats.pearsonr(x,y)   # Pearson's r
        rho = scipy.stats.spearmanr(x,y).correlation   # Spearman's rho
        tau = scipy.stats.kendalltau(x, y)[0]  # Kendall's tau
        m, b = np.polyfit(x, y, 1)

        ccc = c_ccc(x, y)
        
        print(ccc)

        #print(sim+'_k-'+k)
        print("Pearson's r:\t", r)
        print("Spearman's rho:\t", rho)
        print("Kendall's tau:\t", tau)
        #plt.scatter(x, y, alpha=.05, edgecolor='black', linewidth=1)
        plt.hexbin(x, y, gridsize=100, cmap='jet', bins='log')
        #plt.xscale('log')
        #plt.yscale('log')
        plt.plot([i for i in range(0, 16)],[i for i in range(0, 16)], label='y=x\nCCC '+str(round(ccc,2))) 
        plt.colorbar(label='counts')
        if long.split("_")[2] == '13H':
            capture = 'non-exome capture'
        else:
            capture = 'exome capture'
        if short.split("_")[2] == '13H':
            capture_2 = 'non-exome capture'
        else:
            capture_2 = 'exome capture'
        if long.split("_")[1] == 'nanopore':
            platform = 'ONT'
        else:
            platform = 'Illumina'
        if short.split('_')[1] == 'nanopore':
            platform2 = 'ONT'
        else:
            platform2 = 'Illumina'
        plt.xlabel(capture + " " + platform + " transcript-level counts")
        plt.ylabel(capture_2 + " " + platform2 + " transcript-level counts")
        plt.title("lr-kallisto: aligning to "+ref)
        plt.plot(x, m*x + b, c='r', label='Regression line\nPearson R '+str(round(r,2))+'\nSpearman '+r"$\rho$"+' '+str(round(rho,2)))
        plt.legend(loc='upper left')
        plt.savefig("b01_aligned_"+long.split("_")[2]+"_"+long.split("_")[1]+"_vs_"+short.split("_")[1]+"_"+short.split("_")[2]+"_"+ref.replace('/','')+".png", dpi=300)
        plt.clf()

def comparison_oarfish(long, short, ref):
        count = mmread(short+'/matrix.abundance.tpm.mtx')
        #print(count)
        labels = pd.read_csv(short+'/transcripts.txt', header=None, sep='\t')
        #print(np.shape(labels.values))
    
        #print(count.todense())
        count_bus_sr = pd.DataFrame(count.todense().T, columns=['sr_bus_counts'])
        count_bus_sr['transcript_id'] = [labels.values[i][0] for i in range(np.shape(labels.values)[0])]
        # # Don't remove for LRGASP evaluation [labels.values[i][0].split('.')[0] for i in range(np.shape(labels.values)[0])]
        #count_bus.index.name = 'transcript_id'
        count_bus_sr = count_bus_sr[count_bus_sr['sr_bus_counts'] > 0]
        count_bus_sr.to_csv(short+'_bus_lr_init_quant_tcc.tsv', sep="\t", columns=['transcript_id','sr_bus_counts'], header=1, index=0)
        
        count = pd.read_csv(long, sep = '\t', names=['transcript_id', 'length', 'reads_count'], header=0)
        count['reads_count'] = count['reads_count']/np.sum(count['reads_count'])*1000000
        count_bus_lr = count[count['reads_count'] > 0]

        count = count_bus_sr.merge(count_bus_lr, how='outer', on='transcript_id')
        count = count.fillna(0)
        #print(count) 

        x = (count['reads_count'])
        y = (count['sr_bus_counts'])
        x = np.array(np.log2(x+1)).flatten()
        y = np.array(np.log2(y+1)).flatten()
        print("Oarfish: "+short)
        mrd = np.abs(y - x) / y
        print('MRD', np.median(mrd))
        print('nrmse', np.sqrt(np.mean(np.square(x-y),axis=0)) / np.std(y, axis=0))
        #y, x = y/np.sum(y)*1000000 + 1,x/np.sum(x)*1000000 + 1

        r, p = scipy.stats.pearsonr(x,y)   # Pearson's r
        rho = scipy.stats.spearmanr(x,y).correlation   # Spearman's rho
        tau = scipy.stats.kendalltau(x, y)[0]  # Kendall's tau
        m, b = np.polyfit(x, y, 1)

        ccc = c_ccc(x, y)
        
        print(ccc)
        #print(sim+'_k-'+k)
        print("Pearson's r:\t", r)
        print("Spearman's rho:\t", rho)
        print("Kendall's tau:\t", tau)
        #plt.scatter(x, y, alpha=.05, edgecolor='black', linewidth=1)
        plt.hexbin(x, y, gridsize=100, cmap='jet', bins='log')
        #plt.xscale('log')
        #plt.yscale('log')
        plt.plot([i for i in range(0, 16)],[i for i in range(0, 16)], label='y=x\nCCC '+str(round(ccc,2))) 
        plt.colorbar(label='counts')
        if long.split("_")[2] == '13H':
            capture = 'non-exome capture'
        else:
            capture = 'exome capture'
        if short.split("_")[2] == '13H':
            capture_2 = 'non-exome capture'
        else:
            capture_2 = 'exome capture'
        plt.xlabel(capture + " ONT transcript-level counts")
        plt.ylabel(capture_2 + " Illumina transcript-level counts")
        plt.title("Oarfish: aligning to C57BL/6J")
        plt.plot(x, m*x + b, c='r', label='Regression line\nPearson R '+str(round(r,2))+'\nSpearman '+r"$\rho$"+' '+str(round(rho,2)))
        plt.legend(loc='upper left')
        plt.savefig("oarfish_b01_aligned_"+long.split("_")[2]+"_"+long.split("_")[1]+"_vs_"+short.split("_")[1]+".png", dpi=300)
        plt.clf()

def comparison_bambu(long, short, ref):
        count = mmread(short+'/matrix.abundance.tpm.mtx')
        labels = pd.read_csv(short+'/transcripts.txt', header=None, sep='\t')

        count_bus_sr = pd.DataFrame(count.todense().T, columns=['sr_bus_counts'])
        count_bus_sr['transcript_id'] = [labels.values[i][0] for i in range(np.shape(labels.values)[0])]
        count_bus_sr = count_bus_sr[count_bus_sr['sr_bus_counts'] > 0]
        count_bus_sr.to_csv(short+'_bus_sr_quant_tcc.tsv', sep="\t", columns=['transcript_id','sr_bus_counts'], header=1, index=0)

        count = pd.read_csv(long+'/counts_transcript.txt', sep = '\t', names=['transcript_id', 'GENEID', 'reads_count'], header=0)
        count['reads_count'] = count['reads_count']/np.sum(count['reads_count'])*1000000
        count_bus_lr = count#[count['reads_count'] > 0]

        count = count_bus_sr.merge(count_bus_lr, how='inner', on='transcript_id')
        count = count.fillna(0)

        x = (count['reads_count'])
        y = (count['sr_bus_counts'])
        x = np.array(np.log2(x+1)).flatten()
        y = np.array(np.log2(y+1)).flatten()
        print("Bambu: "+short)
        mrd = np.abs(y - x) / y
        print('MRD', np.median(mrd))
        print('nrmse', np.sqrt(np.mean(np.square(x-y),axis=0)) / np.std(y, axis=0))

        r, p = scipy.stats.pearsonr(x,y)   # Pearson's r
        rho = scipy.stats.spearmanr(x,y).correlation   # Spearman's rho
        tau = scipy.stats.kendalltau(x, y)[0]  # Kendall's tau
        m, b = np.polyfit(x, y, 1)

        ccc = c_ccc(x, y)

        print(ccc)
        print("Pearson's r:\t", r)
        print("Spearman's rho:\t", rho)
        print("Kendall's tau:\t", tau)
        plt.hexbin(x, y, gridsize=100, cmap='jet', bins='log')
        plt.plot([i for i in range(0, 16)],[i for i in range(0, 16)], label='y=x\nCCC '+str(round(ccc,2)))
        plt.colorbar(label='counts')
        if long.split("_")[2] == '13H':
            capture = 'non-exome capture'
        else:
            capture = 'exome capture'
        if short.split("_")[2] == '13H':
            capture_2 = 'non-exome capture'
        else:
            capture_2 = 'exome capture'
        plt.xlabel(capture + " ONT transcript-level counts")
        plt.ylabel(capture_2 + " Illumina transcript-level counts")
        plt.title("Bambu: aligning to C57BL/6J")
        plt.plot(x, m*x + b, c='r', label='Regression line\nPearson R '+str(round(r,2))+'\nSpearman '+r"$\rho$"+' '+str(round(rho,2)))
        plt.legend(loc='upper left')
        plt.savefig("bambu_b01_aligned_"+long.split("_")[2]+"_"+long.split("_")[1]+"_vs_"+short.split("_")[1]+".png", dpi=300)
        plt.clf()

def comparison_IQ(long, short, ref):
        count = mmread(short+'/matrix.abundance.tpm.mtx')
        labels = pd.read_csv(short+'/transcripts.txt', header=None, sep='\t')

        count_bus_sr = pd.DataFrame(count.todense().T, columns=['sr_bus_counts'])
        count_bus_sr['transcript_id'] = [labels.values[i][0] for i in range(np.shape(labels.values)[0])]
        count_bus_sr = count_bus_sr[count_bus_sr['sr_bus_counts'] > 0]
        count_bus_sr.to_csv(short+'_bus_sr_quant_tcc.tsv', sep="\t", columns=['transcript_id','sr_bus_counts'], header=1, index=0)

        count = pd.read_csv(long+'/OUT/OUT.transcript_counts.tsv', sep = '\t', names=['transcript_id', 'reads_count'], header=0)
        count['reads_count'] = count['reads_count']/np.sum(count['reads_count'])*1000000
        count_bus_lr = count#[count['reads_count'] > 0]

        count = count_bus_sr.merge(count_bus_lr, how='inner', on='transcript_id')
        count = count.fillna(0)

        x = (count['reads_count'])
        y = (count['sr_bus_counts'])
        x = np.array(np.log2(x+1)).flatten()
        y = np.array(np.log2(y+1)).flatten()
        print("IsoQuant: "+short)
        mrd = np.abs(y - x) / y
        print('MRD', np.median(mrd))
        print('nrmse', np.sqrt(np.mean(np.square(x-y),axis=0)) / np.std(y, axis=0))

        r, p = scipy.stats.pearsonr(x,y)   # Pearson's r
        rho = scipy.stats.spearmanr(x,y).correlation   # Spearman's rho
        tau = scipy.stats.kendalltau(x, y)[0]  # Kendall's tau
        m, b = np.polyfit(x, y, 1)

        ccc = c_ccc(x, y)

        print(ccc)
        print("Pearson's r:\t", r)
        print("Spearman's rho:\t", rho)
        print("Kendall's tau:\t", tau)
        plt.hexbin(x, y, gridsize=100, cmap='jet', bins='log')
        plt.plot([i for i in range(0, 16)],[i for i in range(0, 16)], label='y=x\nCCC '+str(round(ccc,2)))
        plt.colorbar(label='counts')
        if long.split("_")[2] == '13H':
            capture = 'non-exome capture'
        else:
            capture = 'exome capture'
        if short.split("_")[2] == '13H':
            capture_2 = 'non-exome capture'
        else:
            capture_2 = 'exome capture'
        plt.xlabel(capture + " ONT transcript-level counts")
        plt.ylabel(capture_2 + " Illumina transcript-level counts")
        plt.title("IsoQuant: aligning to C57BL/6J")
        plt.plot(x, m*x + b, c='r', label='Regression line\nPearson R '+str(round(r,2))+'\nSpearman '+r"$\rho$"+' '+str(round(rho,2)))
        plt.legend(loc='upper left')
        plt.savefig("IsoQuant_b01_aligned_"+long.split("_")[2]+"_"+long.split("_")[1]+"_vs_"+short.split("_")[1]+".png", dpi=300)
        plt.clf()

comparison("b01_nanopore_13G_casteij_bulk","b01_next2_13G_casteij", 'CAST/EiJ')
comparison("b01_nanopore_13H_casteij_bulk","b01_next2_13H_casteij", 'CAST/EiJ')
#comparison("b01_nanopore_casteij-k-31","b01_next1_13G_casteij")
#comparison("b01_nanopore_casteij-k-31","b01_next2_13G_casteij")
comparison("b01_nanopore_13G","b01_next1_13G", 'C57BL/6J')
comparison("b01_nanopore_13G","b01_next2_13G", 'C57BL/6J')
comparison("b01_nanopore_13G", "b01_nanopore_13H", 'C57BL/6J')
#comparison("b01_nanopore_13G_demux", "b01_next1_13G")

comparison("b01_next1_13G", "b01_next2_13G", 'C57BL/6J')
comparison("b01_next1_13G", "b01_next1_13H", 'C57BL/6J')
comparison("b01_next2_13G", "b01_next2_13H", 'C57BL/6J')
comparison("b01_next1_13H", "b01_next2_13H", 'C57BL/6J')

comparison("b01_nanopore_13H", "b01_next1_13G", 'C57BL/6J')
comparison("b01_nanopore_13G", "b01_next1_13H", 'C57BL/6J')
comparison("b01_nanopore_13H","b01_next1_13H", 'C57BL/6J')
comparison("b01_nanopore_13H","b01_next2_13H", 'C57BL/6J')
comparison_oarfish("b01_nanopore_13H_oarfish.quant", "b01_next1_13H", 'C57BL/6J')
comparison_oarfish("b01_nanopore_13G_oarfish.quant", "b01_next1_13G", 'C57BL/6J')
comparison_bambu("b01_nanopore_13G_bambu", "b01_next1_13G", "C57BL/6J")
comparison_IQ("b01_nanopore_13G_IQ", "b01_next1_13G", "C57BL/6J")
comparison("b01_nanopore_13H_casteij_bulk","b01_next2_13H_casteij", 'CAST/EiJ')
comparison("b01_nanopore_13G_casteij_bulk","b01_next2_13G_casteij", 'CAST/EiJ')
