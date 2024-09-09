import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
mean_CCC = pd.read_csv("metric_meanCCC.tsv", sep='\t')

platform_protocols = ['cDNA_PacBio', 'CapTrap_PacBio', 'dRNA_ONT', 'cDNA_ONT']
indices = range(len(platform_protocols))
# Calculate optimal width
width = np.min(np.diff(indices))/4.

fig, ((ax3, ax4), (ax1, ax2)) = plt.subplots(nrows=2, ncols=2)
#fig = plt.figure()
#ax1 = fig.add_subplot(411)
#ax.bar(indices+2*width,mean_CCC[mean_CCC['tool'] == 'TALON']['meanCCC'],.5*width,color='C5',label='TALON')
ax1.bar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'],.5*width,color='C1',label='lr-kallisto')
ax1.errorbar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC_error'], fmt="o", color="C5")
ax1.bar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'],.5*width,color='C2',label='IsoQuant')
ax1.errorbar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC_error'], fmt="o", color="C5")
ax1.bar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'],.5*width,color='C3',label='Bambu')
ax1.errorbar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC_error'], fmt="o", color="C5")
ax1.bar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'],.5*width,color='C4',label='Oarfish')
ax1.errorbar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC_error'], fmt="o", color="C5")
#tiks = ax.get_xticks().tolist()
#ax1.axes.set_xticklabels(platform_protocols)
ax1.set_xlabel('Protocols and Platforms')
ax1.set_title("LRGASP Human WTC11",fontweight='bold')
ax1.set_xticklabels(platform_protocols, rotation=45);
#ax1.ylim((0,1.5))
ax1.tick_params(bottom = False) 
ax1.set_ylabel('CCC of expression')

ax1.legend(loc='lower center')
#ax1.savefig("ccc_human_WTC11.png", dpi=360)
#plt.show()
#plt.clf()

#fig = plt.figure()
indices = range(4)
# Calculate optimal width
width = np.min(np.diff(indices))/4.
#ax3 = fig.add_subplot(412)
ax3.bar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC'],width/2,color='C1',label='lr-kallisto')
ax3.errorbar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC_error'], fmt="o", color="C5")
ax3.bar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC'],width/2,color='C2',label='IsoQuant')
ax3.errorbar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC_error'], fmt="o", color="C5")
ax3.bar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC'],width/2,color='C3',label='Bambu')
ax3.errorbar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC_error'], fmt="o", color="C5")
ax3.bar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC'],width/2,color='C4',label='Oarfish')
ax3.errorbar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCC_error'], fmt="o", color="C5")
#tiks = ax.get_xticks().tolist()
#plt.axes.set_xticklabels(platform_protocols)
#ax3.set_xlabel('Protocols and Platforms')
#plt.xticks([0,1,2,3], platform_protocols, rotation=90)
ax3.set_title('LRGASP Mouse ES',fontweight='bold')
#plt.xticks([0,1,2,3], platform_protocols, rotation=90);
#ax2.ylim((0,1.5))
ax3.tick_params(bottom = False) 
ax3.set_ylabel('CCC of expression')

ax3.legend(loc='lower center')
#ax2.savefig("ccc_Mouse_ES.png", dpi=360)
#plt.show()
#plt.xticks([0,1,2,3], platform_protocols, rotation=90)

mean_CCC = pd.read_csv("metric_meanCCCofCV2.tsv", sep='\t')

platform_protocols = ['cDNA_PacBio', 'CapTrap_PacBio', 'dRNA_ONT', 'cDNA_ONT']
indices = range(len(platform_protocols))
# Calculate optimal width
width = np.min(np.diff(indices))/4.

#fig = plt.figure()
#ax3 = fig.add_subplot(412)
#ax.bar(indices-width/2,mean_CCC[mean_CCC['tool'] == 'TALON']['meanCCCofCOV2'],width,color='r',label='TALON')
#ax.bar(indices+width/2.,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width,color='b',label='lr-kallisto')
#ax.bar(indices+3*width/2.,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width,color='g',label='Oarfish')
ax2.bar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width/2,color='C1',label='lr-kallisto')
ax2.errorbar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax2.bar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width/2,color='C2',label='IsoQuant')
ax2.errorbar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax2.bar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width/2,color='C3',label='Bambu')
ax2.errorbar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax2.bar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width/2,color='C4',label='Oarfish')
ax2.errorbar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2_error'], fmt="o", color="C5")
#tiks = ax.get_xticks().tolist()
#ax3.axes.set_xticklabels(platform_protocols)
ax2.set_xlabel('Platform and Protocol')
ax2.set_title('LRGASP Human WTC11',fontweight='bold')
plt.xticks([0,1,2,3], platform_protocols, rotation=45);
#ax3.ylim((0,1))
ax2.set_ylabel('CCC of isoform CV$^2$')

ax2.legend(loc='lower center')

#plt.show()
#plt.clf()
#fig = plt.figure()
indices = range(4)
# Calculate optimal width
width = np.min(np.diff(indices))/4.
#ax4 = fig.add_subplot(422)
ax4.bar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'],width/2,color='C1',label='lr-kallisto')
ax4.errorbar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax4.bar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'],width/2,color='C2',label='IsoQuant')
ax4.errorbar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax4.bar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'],width/2,color='C3',label='Bambu')
ax4.errorbar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax4.bar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'],width/2,color='C4',label='Oarfish')
ax4.errorbar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2_error'], fmt="o", color="C5")
#tiks = ax.get_xticks().tolist()
#plt.axes.set_xticklabels(platform_protocols, rotation=90)
#ax4.set_xlabel('Platform and Protocol')
#ax4.ylim((0,1))
ax4.set_title('LRGASP Mouse ES',fontweight='bold')
plt.xticks([0,1,2,3], platform_protocols, rotation=45)

ax4.set_ylabel('CCC of isoform CV$^2$')

ax4.legend(loc='lower center')

plt.show()

