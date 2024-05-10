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
ax1.bar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'],.5*width,color='C1',label='lr-kallisto')
ax1.errorbar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC_error'], fmt="o", color="C5")
ax1.bar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'],.5*width,color='C2',label='IsoQuant')
ax1.errorbar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC_error'], fmt="o", color="C5")
ax1.bar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'],.5*width,color='C3',label='Bambu')
ax1.errorbar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC_error'], fmt="o", color="C5")
ax1.bar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'],.5*width,color='C4',label='Oarfish')
ax1.errorbar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC'], yerr=mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCC_error'], fmt="o", color="C5")
ax1.set_xlabel('Protocols and Platforms')
ax1.set_title("LRGASP Human WTC11")
plt.xticks([0,1,2,3], platform_protocols, rotation=90);
ax1.tick_params(bottom = False) 
ax1.set_ylabel('CCC of expression')

ax1.legend(loc='lower center')

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
ax3.set_title('LRGASP Mouse ES')
ax3.tick_params(bottom = False) 
ax3.set_ylabel('CCC of expression')

ax3.legend(loc='lower center')

mean_CCC = pd.read_csv("metric_meanCCCofCV2.tsv", sep='\t')

platform_protocols = ['cDNA_PacBio', 'CapTrap_PacBio', 'dRNA_ONT', 'cDNA_ONT']
indices = range(len(platform_protocols))
# Calculate optimal width
width = np.min(np.diff(indices))/4.

ax2.bar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width/2,color='C1',label='lr-kallisto')
ax2.errorbar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax2.bar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width/2,color='C2',label='IsoQuant')
ax2.errorbar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax2.bar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width/2,color='C3',label='Bambu')
ax2.errorbar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax2.bar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'],width/2,color='C4',label='Oarfish')
ax2.errorbar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'human_WTC11')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax2.set_xlabel('Platform and Protocol')
ax2.set_title('LRGASP Human WTC11')
ax2.set_ylabel('CCC of isoform CV$^2$')

ax2.legend(loc='lower center')

indices = range(4)
# Calculate optimal width
width = np.min(np.diff(indices))/4.
ax4.bar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'],width/2,color='C1',label='lr-kallisto')
ax4.errorbar(indices-width,mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'lr-kallisto') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax4.bar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'],width/2,color='C2',label='IsoQuant')
ax4.errorbar(indices-width/3.,mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'IsoQuant') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax4.bar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'],width/2,color='C3',label='Bambu')
ax4.errorbar(indices+width/3.,mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'Bambu') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax4.bar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'],width/2,color='C4',label='Oarfish')
ax4.errorbar(indices+width,mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2'], yerr=mean_CCC[(mean_CCC['tool'] == 'Oarfish') & (mean_CCC['sample'] == 'mouse_ES')]['meanCCCofCOV2_error'], fmt="o", color="C5")
ax4.set_title('LRGASP Mouse ES')
plt.xticks([0,1,2,3], platform_protocols, rotation=90)

ax4.set_ylabel('CCC of isoform CV$^2$')

ax4.legend(loc='lower center')

plt.show()
