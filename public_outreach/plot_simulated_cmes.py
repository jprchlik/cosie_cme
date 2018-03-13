import matplotlib
matplotlib.rc('font',size=24)
import numpy as np
from fancy_plot import *
import pandas as pd
import matplotlib.pyplot as plt

kmsolar = 6.96*10.**5. #solar radii to km

#create observation obscured or unobscured Use 30s COSIE cadence
cad = '10S'
scl = float(cad[:-1]) #turn cadence into a float

#write out simulated cmes
cmes = pd.read_csv('simulated_cosie_cme.csv')

#get cme velocity bins
res = 100
bins = np.arange(0,2500,res)
bcme = cmes.groupby(np.digitize(cmes.v,bins))

#get cme obseration duration bins
res1 = 300
bins1 = np.arange(0,15000,res1)
bdur = cmes.groupby(np.digitize(cmes.obs_dur*scl,bins1))

#used bins for plotting
ubin = bins[bcme.obs_dur.mean().index.values-1]+(res/2.)

#bins for histogram plotting cme velocity
#velocity bins
hbin = np.array([[i,i+res] for i in bins[bcme.obs_dur.mean().index.values-1]]).ravel()
hbin = np.concatenate([[hbin[0]],hbin,[hbin[-1]]])
#Percent of total cmes
hplt = np.array([[i,i] for i in 100.*bcme.size()/len(cmes)]).ravel()
hplt = np.concatenate([[0],hplt,[0]])
#Switch to detection percentage based on email for Leon 2018/02/15
hplt = np.array([[i,i] for i in 100.*bcme.two_obs.sum()/bcme.two_obs.count()]).ravel()
hplt = np.concatenate([[0],hplt,[0]])

#bins for histogram plotting cme obs. duration
## Change to CDF 2018/02/13 J. Prchlik
##
hbin1 = cmes.obs_dur.sort_values()*scl
hplt1 = 100.*np.arange(len(hbin1))/float(len(hbin1)-1)

##THIS WAS DONE SO THE DEEP AIA IMAGE MATCHES##
theta1 = np.radians(90) #location of north in the image

fig2, ax2 = plt.subplots(nrows=2,ncols=2,gridspec_kw={'height_ratios':[1,2],'width_ratios':[2,1]},figsize=(12,12))
fig2.subplots_adjust(hspace=0.05,wspace=0.05)

#turn top right figure off
ax2[0,1].axis('off')

#plot scatter points
ax2[1,0].scatter(cmes.v,cmes.obs_dur*scl,color='black',label='Sim. CME')
ax2[1,0].errorbar(ubin,bcme.obs_dur.mean()*scl,yerr=bcme.obs_dur.std()*scl,fmt='s',color='red',label='Mean')
ax2[1,0].set_xlabel('CACTus Velocity [km/s]')
ax2[1,0].set_ylabel('COSIE Obs. Duration [s]')
ax2[1,0].legend(loc='upper right',scatterpoints=1,frameon=False,handletextpad=-0.12)

#plot velocity histogram
ax2[0,0].plot(hbin,hplt,color='black',linewidth=3)
ax2[0,0].set_xticklabels([])
#ax2[0,0].set_ylabel('Occurrence [%]')
#Switch to detection fraction based on email from Leon on 2018/02/15
ax2[0,0].set_ylabel('Detection Rate [%]')
ax2[0,0].set_xlim(ax2[1,0].get_xlim())
#Switch to detection fraction based on email from Leon on 2018/02/15
#ax2[0,0].set_ylim([0.,30.])
ax2[0,0].set_ylim([0.,100.])

#plot observation time histogram
ax2[1,1].plot(100.-hplt1,hbin1,color='black',linewidth=3)
ax2[1,1].set_yticklabels([])
ax2[1,1].set_xlabel('Cumulative COSIE \n Obs. Duration [%]')
ax2[1,1].set_ylim(ax2[1,0].get_ylim())
ax2[1,1].set_xlim([0.,105.])

fancy_plot(ax2[1,0])
fancy_plot(ax2[1,1])
fancy_plot(ax2[0,0])
fig2.savefig('cactus_durat_obs_vs_vel_reverse.png',bbox_pad=.1,bbox_inches='tight')
fig2.savefig('cactus_durat_obs_vs_vel_reverse.eps',bbox_pad=.1,bbox_inches='tight')
fig2.savefig('cactus_durat_obs_vs_vel_reverse.pdf',bbox_pad=.1,bbox_inches='tight')


