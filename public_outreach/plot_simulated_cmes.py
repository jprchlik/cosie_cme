import matplotlib as mpl
mpl.rc('font',size=24)
import numpy as np
from fancy_plot import *
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes,InsetPosition

kmsolar = 6.96*10.**5. #solar radii to km


#Eds requested custom color map
ccmap = mpl.colors.ListedColormap([ 'blue','green','yellow','red'])
#modification to be intensity order
ccmap = mpl.colors.ListedColormap([ 'yellow','green','blue','purple'])
ccmap = mpl.colors.ListedColormap([ 'purple','green','tomato','cyan'])
ccmap = mpl.colors.ListedColormap([ 'purple','red','limegreen','cyan'])
#ccmap.set_over('0.25')
#ccmap.set_under('1.0')

bounds = [np.log10(0.0003),np.log10(0.003),np.log10(0.03),np.log10(0.3),np.log10(3.)]
norm = mpl.colors.BoundaryNorm(bounds, ccmap.N)


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
bins1 = np.arange(-res1+1,15000,res1)
bdur = cmes.groupby(np.digitize(cmes.obs_dur*scl,bins1))


#make 2D histogram 2018/03/13
H,xedges,yedges = np.histogram2d(cmes.v,cmes.obs_dur*scl,bins=(bins,bins1))
H = H.T #transpose for plotting
#set up X,Y values
X, Y = np.meshgrid(xedges, yedges)

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


#normalize by year 2018/03/14 J. Prchlik
obs_sec = (pd.to_datetime(cmes.t0.max())-pd.to_datetime(cmes.t0.min())).total_seconds()
obs_yrs = obs_sec/(60.*60.*24.*365.25) #seconds to years

#plot scatter points
H_plt = np.log10(H/1000./obs_yrs)
#H_plt[np.isfinite(H_plt) == False] = -3.

#Switch to binned data 2018/03/13 J. Prchlik
#ax2[1,0].scatter(cmes.v,cmes.obs_dur*scl,color='black',label='Sim. CME')
#add max and min color values
#plotc = ax2[1,0].pcolormesh(X,Y,H_plt,label=None,cmap=plt.cm.viridis.reversed(),vmax=1.0,vmin=-3.)
#custom color map at Ed's request
plotc = ax2[1,0].pcolormesh(X,Y,H_plt,label=None,cmap=ccmap,norm=norm)
ax2[1,0].errorbar(ubin,bcme.obs_dur.mean()*scl,yerr=bcme.obs_dur.std()*scl,fmt='s',color='black',label='Mean')



#create mask for plotting hatching on images
#bounds = [np.log10(0.0003),np.log10(0.003),np.log10(0.03),np.log10(0.3),np.log10(3.)]
#hatch_list = ['','/','-','\\','x','o']
#for j,i in enumerate(bounds):
#
#    if j == 0: mask = np.ma.masked_greater_equal(H_plt,i)
#    if ((j > 0) & (j < len(bounds)-1)): mask = np.ma.masked_outside(H_plt,bounds[j-1],bounds[j])
#    if j == len(bounds)-1: mask = np.ma.masked_less_equal(H_plt,bounds[j-1])
#
#    ax2[1,0].pcolor(X, Y, mask, hatch=hatch_list[j], alpha=0.)
#

#Add horizontal 90 minute orbit and 500s
ax2[1,0].plot([-100,10000],[3600.,3600.],'--b',linewidth=3)
ax2[1,0].text(1000,3690,'3600s (1 ISS Orbit)',fontsize=18)
ax2[1,0].plot([-100,10000],[500.,500.],'--b',linewidth=3)
ax2[1,0].text(2150,550,'500s',fontsize=18)
#reset xlimit
ax2[1,0].set_xlim([0,2400])


ax2[1,0].set_xlabel('CACTus Velocity [km/s]')
ax2[1,0].set_ylabel('COSIE Obs. Duration [s]')
ax2[1,0].legend(loc='upper right',scatterpoints=1,frameon=False,handletextpad=-0.12)

#set plot limit for 2d histogram 2018/03/14 J. Prchlik
ax2[1,0].set_ylim([-500,9500])

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

#Add annotation of cumlative distribution
ax2[1,1].annotate('92% \n of CMEs\n observed', xy=(92, 1),  xycoords='data',
            xytext=(45,1499), textcoords='data',
            arrowprops=dict(facecolor='blue', shrink=0.0),
            horizontalalignment='right', verticalalignment='top',
            fontsize=18)
ax2[1,1].annotate('50% of CMEs\n observed $>$ 1800s', xy=(50, 1800),  xycoords='data',
            xytext=(55,6000), textcoords='data',
            arrowprops=dict(facecolor='blue', shrink=0.0),
            horizontalalignment='center', verticalalignment='top',
            fontsize=18)

ax2[1,1].set_xlabel('Cumulative COSIE \n Obs. Duration [%]')
ax2[1,1].set_ylim(ax2[1,0].get_ylim())
ax2[1,1].set_xlim([0.,105.])

#turn top right figure off
ax2[0,1].axis('off')
#turn top right figure into color bar
#cbar = fig2.colorbar(plotc,cax=ax2[0,1])
#use inset axis instead

#switch to axis locator 2018/03/21 J. Prchlik
axins = inset_axes(ax2[1,0],
                   width="5%",  # width = 30% of parent_bbox
                   height="40%",  # height : 1 inch
                   loc=9,borderpad=1.0)



cbar = fig2.colorbar(plotc,cax=axins)
cbar.set_label('Log(N$_\mathrm{CMEs}$/year)',fontsize=18)
cbar.set_ticks([-3.0,-2.0,-1.0,0.,1.0])

fancy_plot(ax2[1,0])
fancy_plot(ax2[1,1])
fancy_plot(ax2[0,0])
fig2.savefig('cactus_durat_obs_vs_vel_reverse_2Dhist.png',bbox_pad=.1,bbox_inches='tight')
fig2.savefig('cactus_durat_obs_vs_vel_reverse_2Dhist.eps',bbox_pad=.1,bbox_inches='tight')
fig2.savefig('cactus_durat_obs_vs_vel_reverse_2Dhist.pdf',bbox_pad=.1,bbox_inches='tight')


####################################################
####################################################
#Plot for number of CMEs over a two year period
fig3,ax3 = plt.subplots(figsize=(12,12))

#get cme velocity bins
#res = 500
#bins = np.logspace(np.log10(200.),np.log10(3000.),10)
#Switched to manual bins
bins = np.array([0.,200.,400.,600.,800.,1200.,1800.,2400.])
bcme = cmes.groupby(np.digitize(cmes.v,bins))

#bins for histogram plotting cme velocity
#velocity bins double up so we get 2 points per x value
hbin = np.sort(np.array(bins.tolist()*2))

#Number of CMEs per two year campgain 
hplt = np.array([[i,i] for i in 2.*bcme.two_obs.sum()/11./1000.]).ravel()
hplt = np.concatenate([[0],hplt,[0]])

#plot velocity histogram
ax3.plot(hbin,hplt,color='black',linewidth=3)

ax3.set_xlabel('CACTus Velocity [km/s]')
ax3.set_ylabel('COSIE Observed CMEs per 2 Years [#]')

ax3.set_ylim([0,60])
fancy_plot(ax3)
fig3.savefig('cactus_cme_two_year.png',bbox_pad=.1,bbox_inches='tight')
fig3.savefig('cactus_cme_two_year.eps',bbox_pad=.1,bbox_inches='tight')
fig3.savefig('cactus_cme_two_year.pdf',bbox_pad=.1,bbox_inches='tight')

####################################################
####################################################
#GET NUMBER OF OBITS PER SPEED
#Read in 2011 orbital file
obst = pd.read_csv('../filament/python_f_cosie_obs.csv',sep=',',skiprows=1)

#set the start time for all unobscured observations
obst['obs'] = np.arange(len(obst))

#convert start good observing time to datetime object
obst['start_dt'] = pd.to_datetime(obst.start)
obst.set_index(obst.start_dt,inplace=True)
obst = obst.loc['2011-01-01':'2012-01-01']

#bad time pandas dataframe Use bad time to back fill the orbit preventing over orbit counting
badt = pd.DataFrame()
badt['start_dt'] = obst.start_dt+pd.to_timedelta(obst.end,unit='m')
badt['obs'] = np.arange(len(obst))+1
badt.set_index(badt.start_dt,inplace=True)

#Obscured and unobscured observer times in one array
tott = pd.concat([obst,badt])
tott.sort_index(inplace=True)


#cut cmes to only include observed events
cmes = cmes[cmes.obs_dur > 1.]

#Get orbit number for start of observation
#set cmes index to start time temporarily
cmes.set_index(pd.to_datetime(cmes.start_dt),inplace=True)
cmes['o_start'] = tott.reindex(cmes.index,method='ffill').obs

#set cmes index to end time temporarily
cmes.set_index(pd.to_datetime(cmes.end_dt),inplace=True)
cmes['o_end'] = tott.reindex(cmes.index,method='ffill').obs

#Get the number of orbits from the end orbit number minus the start orbit number times where the cme is observed
cmes['orbits'] = (1+cmes.o_end-cmes.o_start)*(cmes.obs_dur > 1)


#bin up the orbits
#get cme velocity bins
res = 100
bins = np.arange(0,2500,res)
bcme = cmes.groupby(np.digitize(cmes.v,bins))

#Number of CMEs per two year campgain 
hplt = np.array([i for i in bcme.orbits.mean()]).ravel()
herr = np.array([i for i in bcme.orbits.std()]).ravel()
htot = np.array([i for i in bcme.orbits.count()]).ravel()
hbin = np.array([i for i in bcme.v.mean()]).ravel()


#Plot for number of CMEs over a two year period
fig4,ax4 = plt.subplots(figsize=(12,12))

#plot velocity histogram
ax4.scatter(hbin,hplt,color='black')
ax4.errorbar(hbin,hplt,yerr=herr/np.sqrt(htot),color='black',fmt='o')

ax4.set_xlabel('CACTus Velocity [km/s]')
#ax4.set_ylabel('SS Orbits Observed [#]')
ax4.set_ylabel('Single Orbit Observation [\%]')

ax4.set_xlim([0.,2400.])
ax4.set_ylim([0.,110.])

fancy_plot(ax4)
fig4.savefig('cactus_cme_orbits_obs.png',bbox_pad=.1,bbox_inches='tight')
fig4.savefig('cactus_cme_orbits_obs.eps',bbox_pad=.1,bbox_inches='tight')
fig4.savefig('cactus_cme_orbits_obs.pdf',bbox_pad=.1,bbox_inches='tight')


####################################################
