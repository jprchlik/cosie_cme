import matplotlib as mpl
mpl.rc('font',size=24)
import numpy as np
from fancy_plot import *
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes,InsetPosition

kmsolar = 6.96*10.**5. #solar radii to km

#All cactus CMEs
cmea = pd.read_csv('../filament/cme_lz.txt',sep='|',escapechar='#',skiprows=23)
#remove whitespaces from pandas headers
cmea.rename(columns=lambda x: x.strip(),inplace=True)
#cut out region for testing 2011 purposes 2018/02/13 J. Prchlik
cmea.set_index(pd.to_datetime(cmea.t0),inplace=True)
cmea = cmea.loc['2004-10-31':'2015-10-31']
#remove whitespace from halo identifier
cmea['halo?'].replace(r'\s\ \s*','',regex=True,inplace=True)


#get years in cme catalog for determining the uncertainty in the number of fast and slow cmes per year
#Comments from Ed 2018/08/30
uyears = np.unique(cmea.index.year).astype('str')
#Switched to third year frequency because we want to observe 10 CMEs
min_year= 525948.77
spt_year= str(int(round(min_year/3.)))
uyears = pd.date_range('2004-10-31','2015-10-31',freq=spt_year+'T')

#Eds requested custom color map
ccmap = mpl.colors.ListedColormap([ 'blue','green','yellow','red'])
#modification to be intensity order
#ccmap = mpl.colors.ListedColormap([ 'yellow','green','blue','purple'])
#ccmap = mpl.colors.ListedColormap([ 'purple','green','tomato','cyan'])
#ccmap = mpl.colors.ListedColormap([ 'purple','red','limegreen','cyan'])
#ccmap.set_over('0.25')
#ccmap.set_under('1.0')

bounds = [np.log10(0.0003),np.log10(0.003),np.log10(0.03),np.log10(0.3),np.log10(3.)]
norm = mpl.colors.BoundaryNorm(bounds, ccmap.N)


#create observation obscured or unobscured Use 30s COSIE cadence
cad = '10S'
scl = float(cad[:-1]) #turn cadence into a float



good_dur_arr = [30,35,40,45]
good_dur_arr = [30]
for good_dur in good_dur_arr:
    #read simulated cmes
    cmes = pd.read_csv('simulated_cosie_cme_good_dur_{0:1d}.csv'.format(good_dur))
    #cmes = pd.read_csv('simulated_cosie_cme.csv'.format(good_dur))
    #get cme velocity bins
    res = 100
    bins = np.arange(0,2500,res)
    bcme = cmes.groupby(np.digitize(cmes.v,bins))

    #get cme velocity bins for CME number plot
    #Switched back to raw data 2018/08/20 for the total number per year for COSIE
    res2 = 500
    bins2 = np.array([0.,500.,1000.,1500.,2100.])
    #real cme counts
    bcme2_r = cmea.groupby(np.digitize(cmea.v,bins2))
    #simulated cme counts which will turn into a fraction
    bcme2_s = cmes.groupby(np.digitize(cmes.v,bins2))


    #Added 2018/08/30 J. Prchlik per comments from Ed on 2018/08/29
    #get CME number by years 2018/08/30 J. Prchlik
    bcme2_y = []
    #loop over all years
    for j,cme_y in enumerate(uyears):

        #check for ending
        if j == len(uyears)-1:
            cmey = cmea[cme_y:]
        else:
            cmey = cmea[cme_y:uyears[j+1]]
        #real cme counts per third of a year (solved 0 in a bin problem)
        bcme2_yt = np.bincount(np.digitize(cmey.v,bins2),minlength=5)
        #remove 0 index
        hplt2_rt = bcme2_yt[1:]
        #simulated two observation detection fraction
        hplt2_s = np.array([i for i in bcme2_s.two_obs.sum()/bcme2_s.two_obs.count()]).ravel() 
        bcme2_y.append(hplt2_rt*hplt2_s)
    
    
    #convert to numpy array
    bcme2_y = np.array(bcme2_y)

    #Switch to low < 1000 and high < 1000 velocity bins
    unc_num = []
    #low velocity number uncertainty
    unc_num.append(np.std(np.sum(bcme2_y[:,:2],axis=1))/np.sqrt(uyears.size))
    #high velocity number uncertainty                                       
    unc_num.append(np.std(np.sum(bcme2_y[:,2:],axis=1))/np.sqrt(uyears.size))
    #error in mean
    #unc_num = np.std(bcme2_y,axis=0)/np.sqrt(uyears.size)
    #print the velocity bins and the uncertianty in the number in those bins by year
    print('Velocity Bins')
    print(bins2)
    print('Uncertainty in Number per Third Year')
    print(unc_num)

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

    #used bins for plotting
    #ubin2 = bins2[bcme2.obs_dur.mean().index.values-1]+(res2/2.)
    
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

    #Added CME velocity histogram per Ed's comment 2018/08/16
    #bins for histogram plotting cme velocity
    #velocity bins
    hbin2 = np.array([[i,i] for i in bins2[bcme2_s.obs_dur.mean().index.values]]).ravel()
    hbin2 = np.concatenate([[100],[100],hbin2])
    #hbin2[-1] = bins2[-1]
    
    #Simulated CME observation fraction 2018/08/20 J. Prchlik
    hplt2_s = np.array([[i,i] for i in bcme2_s.two_obs.sum()/bcme2_s.two_obs.count()]).ravel()
    hplt2_s = np.concatenate([[0],hplt2_s,[0]])

    #Real CME observation count 2018/08/20 J. Prchlik
    hplt2_r = np.array([[i,i] for i in bcme2_r.v.count()]).ravel()
    hplt2_r = np.concatenate([[0],hplt2_r,[0]])


    #combine real*simulated fraction
    hplt2 = hplt2_r*hplt2_s
    
    #bins for histogram plotting cme obs. duration
    ## Change to CDF 2018/02/13 J. Prchlik
    ##
    hbin1 = cmes.obs_dur.sort_values()*scl
    hplt1 = 100.*np.arange(len(hbin1))/float(len(hbin1)-1)
    
    ##THIS WAS DONE SO THE DEEP AIA IMAGE MATCHES##
    theta1 = np.radians(90) #location of north in the image
    
    #Updated to 2 panel plot based on Ed's request 2018/10/09
    #fig2, ax2 = plt.subplots(nrows=2,ncols=2,gridspec_kw={'height_ratios':[1,2],'width_ratios':[2,1]},figsize=(12,12))
    fig2, ax2 = plt.subplots(nrows=2,ncols=2,gridspec_kw={'height_ratios':[1,2],'width_ratios':[100,1]},figsize=(8,12))
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
    #removed at Ed's request 2018/08/16
    #plotc = ax2[1,0].pcolormesh(X,Y,H_plt,label=None,cmap=ccmap,norm=norm)
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
    #Remove horizontal line 2018/08/20
    ###ax2[1,0].plot([-100,10000],[5400.,5400.],'--',linewidth=3,color='black')
    ###ax2[1,0].text(1000,5200,'5400s (1 ISS Orbit)',fontsize=18)
    ###ax2[1,0].plot([-100,10000],[500.,500.],'--',linewidth=3,color='black')


    #Add crossing time
    cme_vel = np.linspace(100,2200,100,endpoint=True)
    cme_dur = 2.*kmsolar/cme_vel
    ax2[1,0].plot(cme_vel,cme_dur,'--',color='gray',linewidth=3)
    #ax2[1,0].plot(cme_vel,cme_dur*np.sqrt(2.),'--',color='gray',linewidth=3)
    
    ax2[1,0].text(1000,1500,'Crossing Time',color='gray')

    #changed to 50 when dropping p_colormesh 2018/08/16 J. Prchlik
    #Removed at Ed's requestion 2018/08/20
    #ax2[1,0].text(50,550,'500s',fontsize=18)
    #ax2[1,0].text(2150,550,'500s',fontsize=18)
    #reset xlimit
    ax2[1,0].set_xlim([0,2400])
    
    
    ax2[1,0].set_xlabel('CACTus Velocity [km/s]')
    ax2[1,0].set_ylabel('COSIE Obs. Duration [s]')
    ax2[1,0].legend(loc='upper right',scatterpoints=1,frameon=False,handletextpad=-0.12)
    
    #set plot limit for 2d histogram 2018/03/14 J. Prchlik
    #shorter axis range 2018/08/16
    ax2[1,0].set_ylim([-500,4900])
    
    #plot velocity histogram
    ax2[0,0].plot(hbin,hplt,color='black',linewidth=3,label='Detection Rate')
    ax2[0,0].set_xticklabels([])
    #ax2[0,0].set_ylabel('Occurrence [%]')
    #Switch to detection fraction based on email from Leon on 2018/02/15
    ax2[0,0].set_ylabel('Detection Rate [%]')
    ax2[0,0].set_xlim(ax2[1,0].get_xlim())
    #Switch to detection fraction based on email from Leon on 2018/02/15
    #ax2[0,0].set_ylim([0.,30.])
    ax2[0,0].set_ylim([0.,100.])


    #Add CMEs per yer in the top panel 2018/08/16 added fraction times real CMEs 2081/08/20
    ax3 = ax2[0,0].twinx()
    ax3.plot(hbin2,hplt2/obs_yrs,'--',color='blue',linewidth=3,label='# CMEs')
    #dummy plot for legend
    ax2[0,0].plot([-999,-999],[-999,-999],'--',color='blue',linewidth=3,label='# CMEs')
    ax3.set_ylabel('Obs. CMEs [# yr$^{-1}$]',color='blue')
    ax3.yaxis.label.set_color('blue')
    ax3.spines['right'].set_color('blue')
    ax3.tick_params(axis='y', colors='blue')
    #ax3.set_ylim([0,45])
    fancy_plot(ax3)
    ax3.yaxis.set_ticks_position('right')
    ax3.set_yscale('log')
    ax3.set_ylim([1,1E3])
    
    #Add legend
    ax2[0,0].legend(loc='upper right',frameon=False,fontsize=14)
    #plot observation time histogram
   
    #START Removed Cumlative distribution plot at Ed's request 2018/10/09 J. Prchlik
    ##ax2[1,1].plot(100.-hplt1,hbin1,color='black',linewidth=3)
    ##ax2[1,1].set_yticklabels([])
    ##

    ###get the percent observed
    ##thres_val = 360.
    ##obs_index, = np.where(hbin1 > thres_val)
    ##per_obs    = np.max(100.-hplt1[obs_index])
    ##
    ###get where you observe 50%
    ##fif_index, = np.where(np.abs(hplt1-50) == np.min(np.abs(hplt1-50)))
    ##fif_obs    = round(np.mean(hbin1.values[fif_index])/100)*100

    ###Add annotation of cumlative distribution
    ##ax2[1,1].annotate('{0:1.0f}% of    \n  CMEs      \n observed\n for {1:1.0f}s '.format(per_obs,thres_val), xy=(per_obs, thres_val+1),  xycoords='data',
    ##            xytext=(45,1499), textcoords='data',
    ##            arrowprops=dict(facecolor='blue', shrink=0.0),
    ##            horizontalalignment='right', verticalalignment='top',
    ##            fontsize=18)
    ###changed from 6000 to 3000 with resticted time range 2018/08/16 J. Prchlik
    ##ax2[1,1].annotate('50% of CMEs\n observed $>$ {0:1.0f}s'.format(fif_obs), xy=(50, fif_obs),  xycoords='data',
    ##            xytext=(55,3000), textcoords='data',
    ##            arrowprops=dict(facecolor='blue', shrink=0.0),
    ##            horizontalalignment='center', verticalalignment='top',
    ##            fontsize=18)
    ##
    ##ax2[1,1].set_xlabel('Cumulative COSIE \n Obs. Duration [%]')
    ##ax2[1,1].set_ylim(ax2[1,0].get_ylim())
    ##ax2[1,1].set_xlim([0.,105.])
    #removed 2018/08/30 from Ed's request on 2018/08/29
    #ax2[1,1].set_title('            {0:1d} Minutes of \n            Observations \n      per Orbit \n \n \n'.format(good_dur),fontsize=28)
    #END Removed Cumlative distribution plot at Ed's request 2018/10/09 J. Prchlik
    
    #turn top right figure off
    ax2[0,1].axis('off')
    #Turn bottom right figure off
    ax2[1,1].axis('off')
    #turn top right figure into color bar
    #cbar = fig2.colorbar(plotc,cax=ax2[0,1])
    #use inset axis instead
    
    #switch to axis locator 2018/03/21 J. Prchlik
    #removed at Ed's request 2018/08/16
    #axins = inset_axes(ax2[1,0],
    #                   width="5%",  # width = 30% of parent_bbox
    #                   height="40%",  # height : 1 inch
    #                   loc=9,borderpad=1.0)
    #
    #
    #
    #cbar = fig2.colorbar(plotc,cax=axins)
    #cbar.set_label('Log(N$_\mathrm{CMEs}$/year)',fontsize=18)
    #cbar.set_ticks([-3.0,-2.0,-1.0,0.,1.0])
    
    fancy_plot(ax2[1,0])
    #removed figure at Ed's request 2018/10/09
    #fancy_plot(ax2[1,1])
    fancy_plot(ax2[0,0])
    ax2[0,0].yaxis.set_ticks_position('left')
    fig2.savefig('cactus_durat_obs_vs_vel_reverse_2Dhist_{0:1d}.png'.format(good_dur),bbox_pad=.1,bbox_inches='tight')
    fig2.savefig('cactus_durat_obs_vs_vel_reverse_2Dhist_{0:1d}.eps'.format(good_dur),bbox_pad=.1,bbox_inches='tight')
    fig2.savefig('cactus_durat_obs_vs_vel_reverse_2Dhist_{0:1d}.pdf'.format(good_dur),bbox_pad=.1,bbox_inches='tight')
    
   #Not need 2018/08/14 J. Prchlik 
   ##### ####################################################
   ##### ####################################################
   ##### #Plot for number of CMEs over a two year period
   ##### fig3,ax3 = plt.subplots(figsize=(12,12))
   ##### 
   ##### #get cme velocity bins
   ##### #res = 500
   ##### #bins = np.logspace(np.log10(200.),np.log10(3000.),10)
   ##### #Switched to manual bins
   ##### bins = np.array([0.,200.,400.,600.,800.,1200.,1800.,2400.])
   ##### bcme = cmes.groupby(np.digitize(cmes.v,bins))
   ##### 
   ##### #bins for histogram plotting cme velocity
   ##### #velocity bins double up so we get 2 points per x value
   ##### hbin = np.sort(np.array(bins.tolist()*2))
   ##### 
   ##### #Number of CMEs per two year campgain 
   ##### hplt = np.array([[i,i] for i in 2.*bcme.two_obs.sum()/11./1000.]).ravel()
   ##### hplt = np.concatenate([[0],hplt,[0]])
   ##### 
   ##### #plot velocity histogram
   ##### ax3.plot(hbin,hplt,color='black',linewidth=3)
   ##### 
   ##### ax3.set_xlabel('CACTus Velocity [km/s]')
   ##### ax3.set_ylabel('COSIE Observed CMEs per 2 Years [#]')
   ##### 
   ##### ax3.set_ylim([0,60])
   ##### fancy_plot(ax3)
   ##### fig3.savefig('cactus_cme_two_year.png',bbox_pad=.1,bbox_inches='tight')
   ##### fig3.savefig('cactus_cme_two_year.eps',bbox_pad=.1,bbox_inches='tight')
   ##### fig3.savefig('cactus_cme_two_year.pdf',bbox_pad=.1,bbox_inches='tight')
   ##### 
   ##### ####################################################
   ##### ####################################################
   ##### #GET NUMBER OF OBITS PER SPEED
   ##### #Read in 2011 orbital file
   ##### obst = pd.read_csv('../filament/python_f_cosie_obs.csv',sep=',',skiprows=1)
   ##### 
   ##### #set the start time for all unobscured observations
   ##### obst['obs'] = np.arange(len(obst))
   ##### #parameter to fix over counting orbits
   ##### obst['fix'] = 0
   ##### 
   ##### #convert start good observing time to datetime object
   ##### obst['start_dt'] = pd.to_datetime(obst.start)
   ##### obst.set_index(obst.start_dt,inplace=True)
   ##### obst = obst.loc['2011-01-01':'2012-01-01']
   ##### 
   ##### #bad time pandas dataframe Use bad time to back fill the orbit preventing over orbit counting
   ##### badt = pd.DataFrame()
   ##### badt['start_dt'] = obst.start_dt+pd.to_timedelta(obst.end,unit='m')
   ##### badt['obs'] = np.arange(len(obst))+1
   ##### #parameter to fix over counting orbits
   ##### badt['fix'] = -1
   ##### badt.set_index(badt.start_dt,inplace=True)
   ##### 
   ##### #Obscured and unobscured observer times in one array
   ##### tott = pd.concat([obst,badt])
   ##### tott.sort_index(inplace=True)
   ##### 
   ##### 
   ##### #cut cmes to only include observed events
   ##### cmes = cmes[cmes.obs_dur > 1.]
   ##### 
   ##### #Get orbit number for start of observation
   ##### #set cmes index to start time temporarily
   ##### cmes.set_index(pd.to_datetime(cmes.start_dt),inplace=True)
   ##### cmes['o_start'] = tott.reindex(cmes.index,method='ffill').obs
   ##### 
   ##### #set cmes index to end time temporarily
   ##### cmes.set_index(pd.to_datetime(cmes.end_dt),inplace=True)
   ##### cmes['o_end'] = tott.reindex(cmes.index,method='ffill').obs
   ##### cmes['fix'] = tott.reindex(cmes.index,method='ffill').fix
   ##### 
   ##### #Get the number of orbits from the end orbit number minus the start orbit number times where the cme is observed
   ##### cmes['orbits_t'] = (1+cmes.o_end-cmes.o_start+cmes.fix)
   ##### 
   ##### #Count the 1 orbit fraction
   ##### cmes['orbits'] = cmes.orbits_t*(cmes.orbits_t < 1.1)
   ##### 
   ##### 
   ##### #bin up the orbits
   ##### #get cme velocity bins
   ##### res = 100
   ##### bins = np.arange(0,2500,res)
   ##### bcme = cmes.groupby(np.digitize(cmes.v,bins))
   ##### 
   ##### #Number of CMEs per two year campgain 
   ##### hplt = np.array([i for i in bcme.orbits.sum()]).ravel()
   ##### herr = np.array([i for i in bcme.orbits.std()]).ravel()
   ##### htot = np.array([i for i in bcme.orbits.count()]).ravel()
   ##### hbin = np.array([i for i in bcme.v.mean()]).ravel()
   ##### 
   ##### 
   ##### #Plot for number of CMEs over a two year period
   ##### fig4,ax4 = plt.subplots(figsize=(8,8))
   ##### 
   ##### #plot velocity histogram
   ##### ax4.scatter(hbin,hplt/htot*100.,color='black')
   ##### #ax4.errorbar(hbin,hplt/htot*100.,yerr=herr/np.sqrt(htot),color='black',fmt='o')
   ##### 
   ##### ax4.set_xlabel('CACTus Velocity [km/s]')
   ##### #ax4.set_ylabel('SS Orbits Observed [#]')
   ##### ax4.set_ylabel('Single Orbit Obs. Fraction [%]')
   ##### 
   ##### ax4.set_xlim([0.,2400.])
   ##### ax4.set_ylim([-5.,110.])
   ##### 
   ##### fancy_plot(ax4)
   ##### fig4.savefig('cactus_cme_orbits_obs.png',bbox_pad=.1,bbox_inches='tight')
   ##### fig4.savefig('cactus_cme_orbits_obs.eps',bbox_pad=.1,bbox_inches='tight')
   ##### fig4.savefig('cactus_cme_orbits_obs.pdf',bbox_pad=.1,bbox_inches='tight')
   ##### 
   ##### 
   ##### ####################################################
