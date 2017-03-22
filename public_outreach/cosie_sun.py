import matplotlib
matplotlib.rc('font',size=24)
import numpy as np
from fancy_plot import *
import matplotlib.pyplot as plt
from astropy.io import ascii
from sunpy.net import hek
from datetime import datetime,timedelta
import sunpy.sun
import matplotlib.image as mpimg
from math import atan2,cos,sin

kmsolar = 6.96*10.**5. #solar radii to km
#Calculate the allowed dv for LASCO to CME using available position angles 
def calc_min_v(i,solarad,theta1,res=0.5,detrad=1.5):
#LASCO X and Y values
#current solar radius at time of obs in solar arcseconds
    rad = detrad*solarrad #arcsec lasco C2 inner edge

#    paa = np.arange(pa-dpa/2.,pa+dpa/2.,res)
    LX,LY,paa = cal_ind_v(i['pa'],theta1,detrad,solarrad)
#total obs time
    obstime = (i['end_dt']-i['start_dt']).seconds/60. #minutes 
    dx = (LX-i['X'])/obstime #arcsec per min
    dy = (LY-i['Y'])/obstime #arcsec per min

    v = np.sqrt(dx**2+dy**2)

    return dx,dy,v
    

#Calculate Vx an Vy components of CME using position of detected CME in CaCTUS
def return_obs_vel_comp(i,solarrad,theta1,res=0.5,detrad=1.5):
    global kmsolar

    
    paa = np.arange(i['pa']-i['da']/2.,i['pa']+i['da']/2.,res)

    LX,LY,paa = cal_ind_v(paa,theta1,detrad,solarrad)
    LX,LY,paa = cal_ind_v(i['pa'],theta1,detrad,solarrad)


    dx,dy = float(LX-i['X']),float(LY-i['Y'])
    theta_obs = atan2(dy,dx)

    vx,vy = i['v']*cos(theta_obs),i['v']*sin(theta_obs) #km/s


    v = np.sqrt(vx**2+vy**2)

    return vx,vy


##CALCULATE X,Y VALUES IN ARCSEC
def cal_ind_v(paa,theta1,detrad,solarrad):

    LX,LY = solarrad*detrad*np.cos(np.radians(paa)+theta1),solarrad*detrad*np.sin(np.radians(paa)+theta1) #arcseconds

    return LX,LY,paa


#calculates datetime of column in array
def calc_dt(tab,strtime,fmt):

    for i in strtime: tab[i+'_dt'] = [datetime.strptime(date,fmt) for date in tab[i]]

    return tab

#client = hek.HEKClient()

#2011 events observed with CaCTUS and Filament McCatalog
cmes = ascii.read('../filament/cosie_cat.dat')
#start and end of unobscured observing times for hypothetical COSIE in 2011
obst = ascii.read('../filament/python_f_cosie_obs.csv',delimiter=',')

fmt ='%Y/%m/%dT%H:%M:%S'

#turns string times into arrays and adds to astropy table
cmes = calc_dt(cmes,['start','end'],fmt)
fmt = fmt.replace('T',' ')
obst = calc_dt(obst,['start'],fmt)

#no start and endtime in obst table so calc endtime manually
obst['end_dt'] = [obst['start_dt'][j]+timedelta(minutes=i) for j,i in enumerate(obst['end'])]


timeres = 0.5 #half minute is fast COSIE cadence (30s)

eyetest= False

fig1,ax1 = plt.subplots(figsize=(24,24))
#prevent annoying white box
fig1.subplots_adjust(left=0, bottom=0, right=1, top=1)

fig2,ax2 = plt.subplots(figsize=(24,24))
#prevent annoying white box
fig2.subplots_adjust(left=0, bottom=0, right=1, top=1)
##FANCY OFF LIMB DEEP AIA EXPOSURES##
#image = mpimg.imread("../photos/cosie_fov_color_ring.png")
#ax1.imshow(image,extent=[-3,3,-3,3])

#Add Radial Circles
r1 = plt.Circle((0,0),radius=1.0,color='gray',fill=True,linewidth=5,zorder=0)
r2 = plt.Circle((0,0),radius=2.0,color='black',fill=False,linewidth=5)
r3 = plt.Circle((0,0),radius=3.0,color='black',fill=False,linewidth=5)
#Add gray sun
s1 = plt.Rectangle((-3.3,-3.3),6.6,6.6,color='black',fill=False,linewidth=15)
ax1.add_patch(r1)
ax1.add_patch(r2)
ax1.add_patch(r3)
ax1.add_patch(s1)
#ax2.add_patch(r1)
#ax2.add_patch(r2)
#ax2.add_patch(r3)

#Add R2 and R3 labels
theta = np.radians(315)
x0,y0 = np.cos(theta),np.sin(theta)
ax1.text(2.*x0,2.*y0,'   2R$_\odot$',color='black',fontsize=34,ha='left',weight='bold')
ax1.text(3.*x0,3.*y0,'   3R$_\odot$',color='black',fontsize=34,ha='left',weight='bold')
ax1.text(-3.22,-3.22,'COSIE FOV',color='black',fontsize=34,ha='left',weight='bold')


#ax2.text(2.*x0,2.*y0,'   2R$_\odot$',color='black',fontsize=34,ha='left')
#ax2.text(3.*x0,3.*y0,'   3R$_\odot$',color='black',fontsize=34,ha='left')
#ax2.text(-2.92,-2.92,'COSIE FOV',color='black',fontsize=34,ha='left')
#

##THIS WAS DONE SO THE DEEP AIA IMAGE MATCHES##
theta1 = np.radians(97) #location of north in the image
#theta1 = np.radians(90) #location of north no image
##PUT A DOT ON NORTH FOR TESTING PURPOSES##
#x1,y1 = np.cos(0.+theta1),np.sin(0+theta1)
#ax1.scatter(x1,y1,color='red',s=75)


#CUT SHORT for testing purposes
#cmes = cmes[:3]


for ff,i in enumerate(cmes):
    print 'NEW CME [km/s]'

#initalize variables before while loop
    base = i['start_dt']
    neval = base
    obstest, = np.where((obst['start_dt'] <= neval ) & (obst['end_dt'] >= neval))

#Find where COSE Eclipses would occur during test 2011 observations
    obsed = []
    timea = []
    timea.append(neval)
    obsed.append(obstest.size)
#calc smallest velocity
    solarrad = sunpy.sun.solar_semidiameter_angular_size(i['start_dt']).value
    dx,dy,minv = calc_min_v(i,solarrad,theta1)
    print dx*kmsolar/60./solarrad,dy*kmsolar/60./solarrad
    
    vx,vy = return_obs_vel_comp(i,solarrad,theta1)

#    print i['v']
#    print vx,vy,np.sqrt(vx**2+vy**2)

    dx,dy = (np.array([vx,vy])/kmsolar)*(60.*solarrad) #(km/s)*(Rsun/km)*(s/min)*(arcsec/Rsun) = arcsec/min
    print 'AIA to CaCTUS velocity (Vx,Vy) = ({0:5.2f},{1:5.2f}) [km/s]'.format(dx*kmsolar/60./solarrad,dy*kmsolar/60./solarrad)
    print 'CaCTUS velocity (Vx,Vy) = ({0:5.2f},{1:5.2f}) [km/s]'.format(vx,vy)

#initalize x and y arrays
    xs = [i['X']]
    ys = [i['Y']]

  #loop variable 
    k=0 
#    while neval <= i['end_dt']:
    while ((np.abs(xs[-1]) < 3.3*solarrad) & (np.abs(ys[-1]) < 3.3*solarrad)):
        neval = neval+timedelta(minutes=timeres) 
        timea.append(neval)
        obstest, = np.where((obst['start_dt'] <= neval ) & (obst['end_dt'] >= neval))
        obsed.append(obstest.size)
 
        xs.append(i['X']+dx*k*timeres)
        ys.append(i['Y']+dy*k*timeres)

        k+=1

#convert x and y values into arrays
    xs = np.array(xs)/solarrad
    ys = np.array(ys)/solarrad

#separate array values into good and bad
    obsed = np.array(obsed)
    timea = np.array(timea)
    good = obsed == 1
    badd = obsed != 1

    cgood, = np.where(good)
    print "{0:4d} good obseravtions with {1:4.2f} min cad".format(cgood.size,timeres)

    LX,LY,paa = cal_ind_v(i['pa'],theta1,1.5,solarrad)


    ax1.scatter(xs[good],ys[good],color='blue',label=None)
    ax1.scatter(xs[badd],ys[badd],marker='x',color='red',label=None)


    ax2.scatter(xs[good],ys[good],color='blue',label=None)
    ax2.scatter(xs[badd],ys[badd],marker='x',color='red',label=None)
#    ax1.text(xs[-7],ys[-7],'{0:2d}'.format(ff))

#show where LASCO picks up CME
#   ax1.scatter(LX/solarrad,LY/solarrad,marker='s',color='green',s=100)

#The CME breakupfinder passes the eye teest
    if eyetest:
        fig,ax =plt.subplots()
        use, = np.where((obst['start_dt'] >= base) & (obst['end_dt'] <= i['end_dt']))
    
        ptest = np.vstack([obst['start_dt'][use],obst['end_dt'][use]])
        uvals = np.ones(ptest.shape)
    
    
        ax.scatter(timea[good],obsed[good],color='blue')
        ax.scatter(timea[badd],obsed[badd],color='red')
        for j in use:
            ax.plot([obst['start_dt'][j],obst['end_dt'][j]],[1,1],color='green')
    
        plt.show() 


##PUT POINTS OFF THE PLOT FOR LEGEND##
ax1.scatter(-1111,1111,color='blue',label='Observable',s=175)
ax1.scatter(-1111,1111,marker='x',color='red',label='Eclipsed',s=175)

##FORMAT LEGEND SO IT LOOKS NICE##
leg = ax1.legend(bbox_to_anchor=(0.00, 0.13),
           scatterpoints=1,loc=2,
           frameon=False,handletextpad=-.1,fontsize=36)


#SETS LEGEND COLOR TO WHITE FOR WHEN YOU USE THE AIA IMAGE
#for text in leg.get_texts(): text.set_color('w')

ax1.set_axis_off()
ax1.set_xlim([-3.330,3.330])
ax1.set_ylim([-3.330,3.330])
ax1.set_xlabel('X [R$_\odot$]',fontsize=44)
ax1.set_ylabel('Y [R$_\odot$]',fontsize=44)
#fancy_plot(ax1)
fig1.savefig('COSIE_CME_plot.png',bbox_pad=None,bbox_inches=0,dpi=200)
#plt.show()
    
         


    

#event_type='CME'
#for i in np.arange(cmes['start'].size):
#    result = client.query(hek.attrs.Time(cmes['start'][i].replace('T',' '),cmes['end'][i].replace('T',' ')),hek.attrs.EventType(event_type))
#    print result
