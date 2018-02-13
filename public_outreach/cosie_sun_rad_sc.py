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
import pandas as pd

kmsolar = 6.96*10.**5. #solar radii to km
#Calculate the allowed dv for LASCO to CME using available position angles 
def calc_min_v(i,solarad,theta1,res=0.5,detrad=2.5):
#LASCO X and Y values
#current solar radius at time of obs in solar arcseconds
    rad = detrad*solarrad #arcsec lasco C2 inner edge from http://adsabs.harvard.edu/abs/1995SoPh..162..357B

    LX,LY,paa = cal_ind_v(i['pa'],theta1,detrad,solarrad)
#total obs time
    obstime = float((i['end_dt']-i['start_dt']).seconds/60.) #minutes 
    dx = (LX-i['X'])/obstime #arcsec per min
    dy = (LY-i['Y'])/obstime #arcsec per min

    v = np.sqrt(dx**2+dy**2)

    return dx,dy,v
    

#Calculate Vx an Vy components of CME using position of detected CME in CaCTUS
def return_obs_vel_comp(i,solarrad,theta1,res=0.5,detrad=2.5):
    global kmsolar

    
    paa = np.arange(i['pa']-i['da']/2.,i['pa']+i['da']/2.,res)


#    LX,LY,paa = cal_ind_v(paa,theta1,detrad,solarrad)
#get the X and Y position from CACTus and return position angle
    LX,LY,paa = cal_ind_v(i['pa'],theta1,detrad,solarrad)


#Change from the CACTus C2 from initial 
    dx,dy = float(LX-i['X']),float(LY-i['Y'])
    theta_obs = atan2(dy,dx)

    vx,vy = i['v']*cos(theta_obs),i['v']*sin(theta_obs) #decompose velocity from CACTus catalog [km/s]

#Abosolute velocity
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

#All cactus CMEs
cmes = pd.read_csv('../filament/cme_lz.txt',sep='|',escapechar='#',skiprows=23)

#remove whitespaces from pandas headers
cmes.rename(columns=lambda x: x.strip(),inplace=True)


#Read in 2011 orbital file
obst = pd.read_csv('../filament/python_f_cosie_obs.csv',sep=',',skiprows=1)

#set the start time for all unobscured observations
obst['obs'] = 1

#convert start good observing time to datetime object
obst['start_dt'] = pd.to_datetime(obst.start)
obst.set_index(obst.start_dt,inplace=True)

#cut out region for testing 2011 purposes 2018/02/13 J. Prchlik
cmes.set_index(pd.to_datetime(cmes.t0),inplace=True)
cmes = cmes.loc['2011-01-01':'2012-01-01']
obst = obst.loc['2011-01-01':'2012-01-01']

#bad time pandas dataframe
badt = pd.DataFrame()
badt['start_dt'] = obst.start_dt+pd.to_timedelta(obst.end,unit='m')
badt['obs'] = 0
badt.set_index(badt.start_dt,inplace=True)


#drop start and end columns so I can stack and sort 
#Not needed
#obst.drop(['start','end'],inplace=True)


#Obscured and unobscured observer times in one array
tott = pd.concat([obst,badt])
tott.sort_index(inplace=True)


#create observation obscured or unobscured Use 30s COSIE cadence
cad = '30S'
relt = pd.DataFrame(index=pd.date_range(tott.index.min(),tott.index.max(),freq=cad),columns=['obs'])

#interpolate total observation time onto new grid forward filling all 1s and 0s
tott = tott.reindex(relt.index,method='ffill')



#remove whitespace from halo identifier
#cmes['halo?'].replace(r'\s\ \s*','',regex=True,inplace=True)

#remove all cmes without the Halo distinction
cmes = cmes[cmes['halo?'] != '']

d_fmt ='%Y/%m/%d %H:%M'

#turns string times into arrays and adds to pandas array
cmes['cactus_dt'] = pd.to_datetime(cmes.t0)

#inner size of lasco C2 coronagraph
#http://adsabs.harvard.edu/cgi-bin/nph-bib_query?1995SoPh..162..357B&db_key=AST
r_lasco = 2.2 #Rsun (assume at minimum distance)


#end cosie radius
r_cosie = 3.0

#travel time between lasco cme detection and limb
cmes['dt_limb'] = pd.to_timedelta(((2.2-1.)*kmsolar)/cmes.v,unit='s')
cmes['dt_cosf'] = pd.to_timedelta(((3.-2.2)*kmsolar)/cmes.v,unit='s')


#start and end times for cme assuming constant velocity
cmes['start_dt'] = cmes.cactus_dt-cmes.dt_limb
cmes['end_dt']   = cmes.cactus_dt+cmes.dt_cosf


#observed duration (get sum all observation times [0 when none, 1 when observed] from 2011 orbit)
cmes['obs_dur'] = cmes.apply(lambda x: tott.loc[((tott.index >= x.start_dt) & (tott.index <= x.end_dt)),'obs'].sum(),axis=1)


##THIS WAS DONE SO THE DEEP AIA IMAGE MATCHES##
theta1 = np.radians(90) #location of north in the image

fig1, ax1 = plt.subplots()

ax1.scatter(cmes.v,cmes.obs_dur,color='black')
ax1.set_xlabel('CACTus Velocity [km/s]')
ax1.set_ylabel('Number of Frames with CME [\#]')

fancy_plot(ax1)

fig1.savefig('cactus_frame_obs_vs_vel.png',bbox_pad=.1,bbox_inches='tight')


#CUT SHORT for testing purposes
#cmes = cmes[:3]
#Do not apply cacl. for now 2018/02/13 J. Prchlik
#cmes['solar_ang'] = cmes.t0.apply( sunpy.sun.solar_semidiameter_angular_size) #sunpy.sun.solar_semidiameter_angular_size(cmes.cactus_dt).value

