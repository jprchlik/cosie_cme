from astropy.io import ascii
from download_lasco_aia_images import dl_event
from datetime import datetime

#calculates datetime of column in array
def calc_dt(tab,strtime,fmt):
    for i in strtime: tab[i+'_dt'] = [datetime.strptime(date,fmt) for date in tab[i]]
    return tab

#2011 events observed with CaCTUS and Filament McCatalog
cmes = ascii.read('../filament/cosie_cat.dat')

fmt ='%Y/%m/%dT%H:%M:%S'

#turns string times into arrays and adds to astropy table
cmes = calc_dt(cmes,['start','end'],fmt)
fmt = fmt.replace('T',' ')

event = dl_event(cmes['start_dt'][0],cmes['end_dt'][0],x0=cmes['X'][0],y0=cmes['Y'][0])
event.loop_download()


