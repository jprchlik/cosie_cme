from astropy.io import ascii
from download_lasco_aia_images import dl_event
from datetime import datetime,timedelta

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

h0,w0 = 3000,3000
xoff, yoff = h0/3.,w0/3.
advstart = timedelta(minutes=30)

#loop and download all movies
for i in range(cmes['X'].size):

    if cmes['X'][i] < 0.: xoff = -xoff
    if cmes['Y'][i] < 0.: yoff = -yoff
    
    event = dl_event(cmes['start_dt'][i]-advstart,cmes['end_dt'][i],x0=cmes['X'][i]+xoff,y0=cmes['Y'][i]+yoff,h0=h0,w0=w0,res=1.2)
    event.loop_download()
    event.make_movie_files()


