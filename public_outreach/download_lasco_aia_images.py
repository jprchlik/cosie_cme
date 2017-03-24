from sunpy.net.helioviewer import HelioviewerClient
import matplotlib.pyplot as plt
hv = HelioviewerClient()
from datetime import timedelta,datetime
import os




class dl_event:

    def __init__(self,start,end,cadence=60.,outfmt='%Y/%m/%d/',aia_wav='304',bdir='cme_images/',x0=0,y0=0,h0=1024,w0=1024,res=1.2,pad=20,water=True):
        self.start = start
        self.end   = end
        self.cad   = timedelta(seconds=cadence)
        self.outfmt= outfmt
        self.aia   = '{0:3d}'.format(int(float(aia_wav)))
        self.res   = res
        self.pad   = pad #number of extra observations to pad
        self.water = water


#set up directory structure
        self.bdir  = bdir
        if self.bdir[-1] != '/': self.bdir = self.bdir+'/'
        self.edir  = self.start.strftime(self.outfmt)
        if self.edir[-1] != '/': self.edir = self.bdir+'/'
        self.edir = self.edir+self.aia+'/'
        
        self.x0 = x0
        self.y0 = y0
        self.h0 = h0
        self.w0 = w0

        if os.path.exists(self.bdir+self.edir):
           print 'Path Already Exists...Proceeding'
        else:
           os.makedirs(self.bdir+self.edir) 
       
#loop over all times and download
    def loop_download(self):
        length = float((self.start-self.end).seconds)
        number = int(length/self.cad.seconds+1.+self.pad) #correct for python indexing and add 1 extra time step

        for i in range(number): self.download_files(i)
   
#download files
    def download_files(self,i):
        file = hv.download_png(self.start+i*self.cad,self.res, 
                               "[SDO,AIA,AIA,{0},1,100],[SOHO,LASCO,C2,white-light,1,100]".format(self.aia),
                               directory=self.bdir+self.edir,
                               x0=self.x0, y0=self.y0, width=self.h0, height=self.w0,watermark=self.water)
