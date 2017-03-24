from sunpy.net.helioviewer import HelioviewerClient
import matplotlib.pyplot as plt
hv = HelioviewerClient()
from datetime import timedelta,datetime




class dl_event:

    def __init__(self,start,end,cadence=60.,outfmt='%Y/%m/%d/',aia_wav='304',bdir='cme_images/',x0=0,y0=0,h0=1024,w0=1024):
        self.start = start
        self.end   = end
        self.cad   = cadence
        self.outfmt= outfmt
        self.aia   = '{0:3d}'.format(round(float(aia_wav)))
        self.bdir  = bdir
        
        self.x0 = x0
        self.y0 = y0
        self.h0 = h0
        self.w0 = w0
       


   
    def download_files(self)

        file = hv.download_png(, 6, "[SDO,AIA,AIA,304,1,100],[SOHO,LASCO,C2,white-light,1,100]",directory=self.bdir+, x0=xAIA, y0=yAIA, width=h, height=w)
