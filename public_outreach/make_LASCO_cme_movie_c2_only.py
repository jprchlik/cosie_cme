import matplotlib
matplotlib.use('Agg')
from skimage import io, transform
import astropy.units as u
import matplotlib as mpl
#Noise filter
import scipy.signal as signal 
#sharpening
import scipy.ndimage as ndimage
from skimage.restoration import denoise_tv_chambolle, denoise_bilateral

from skimage.filters.rank import enhance_contrast
from multiprocessing import Pool
from skimage.filters.rank import mean_bilateral
from skimage.morphology import disk,reconstruction
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import glob
import pyfits
import sunpy.cm as cm
import sunpy.map as map
import sunpy.sun as sun


def update_movie(num):
   try:
#read fits files into maps
        ifile = readc2[num]
        mymap = map.Map(ifile)
#        olmap = map.Map(readc3[j-1],readc2[j-1],composite=True)
#make sure there are both c2 and c3 maps in the map file
#        print '###########################'
#        print 'frame={0:4d}, file={1}'.format(num,ifile) 
#        print '###########################'

#do a first reduction of the image
        mymap,img2,l2e,c2c,goodl2,badl2 = reduce_mymap(mymap,ifile,diag=True,first=True)
#actually plot the image
        dpi = 100
        fig, ax = plt.subplots(figsize=(float(sc_w2)/dpi,float(sc_h2)/dpi))
        fig.set_dpi(dpi)
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
        ax.set_axis_off()
# fix center from being white
        c2c = new_fix_white(c2c)
        im2 = ax.imshow(c2c,extent=l2e,origin='lower',zorder=1)
#savefigure to file
        fig.savefig('data/lasco_new_power/frame_{0:5d}.png'.format(num).replace(' ','0'),bbox_inches=0)
        plt.clf()
        plt.close()
        return 
   except:
        plt.clf()
        plt.close()
        return 
    


#smooths over bad pixels
def fix_valleys(img):
    seed = np.copy(img)
    seed[1:-1,1:-1] = img.max()
    mask = img
    fill = reconstruction(seed,mask,method='erosion')
    return fill
#replace white space in image
def fix_peaks(img):
    seed[1:-1,1:-1] = img.max()
    mask = img
    fill = reconstruction(seed,mask,method='erosion')
    return fill




#most replace white (hot) pixels with black pixels
def new_fix_white(c2c):

    wc2bad = c2c[:,:,3] == 0
#    wc2bad = ((wc2bad[:,:,0]) & (wc2bad[:,:,1]) & (wc2bad[:,:,2]) & (wc2bad[:,:,3]))  
    c2c[wc2bad,0] = 0.0
    c2c[wc2bad,1] = 0.0
    c2c[wc2bad,2] = 0.0
    c2c[wc2bad,3] = 1.0

    return c2c


#put circles around LASCO data
def reduce_mymap(mymap,ifile,diag=False,first=False):
# just reuse the good and bad pixel values from the first image
    if first == False:
        global badl2,goodl2


#get the solar radius in arceconds
#    ttime = ifile[10:14]+'-'+ifile[14:16]+'-'+ifile[16:18]+' '+ifile[19:21]+':'+ifile[21:23]+':'+ifile[23:25]
    mymap.meta['date-obs'] = mymap.meta['date-obs'][:-1].replace('T',' ')
    
    rsun_arcseconds2 = sun.solar_semidiameter_angular_size(t=mymap.date).value

    img2  = mymap

#    remove_bad(mymap,c2,c3,img2,img3)

    #extent of image in arcseconds
#j   l2e = img2.xrange.tolist()+img2.yrange.tolist()
    l2e = [img2.xrange[0].value,img2.xrange[1].value,img2.yrange[0].value,img2.yrange[1].value]


    
#create array of the maximum and minimum radii for lasco2 and lasco3
    pltr2 = np.array([r2min,r2max])
    
    #find extent in solar radii
    s2e = np.array(l2e)/rsun_arcseconds2

    
    
    if first:
        #get data coordinates in x(arcsec), y(arcsec)
        mesh_x, mesh_y = np.meshgrid(np.arange(mymap.data.shape[0]),np.arange(mymap.data.shape[1]))
        mesh_x2, mesh_y2 = mymap.pixel_to_data(mesh_x.T*u.pix,mesh_y.T*u.pix)
    
       #covert arcseconds to solar radii and find the radial distance  
        r2 = np.sqrt((mesh_x2)**2+(mesh_y2)**2)/rsun_arcseconds2
        
        #find the between a filtered solar radii 
        goodl2 = ((r2 > r2min*u.arcsec) & (r2 < r2max*u.arcsec))
        
    	#make mask for the data to create black strips 
        badl2 = goodl2 == False

#use image2 and img3 as pixel arrays
    img2 = img2.data

#fix low.dead pixels
    img2 = fix_valleys(img2)

# Define color image normalization
    bot2,top2 = 2.5,100
    fit2 = img2 > 0
    power2 = .25


    timg2 = np.log10(img2)
    timg2 = np.arcsinh(img2)

# get percentile values for lasco2 and lasco3 images (*1e10 to make number greater than 1 for convience)
#    v2_min,v2_max = -2.7972763e-10,3.9886609e-9
    v2_min,v2_max = np.percentile((timg2[fit2]),(bot2,top2))
# hand set for stability
#    v2_max = np.log10(1.8e-7)
    v2_max = np.arcsinh(1.8e-7)
#enhance contrast
#create color normalization object for c2 data
    normc2 = mpl.colors.Normalize(vmin=v2_min,vmax=v2_max+0.4)#*(v2_max-v2_min))
    c2c = cm.cm.soholasco2(normc2(timg2))
#set alpha to zero so you cannot see the l2 data
    c2c[:,:,3][badl2] = 0.0

    if first:
	    return mymap,img2,l2e,c2c,goodl2,badl2
	#return parameters for plotting
    if diag:
        return mymap,img2,l2e,c2c
    else:
        return mymap

#grad all LASCO files
startd = '/Volumes/Pegasus/jprchlik/LASCO/'
startd = 'level1/'
c2 = glob.glob(startd+'C2*fts')

#replace all characters to make a time date array
c2d = np.array([float(w.replace('_','').replace('C2','').replace('.fts','').replace(startd,'')) for w in c2])

#sort c2 data by time
idx = np.argsort(c2d)
c2 = np.array(c2)[idx]
c2d = c2d[idx]


dat2 = pyfits.open(c2[0])

img2 = dat2[0].data

#####find new size
sc_w2 = np.shape(img2)[0]
sc_h2 = np.shape(img2)[1]
#solar radii to cut
r2min = 2.4
r2max = 6.0
r3min = 5.0
r3max = 128.0

readc2 = c2



#num = 0
#readc3 = readc3[:]
ain = np.arange(2020)+3000
ain = np.arange(10)#3000
pool = Pool(processes=10)
out = pool.map(update_movie,ain)

#for j,ifile in enumerate(readc3):
#for j in ain:
#    num = update_movie(readc2[j],num)
#
