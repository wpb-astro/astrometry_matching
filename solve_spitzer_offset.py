"""
solve for the astrometric solution to convert Spitzer --> HLF v2.0

- open and plot the image in aplpy
- do source detection in SEP
- grab the brightest N objects
- plot them on the image
- save a table of the WCS positions

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

import sep_demo

def plot_image(filename='../pointings/hlsp_xdf_hst_wfc3ir-60mas_hudf_f160w_v1_sci.fits',
               invert=True, **kwargs):
    f1 = aplpy.FITSFigure(filename)
    f1.set_title(filename.split('/')[-1])
    f1.show_grayscale(invert=invert, **kwargs)
    return f1

def get_cutout_recentering():
    # recenter the images based on the Spitzer cutout
#    return {'x':53.16854000, 'y':-27.79702000, 'radius':32/3600.}
    return {'x':53.16249000, 'y':-27.79141000,'radius':150/3600.}

def plot_and_extract(fitsfile, n=False, savetablename=False, savefigs=False):
    '''
    
    '''
    objects, vminmax_dict = sep_demo.sep_run(fitsfile)
    fig1 = plt.gcf()
    fig2 = plot_image(filename=fitsfile, invert=True, **vminmax_dict)
    fig2.recenter( **get_cutout_recentering() )
    # plot positions of the detections
    circles_kwargs = {'radius': 2/3600.}
    coords = fig2.pixel2world(objects['x'], objects['y'])
    coords = np.hstack(coords).reshape((len(coords[0]),2),order='F')
    xyrad = get_cutout_recentering()
    sel =  ( np.abs(coords[:,0] - xyrad['x']) <= xyrad['radius'])
    sel *= ( np.abs(coords[:,1] - xyrad['y']) <= xyrad['radius'])
    coords = coords[sel]
    objects = objects[sel]

    for i, (xw, yw) in enumerate(coords): #zip(*coords)):
        fig2.show_circles(xw,yw,color='y',lw=1.1, **circles_kwargs)

    if n is not False:
        objects = sep_demo.get_pixels_brightsort(objects, n=n, cpeak=True)
        coords = fig2.pixel2world(objects['x'], objects['y'])
        coords = np.hstack(coords).reshape((len(coords[0]),2),order='F')
        for i, (xw, yw) in enumerate(coords): #zip(*coords)):
            fig2.show_circles(xw,yw,color='r',lw=2., **circles_kwargs)

    t=Table(coords, names=['ra','dec'])
    if savetablename is not False: 
        t.write( savetablename, format='ascii')
    if savefigs:
        fig1.savefig('sep_fig.png',dpi=150)
        fig2.savefig('aplpy_fig.png',dpi=150)

    return t


