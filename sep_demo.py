"""
example using SEP
follows demo from SEP readthedocs page

"""

import numpy as np
import sep
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Ellipse

rcParams['figure.figsize'] = [10., 8.]

fitsfile = '../pointings/spitzer/Spitzer_IRAC_DR3/0001_53.16854000_-27.79702000_s_irac_1_s1_v0.30_sci.fits'
#fitsfile = '../pointings/v1/hlsp_xdf_hst_wfc3ir-60mas_hudf_f160w_v1_sci.fits'

def sep_run(fitsfile, add_circles=False):
    '''
    get the detection object
    save the vmin, vmax params
    '''
    data=fits.open(fitsfile)[0].data
    # show the image
    m, s = np.mean(data), np.std(data)
    vminmax_dict = {'vmin':m-s, 'vmax':m+s}

    if False:
        plt.imshow(data, interpolation='nearest', cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
        plt.colorbar()

    # measure a spatially varying background on the image
    try:
        bkg = sep.Background(data)
    except ValueError:
        bkg = sep.Background( data.byteswap().newbyteorder() )

    # evaluate background as 2-d array, same size as original image
    bkg_image = np.array(bkg)
    if False:
        plt.figure()
        plt.imshow(bkg_image, interpolation='nearest', cmap='gray', origin='lower')
        plt.colorbar()

    # subtract the background
    data_sub = data - bkg

    # object detection - threshold = 1.5 * err
    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)

    # plot
    from matplotlib.patches import Ellipse

    # plot background-subtracted image
    if True:
        fig, ax = plt.subplots()
        m, s = np.mean(data_sub), np.std(data_sub)
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                       vmin=m-s, vmax=m+s, origin='lower')
    #    plt.colorbar(cax=ax)

    if add_circles:
        # plot an ellipse for each object
        for i in range(len(objects)):
            e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                        width=3*objects['a'][i],
                        height=3*objects['b'][i],
                        angle=objects['theta'][i] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)

    # apphot - 3 pixel radius
    flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'],
                                         r=3.0, err=bkg.globalrms, gain=1.0)

    return objects, vminmax_dict

def get_pixels_brightsort(objects, n=20, peak=False, cpeak=False):
    if peak is cpeak:
        print('identify the column')
        return
    if peak:
        col = 'peak'
    if cpeak:
        col = 'cpeak'
    o = objects[ np.argsort(objects[col])][::-1]
    if n < len(o):
        o = o[:n]
    return o 

