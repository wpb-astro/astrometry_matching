import numpy as np
import os
import aplpy
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import LambdaCDM
import astropy.units as u
from astropy.coordinates import SkyCoord

from basedir import basedir

import seaborn as sns
sns.set_context("talk") # options include: talk, poster, paper
sns.set_style("ticks") #whitegrid")
sns.set_style({"xtick.direction": "in","ytick.direction": "in",
               "xtick.top":True, "ytick.right":True,
               "xtick.major.size":12, "xtick.minor.size":4,
               "ytick.major.size":12, "ytick.minor.size":4,
               })
### color palettes
colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
colors += ["cloudy blue", "browny orange", "dark sea green"]
sns.set_palette(sns.xkcd_palette(colors))
orig_palette = sns.color_palette()

def plot_image(filename='../pointings/hlsp_xdf_hst_wfc3ir-60mas_hudf_f160w_v1_sci.fits',
               invert=True, **kwargs):
    f1 = aplpy.FITSFigure(filename)
    f1.set_title(filename.split('/')[-1])
    f1.show_grayscale(invert=invert, **kwargs)
    return f1

def plot_rgb_image(filename='../pointings/xdf_814_105_140_rgb-2.png'):
    f1 = aplpy.FITSFigure(filename)
    f1.show_rgb()
    return f1

def add_all_fpho_positions(aplpy_FITSFigure_obj, add_circle=0.5,
                           show_id=True, rgb_image=False,
                           save_all_images=False, save_radius=2.,
                           override_objlist=None):
    """
    plot positions of forcepho objects on an existing image
    :param add_circle: plotting radius
    :param show_id: if True, print the forcepho id on the image
    :param rgb_image: governs colors of text (show_id)
    :param save_all_images: recenter and save around all sources?
    :param save_radius: FOV radius, if saving
    :param override_objlist: list of forcepho ids (must match fpho_id)
    """
    # kwargs for the fpho sources
    circle_kwargs = {'linewidth':1.5, 'color':'firebrick'}
    text_kwargs = {'size':12, 'color':'c', 'verticalalignment':'top'}
    if not rgb_image:
        text_kwargs['color'] = 'c'
    else:
        text_kwargs['color'] = orig_palette[2]

    f1 = aplpy_FITSFigure_obj

#    t = Table.read('%scatalogs/%s' % (basedir, 'fpho_positions.dat'), format='ascii' )
    t = Table.read('../catalogs/fpho_positions.dat',format='ascii')
    ras, decs, pids = t['ra'], t['dec'], t['fpho_id']
    if override_objlist is None:
        override_objlist = pids

    for i, (ra,dec,pid) in enumerate(zip(ras, decs, pids)):
        f1.show_circles(ra, dec, radius=add_circle/3600., **circle_kwargs)
        if show_id:
            f1.add_label(ra, dec+.99*add_circle/3600., pid, **text_kwargs)
    if save_all_images:
        for i, (ra,dec,pid) in enumerate(zip(ras, decs, pids)):
            if pid in override_objlist:
                f1.recenter(ra, dec, radius=save_radius/3600.)
                plt.tight_layout()
                plt.savefig('direct_%s.png' % pid, dpi=150)


def add_specz_positions(aplpy_FITSFigure_obj, add_circle=0.5,
                        show_id=True, add_muse_pos=False, add_muse_z=False,
                        rgb_image=False):
    '''
    add positions and circles of objects in XDF with spec-z's

    :param aplpy_FITSFigure_obj:
        aplpy fitsfigure object
    :param add_circle:
        radius of circle (arcseconds)
    :param show_id:
        bool; show forcepho patch_sourceid (zero-index)
    :param add_muse_pos:
        bool; add circles around MUSE sources that are matched to fpho source
    :param add_muse_z:
        bool; add label showing the MUSE redshift
    :param rgb_image:
        bool; if plotting on color image, change color of muse labelS
    '''
    # kwargs for the fpho sources
    circle_kwargs = {'linewidth':1.5, 'color':'firebrick'}
    text_kwargs = {'size':12, 'color':'c', 'verticalalignment':'top'}

    # kwargs for the muse counterparts
    muse_circle_kwargs = {'linewidth':2, 'radius':0.2/3600.}
    muse_text_kwargs = {'size':12, 'verticalalignment':'top'}
    if not rgb_image:
        text_kwargs['color']='c'
        muse_circle_kwargs['color'] = 'b'
    else:
        text_kwargs['color']=orig_palette[2]
        muse_circle_kwargs['color'] = orig_palette[1]
    muse_text_kwargs['color'] = muse_circle_kwargs['color'] 
    muse_text_fmt = r'$z=$%0.1f'

    f1 = aplpy_FITSFigure_obj

#    t = Table.read('%scatalogs/%s' % (basedir, 'k3d_xdf.fits') )
#    t = Table.read('%scatalogs/%s' % (basedir, 'muse_overlap_pos_0.5arcsec.dat'), format='ascii' )
    t = Table.read('../catalogs/%s' % ('muse_overlap_pos_0.5arcsec.dat'), format='ascii' )
    ras, decs, pids = t['ra'], t['dec'], t['fpho_id']
    muse_ra, muse_dec, muse_z = t['ra_muse'], t['dec_muse'], t['z_muse']

    for i, (ra,dec,pid) in enumerate(zip(ras, decs, pids)):
        f1.show_circles(ra, dec, radius=add_circle/3600., **circle_kwargs)
        if show_id:
            f1.add_label(ra, dec+.99*add_circle/3600., pid, **text_kwargs)
        if add_muse_pos:
            f1.show_circles(muse_ra[i], muse_dec[i], **muse_circle_kwargs)
        if add_muse_z:
            lbl = muse_text_fmt % muse_z[i]
            f1.add_label(muse_ra[i], muse_dec[i]-1.01*.2/3600., lbl, **muse_text_kwargs)

def add_muse_sources(aplpy_FITSFigure_obj, add_circle=.2, rgb_image=False):
    '''
    add positions and circles of objects in XDF with MUSE spec-z's
    MUSE fiber size: 0.2 arcseconds

    '''
    circle_kwargs = {'linewidth':2, 'color':'b'}
    if rgb_image:
        circle_kwargs['color'] = orig_palette[1]

    f1 = aplpy_FITSFigure_obj
#    t = Table.read('%scatalogs/%s' % (basedir, 'muse_sources.dat'), format='ascii' )
    t = Table.read('../catalogs/muse_sources.dat',format='ascii')
    ras, decs = t['ra'], t['dec']

    for ra,dec in zip(ras, decs):
        f1.show_circles(ra, dec, radius=add_circle/3600., **circle_kwargs)

def save_cutouts(only_specz=True, rgb_image=True):
    f1 = plot_rgb_image()
    plt.gcf().set_size_inches(7,6,forward=True)
    add_muse_sources(f1, rgb_image=rgb_image)
    add_specz_positions(f1, add_circle=0.,
                        show_id=False, add_muse_pos=True, 
                        add_muse_z=True, rgb_image=rgb_image)
    if only_specz:
#        t = Table.read('%scatalogs/%s' % (basedir, 'muse_overlap_pos_0.5arcsec.dat'), format='ascii' )
        t = Table.read('../catalogs/%s' % ('muse_overlap_pos_0.5arcsec.dat'), format='ascii' )
        objlist = t['fpho_id']
    else:
        objlist = None

    add_all_fpho_positions(f1, add_circle=0.5,
                           show_id=True, rgb_image=rgb_image,
                           save_all_images=True, save_radius=2.,
                           override_objlist=objlist)
    return f1

