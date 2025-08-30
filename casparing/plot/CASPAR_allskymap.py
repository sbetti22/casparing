from matplotlib import patheffects
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes.frame import EllipticalFrame

import .CASPAR_plotting as cplot
import warnings
warnings.filterwarnings("ignore")


def insets(ax, a, wcs, ax_pos, bounds, manual_set=False):
    axins = ax.inset_axes(ax_pos)
    cmap = cm.binary_r
    if manual_set:
        new_array = np.ma.array(a, mask=a==0)
        cmap.set_bad('white',1.)
    else:
        new_array = a
    axins.imshow(new_array, vmin=-6.1, vmax=11., origin='lower', cmap=cmap, aspect='auto')
    x1, x2, y1, y2 = bounds
    ls = 5
    if manual_set:
        sky = SkyCoord(b=[-2, -16]*u.degree,l=[175,170]*u.degree, frame='galactic')
        x,y = wcs.world_to_pixel(sky)
        axins.set_xticks(x)
        axins.set_yticks(y)
        axins.set_yticklabels(['$-$2', '$-$16'])
        axins.set_xticklabels([175, 170])
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.minorticks_on()
        axins.tick_params(which='major', labelsize=ls, length=2, direction='in', top=True, right=True)
    else:
        axinsY = axins.twinx().twiny()
        axinsX = axins.twinx()
        gl, gb = wcs.wcs_pix2world([x1,x2], [y1, y2], 1)#, [x2, y2])
        print('gl', gl, 'gb', gb)
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axinsY.set_xlim(gl[0], gl[1])
        axinsY.set_ylim(gb[0], gb[1])
        axins.set_xticks([])
        axins.set_yticks([])
        axinsX.set_ylim(gb[0], gb[1])
        axinsY.xaxis.set_ticks_position("bottom")
        axinsY.set_yticks([])
        axinsX.yaxis.set_ticks_position("left")
        axinsX.tick_params(which='major', labelsize=ls, length=2, direction='in', top=True, right=True)
        axinsY.tick_params(which='major', labelsize=ls, length=2, direction='in', top=True, right=True)
    return axins
    
def allsky(df_caspar, akari_file, **kwargs):
    '''kwargs: savefig'''
#    a = fits.getdata('../akari_mollweide_65_1_4096.fits')
#    hdu = fits.open('../akari_mollweide_65_1_4096.fits')[1]
    a = fits.getdata(akari_file)
    hdu = fits.open(akari_file)[1]
    wcs = WCS(hdu.header)
    plt.figure(dpi=150)
    ax = plt.subplot(projection=wcs, frame_class=EllipticalFrame)

    path_effects=[patheffects.withStroke(linewidth=3, foreground='black')]
    ax.coords.grid(color='white')
    ax.coords['glon'].set_ticklabel(color='white', path_effects=path_effects)

    im = ax.imshow(a, vmin=-6.1, vmax=11., origin='lower', cmap='binary_r')

    # Clip the image to the frame
    im.set_clip_path(ax.coords.frame.patch)
    
    axins1 = insets(ax, a, wcs, [-0.03, -0.6, 0.3, 0.4], [10, 200,650, 1000], manual_set=True) # Taurus
    axins2 = insets(ax, a, wcs, [0.37, -0.6, 0.3, 0.4], [2777,2788,716,728]) # n Cha 
    axins3 = insets(ax, a, wcs, [0.75, -0.6, 0.3, 0.4], [3720, 3780, 750, 810]) #sigma ori
    axins4 = insets(ax, a, wcs, [-0.03, 1.2, 0.3, 0.4], [1973, 1983, 995, 1015]) # Lagoon 
    axins5 = insets(ax, a, wcs, [0.37, 1.2, 0.3, 0.4], [2240,2350,1090,1350]) # Lupus
    axins6 = insets(ax, a, wcs, [0.75, 1.2, 0.3, 0.4], [2739,2760,796,822]) # Cha I
 
    cluster, marker= cplot.plot_cluster()

    df_caspar['Association'] = df_caspar['Association'].fillna('none')
    color = cm.jet(np.linspace(0.2, 1, len(cluster)-1))
    for i, name in enumerate(cluster):
        df_cluster = df_caspar[df_caspar['Association'].str.contains(name)]
        if len(df_cluster)>0:
            c = SkyCoord(ra=df_cluster['RA (J2000.0)']*u.degree, 
                         dec=df_cluster['Dec (J2000.0)']*u.degree, frame='icrs').galactic
            ra, dec = c.l, c.b
            if name == 'p Oph':
                name = 'ρ Oph'
            if name == 'n Chamaeleontis':
                name = 'η Chamaeleontis'
            if name == 'B Pictoris':
                name = 'β Pictoris'
            if name == 'Field':
                name = 'Field/unknown'
                col='k'
            else:
                col=color[i]
            x, y = wcs.world_to_pixel(c)
            ax.scatter(ra, dec, facecolor=col, zorder=10, s=25, edgecolor='k',marker=marker[i], 
                      label=name, alpha=0.7, transform=ax.get_transform('world'))

            for axi in [axins1, axins2, axins3, axins4, axins5,axins6]:
                axi.scatter(x,y, facecolor=col, zorder=1000, s=20, edgecolor='k',marker=marker[i], 
                  alpha=.7)


    leg = ax.legend(loc='upper left', bbox_to_anchor=(1.1, 1.05), fontsize=4)
    ax.tick_params(which='major', labelsize=7)
    ax.indicate_inset_zoom(axins1, edgecolor="black")
    ax.indicate_inset_zoom(axins2, edgecolor="black")
    ax.indicate_inset_zoom(axins3, edgecolor="black")
    ax.indicate_inset_zoom(axins4, edgecolor="black")
    ax.indicate_inset_zoom(axins5, edgecolor="black")
    ax.indicate_inset_zoom(axins6, edgecolor="black")
    for axi in [axins1, axins2, axins3, axins4, axins5,axins6]:
        axi.minorticks_on()
        axi.tick_params(which='both', direction='in', top=True, right=True)
    plt.tight_layout()
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.show()
    
