#!/usr/bin/env python
"""
Plot PDSI for NSW using data from
http://www.cgd.ucar.edu/cas/catalog/climind/pdsi.html

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (07.04.2016)"
__email__ = "mdekauwe@gmail.com"

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import sys
import datetime as dt
import calendar
import pandas as pd
import brewer2mpl
from brewer2mpl import diverging
import brewer2mpl
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import AxesGrid
import glob
import os
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sns
sns.set(style="white")

def colorbar_index(cax=None, ncolours=None, cmap=None, orientation=None):
    cmap = cmap_discretize(cmap, ncolours)
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolours+0.5)
    colorbar = plt.colorbar(mappable, cax=cax, orientation=orientation)
    colorbar.set_ticks(np.linspace(0, ncolours, ncolours))

    return colorbar

def cmap_discretize(cmap, N):
    """
    Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)

    """
    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in range(N + 1)]
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)


"""
fname = "pdsi.monthly.maps.1870-2005.fawc=1.r2.5x2.5.nc"
f = nc.Dataset(fname, 'r')
lats = f.variables['lat'][:]
lons = f.variables['lon'][:]
time = f.variables['time'][:]
# 2000-2005
pdsi = f.variables['PDSI'][:,:,:]
"""

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.sans-serif'] = "Helvetica"

fig = plt.figure(figsize=(14, 10))
grid = AxesGrid(fig, [0.05,0.05,0.9,0.9], nrows_ncols=(10,12), axes_pad=0.1,
                cbar_mode='single', cbar_pad=0.2, cbar_size="3%",
                cbar_location='bottom', share_all=True)



path = "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/"
case1 = "GSWP3_gw_on_no_aquifer_influence_unify_para/" #"GSWP3_gw_on_recharge=0_unify_para/"
# case2 = "GSWP3_gw_on_unify_para/"


var_name  = 'Evap' # 'TVeg', 'ESoil', 'SoilMoist'
var_unit  = '(mm/d)' # '(m3/m3)'
plot_type = 'Anomaly' # 'Difference'
# Range on colourbar
tag       = 'norecharge'
ncolours = 11
vmin = -0.5
vmax = 0.5

counter = 0

for year in range(2001, 2011, 1):

    f1 = nc.Dataset('%s%sOutputs/gw_on/cable_out_%s.nc' %(path, case1, year), 'r')

    if var_name == 'SoilMoist':
        if year == 2001:
            var1_aclm = np.average(((f1.variables['SoilMoist'][:,0,:,:]*0.022+f1.variables['SoilMoist'][:,1,:,:]*0.058 \
        	      + f1.variables['SoilMoist'][:,2,:,:]*0.154+f1.variables['SoilMoist'][:,3,:,:]*0.409 \
        	      + f1.variables['SoilMoist'][:,4,:,:]*1.085+f1.variables['SoilMoist'][:,5,:,:]* 2.872)/4.6),axis = 0)
        else:
        	var1_aclm = var1_aclm + np.average(((f1.variables['SoilMoist'][:,0,:,:]*0.022+f1.variables['SoilMoist'][:,1,:,:]*0.058 \
        	      + f1.variables['SoilMoist'][:,2,:,:]*0.154+f1.variables['SoilMoist'][:,3,:,:]*0.409 \
        	      + f1.variables['SoilMoist'][:,4,:,:]*1.085+f1.variables['SoilMoist'][:,5,:,:]* 2.872)/4.6),axis = 0)
    else:
        if year == 2001:
            var1_aclm = np.average(f1.variables[var_name][:,:,:],axis = 0)*3600*24
        else:
            var1_aclm = var1_aclm + np.average(f1.variables[var_name][:,:,:],axis = 0)*3600*24

var1_avg = var1_aclm/10.

for year in range(2001, 2011, 1):

    f1 = nc.Dataset('%s%sOutputs/gw_on/cable_out_%s.nc' %(path, case1, year), 'r')
    lats = f1.variables['latitude'][:,0]
    lons = f1.variables['longitude'][0,:]
    time = f1.variables['time'][:]

    if var_name == 'SoilMoist':
        var1 = (f1.variables['SoilMoist'][:,0,:,:]*0.022+f1.variables['SoilMoist'][:,1,:,:]*0.058 \
             	+ f1.variables['SoilMoist'][:,2,:,:]*0.154+f1.variables['SoilMoist'][:,3,:,:]*0.409 \
        	    + f1.variables['SoilMoist'][:,4,:,:]*1.085+f1.variables['SoilMoist'][:,5,:,:]* 2.872)/4.6
    else:
        var1 = f1.variables[var_name][:,:,:]

    m = Basemap(projection='cyl', llcrnrlon=lons[0], llcrnrlat=lats[0], \
            urcrnrlon=lons[-1], urcrnrlat=lats[-1], resolution='c')

    for mon in range(12):
        var1[mon,:,:] = var1[mon,:,:] - var1_avg
        bmap = sns.blend_palette(["red", "white", "blue"], ncolours, as_cmap=True)
        ax = grid[counter]
        m.ax = ax
        m.drawcoastlines(linewidth=0.1, color='k')
        m.drawcountries(linewidth=0.1, color='k')
        ## Here !
        image = m.imshow(var1[mon,:,:], bmap,
                     colors.Normalize(vmin=vmin, vmax=vmax, clip=True),
                     interpolation='nearest')

        cbar = colorbar_index(cax=grid.cbar_axes[0], ncolours=ncolours, cmap=bmap,
                          orientation='horizontal')
        cbar.set_ticklabels(np.linspace(vmin, vmax, num=ncolours))
        cbar.set_label("%s-%s %s" %(var_name, plot_type, var_unit), fontsize=16)
        ax.set_xlim(140.5, 154)
        ax.set_ylim(-38, -28)

        if counter == 0:
            yr = 2001
        elif counter == 12:
            yr = 2002
        elif counter == 24:
            yr = 2003
        elif counter == 36:
            yr = 2004
        elif counter == 48:
            yr = 2005
        elif counter == 60:
            yr = 2006
        elif counter == 72:
            yr = 2007
        elif counter == 84:
            yr = 2008
        elif counter == 96:
            yr = 2009
        elif counter == 108:
            yr = 2010

        if counter == 0 or counter == 12 or counter == 24 or counter == 36 or counter == 48 or\
           counter == 60 or counter == 72 or counter == 84 or counter == 96 or counter == 108:
            textstr='%d' % (yr)
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(-0.5, 0.6, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)

        if counter == 0:
            textstr='Jan'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 1:
            textstr='Feb'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 2:
            textstr='Mar'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 3:
            textstr='Apr'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 4:
            textstr='May'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 5:
            textstr='Jun'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 6:
            textstr='Jul'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 7:
            textstr='Aug'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 8:
            textstr='Sep'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 9:
            textstr='Oct'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 10:
            textstr='Nov'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)
        elif counter == 11:
            textstr='Dec'
            props = dict(boxstyle='round', facecolor='white', alpha=1.0, ec="white")
            ax.text(0.3, 1.3, textstr, transform=ax.transAxes, fontsize=12,
                     verticalalignment='top', bbox=props)

        counter = counter + 1

fig.savefig("%s-%s-%s.png" %(var_name, plot_type, tag), bbox_inches='tight', pad_inches=0.1, dpi=300)
plt.show()
