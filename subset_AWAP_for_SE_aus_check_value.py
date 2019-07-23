#!/usr/bin/env python
"""
Subset NVIS file for SE Aus.
Original data from the netcdf Anna made
/srv/ccrc/data04/z3509830/LAI_precip_variability/GRID_NVIS4_2_AUST_EXT_MVG
That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (27.06.2019)"
__email__ = "mdekauwe@gmail.com"


import xarray as xr
import matplotlib.pyplot as plt

#fname = "/short/w35/mm3972/data/AWAP_to_netcdf/Tair/AWAP.Tair.3hr.2000.nc"
#fname = "/short/w35/mm3972/cable/src/CABLE-AUX/offline/MD_elev_orig_std_avg-sand_AWAP_AU_mask.nc"
#fname = "/short/w35/mm3972/cable/src/CABLE-AUX/offline/MD_elev_orig_std_avg-sand_AWAP_AU_landmask.nc"
#fname = "/short/w35/mm3972/data/AWAP_to_netcdf/Wind/AWAP.Wind.3hr.2000.nc"
fname = "/short/w35/mm3972/cable/src/CABLE-AUX/offline/gswp/Wind/GSWP3.BC.Wind.3hrMap.2000.nc"
#out_fname = "data/SE_aus_awap_grid.nc"

ds = xr.open_dataset(fname)

lat_bnds, lon_bnds = [-34.755,-34.745], [145.245,145.255]
#ds = ds.ssat_vec.sel(latitude=slice(*lat_bnds), longitude=slice(*lon_bnds))
#ds = ds.landsea.sel(lat=slice(*lat_bnds), lon=slice(*lon_bnds))
ds = ds.Wind.sel(lat=slice(*lat_bnds), lon=slice(*lon_bnds))
print(ds[0,:,:])

#fig = plt.figure()

#plt.imshow(ds[:,:], origin='upper')
#plt.colorbar()
#plt.show()

#fig.savefig("landsea", bbox_inches='tight', pad_inches=0.1)

#ds.to_netcdf(out_fname)
