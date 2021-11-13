"""
Simple script to download etopo1 topography/bathymetry data from
    http://www.ngdc.noaa.gov/mgg/global/global.html

The etopo1 data has 1-arcminute resolution, but you can request coarsening.
E.g. set resolution = 4./60. for 4-arcminute resolution.

"""

import os
from clawpack.geoclaw import topotools

plot_topo = True

# Set the limits of the domain and resolution:

resolution = 1./60.   # degrees

# If user environment variable ETOPO is set to a valid path to a directory, 
# downloaded data will be stored there.  Allows sharing same data among
# various projects.

try:
    etopo_dir = os.environ['ETOPO']
    os.chdir(etopo_dir)  # make sure it's a valid directory
except:
    print("ETOPO not set or invalid directory, setting etopo_dir='.'")
    etopo_dir = '.'

extent = [-125,-123,46,48]
coarsen = 1
topo = topotools.read_netcdf('etopo1', extent=extent, 
                             coarsen=coarsen, verbose=True)

name = 'etopo1_%s_%s_%s_%s_%smin' % tuple(extent + [coarsen])
print('name = ',name)
fname = os.path.join(etopo_dir, name + '.asc')
topo.write(fname, topo_type=3, header_style='asc', 
                     grid_registration='llcorner', Z_format='%.0f')



if plot_topo:
    # plot the topo and save as a png file...
    import matplotlib.pyplot as plt
    topo.plot()
    plt.title('Topo file %s' % name)
    fname = name + '.png'
    fname = os.path.join(etopo_dir, fname)
    plt.savefig(fname)
    print('Created %s' % fname)


