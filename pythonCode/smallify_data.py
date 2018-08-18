# make small boxes

numprocs = 4

box = 8


#---------------------------------------------------------------
#---------------------------------------------------------------
import numpy as np
import yt

no_data_drive = True # stuff on computer
base_directory = 'L25n128TNG'
nfiles = 50

from sys import path # all paths to look for computer selector
path.append('/n/home01/jnaiman/ArepoCodePython/')
path.append('/Users/jillnaiman1/ArepoCodePython/')
path.append('/Users/jillnaiman/ArepoCodePython/')
path.append('/n/home01/jnaiman/ArepoCodePython/rprocess/')
path.append('/Users/jillnaiman1/ArepoCodePython/rprocess/')
path.append('/Users/jillnaiman/ArepoCodePython/rprocess/')
from computer_selector import select_computer
sysflag, basedir, savefigbase, \
    codedir, matplotlib, \
    plt, pylab = \
                 select_computer(path, False, paper='mwalphas') # import the proper paths of things

from chemutils import element_info, elements_to_plot

import numpy as np
import matplotlib.cm as cm

from rebin import rebin

if no_data_drive and sysflag == 'HOME':
    basedir = '/Users/jillnaiman1/illustrisData/'

# import illustris python stuffs
import illustris_python as il
from haloutils import grab_snapnum, grab_mwsubhalos, list_full_snapshots
from mpl_toolkits.axes_grid1 import make_axes_locatable

outputdir = '/Output'
# just to check
if sysflag == 'ODYSSEY' or sysflag == 'DRACO':
    outputdir = '/output'


#---------------------------------------------------------------
#---------------------------------------------------------------

# grab hids and whatnot
snapnum = grab_snapnum(basedir, base_directory, outputdir, no_data_drive=no_data_drive)
header = il.groupcat.loadHeader(basedir+base_directory+outputdir+'/',snapnum)
full_snaps = list_full_snapshots(base_directory)

my_storage = {}
ns_array = np.linspace(0, snapnum, snapnum+1, dtype=int)
    
import h5py # gotta do this by hand... sigh.

from scipy.signal import resample_poly
# resample_poly(x, up, down, axis=0, window=('kaiser', 5.0))

nbox = np.linspace(0,box,box+1,dtype=int)

for n in ns_array:
    print(' n = ' + str(n))
    dm = il.snapshot.loadSubset(basedir+base_directory+outputdir+'/',n,'dm',['Coordinates','Velocities'])
    coords = dm['Coordinates']
    #vel = dm['Velocities']
    #vel = (vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2)**0.5

    coords[:,0] *= box/coords[:,0].max()
    coords[:,1] *= box/coords[:,1].max()
    coords[:,2] *= box/coords[:,2].max()

    xx = np.histogramdd(coords,[nbox,nbox,nbox])

    coords_new = xx[0]

    fname = basedir + base_directory + '/ledcube/coords_' \
            + str(n).zfill(3) + '_box' + str(box).zfill(2) + '.hdf5'

    f = h5py.File(fname, 'w')
    dset1 = f.create_dataset("coords", data=coords_new, dtype='f4')
    f.flush()
    f.close()
    
 
