# make small boxes

numprocs = 4

box = 16


#---------------------------------------------------------------
#---------------------------------------------------------------
import numpy as np

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

import numpy as n
import scipy.interpolate
import scipy.ndimage

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = n.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None


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
    
 
