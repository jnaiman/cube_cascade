# this will plot how the cube of data should look

# what is the size of the box
box = 8 # box X box X box size

data_files = '/Users/jillnaiman1/illustrisData/L25n128TNG/ledcube/coords*box08*hdf5'

#------------------------------------------

import glob as glob
import h5py

from sys import path
path.append('/Users/jillnaiman1/pulsarPython/')
#path.append('/Users/jillnaiman/pulsarPython/')
from computer_selector import select_computer
sysflag, basedir, savefigbase, codedir, matplotlib, plt, \
  pylab, datadir, runsdir, runsdir_large = select_computer() # import the proper paths of things

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

plt.ion() # interactive plots

files = glob.glob(data_files)

# sort all the files
files.sort()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.linspace(0,box-1,box)
y = np.linspace(0,box-1,box)
z = np.linspace(0,box-1,box)

nn = 255
x = np.arange(nn)
ys = [i+x+(i*x)**2 for i in range(nn)]
#colors = cm.rainbow(np.linspace(0, 1, len(ys)))
colors = cm.Reds(np.linspace(0, 1, len(ys)))

for i in xrange(0,len(files)):
    f = h5py.File(files[i], 'r')
    coords = f['coords'].value
    **need to unloop colors**
    if i == 0:
        ax.scatter(x,y,z, c=colors[ind], marker='o')
        ax.set_xlabel('X Label')
        aqx.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
    else:
        ax.lines[0].set_data(x,y,z,c=colors[ind])




plt.show()
