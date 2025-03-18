import scipy as sp
import numpy as np

boundarylayer_udata = np.loadtxt("data/boundarylayer_udata.txt", skiprows=1)
turbulencegrid_udata = np.loadtxt("data/turbulencegrid_udata.txt", skiprows=1)

boundrylayer_t = boundarylayer_udata[:,0]
boundrylayer_u = boundarylayer_udata[:,1]

turbulencegrid_t = turbulencegrid_udata[:,0]
turbulencegrid_u = turbulencegrid_udata[:,1]

boundrylayer_u = np.mean(boundrylayer_u)
turbulencegrid_U = np.mean(turbulencegrid_u)
