import scipy as sp
import numpy as np
import scipy.stats as sps
from scipy.stats import norm
import matplotlib.pyplot as plt

boundaryLayer_udata = np.loadtxt("data/boundarylayer_udata.txt", skiprows=1)
turbulenceGrid_udata = np.loadtxt("data/turbulencegrid_udata.txt", skiprows=1)

boundryLayer_t = boundaryLayer_udata[:,0]
boundryLayer_u = boundaryLayer_udata[:,1]

turbulenceGrid_t = turbulenceGrid_udata[:,0]
turbulenceGrid_u = turbulenceGrid_udata[:,1]

boundryLayer_U = np.mean(boundryLayer_u)
turbulenceGrid_U = np.mean(turbulenceGrid_u)

boundryLayer_u_prime = boundryLayer_u - boundryLayer_U
turbulenceGrid_u_prime = turbulenceGrid_u - turbulenceGrid_U

boundryLayer_std = np.std(boundryLayer_u_prime)
turbulenceGrid_std = np.std(turbulenceGrid_u_prime)

#Task a

# boundryLayer_turbulence_intensity = boundryLayer_u_prime/boundryLayer_U
# turbulenceGrid_turbulence_intensity = turbulenceGrid_u_prime/turbulenceGrid_U

boundryLayer_turbulence_intensity = boundryLayer_std/boundryLayer_U
turbulenceGrid_turbulence_intensity = turbulenceGrid_std/turbulenceGrid_U

print(f"This is the turbulence intensity for both data sets {np.mean(boundryLayer_turbulence_intensity):.3g} and {np.mean(turbulenceGrid_turbulence_intensity):.3g}\n")
print
#Task b

boundryLayer_skewness = sps.skew(boundryLayer_u)
turbulenceGrid_skewness = sps.skew(turbulenceGrid_u)

print(f"This is the skewness for both data sets {boundryLayer_skewness:.3g} and {turbulenceGrid_skewness:.3g}\n")

#Task c
boundryLayer_kurtosis= sps.kurtosis(boundryLayer_u)
turbulenceGrid_kurtosis= sps.kurtosis(turbulenceGrid_u)

print(f"This is the kurtosis for both data sets {boundryLayer_kurtosis:.3g} and {turbulenceGrid_kurtosis:.3g}")

#Task d

# boundryLayer_probability = norm.cdf(boundryLayer_u ,loc= boundryLayer_U , scale= np.std(boundryLayer_u))
# boundryLayer_x_axis = np.linspace(0,1,len(boundryLayer_probability))

# plt.figure(figsize=(10, 6))
# plt.plot(boundryLayer_probability,boundryLayer_x_axis)
# plt.show()


boundryLayer_normalized = boundryLayer_u_prime / boundryLayer_std
turbulenceGrid_normalized = turbulenceGrid_u_prime / turbulenceGrid_std

x_range = np.linspace(-4, 4, 1000)
gaussian_pdf = sps.norm.pdf(x_range)

boundryLayer_kde = sps.gaussian_kde(boundryLayer_normalized)
turbulenceGrid_kde = sps.gaussian_kde(turbulenceGrid_normalized)

plt.figure(figsize=(10, 6))

plt.plot(x_range, boundryLayer_kde(x_range), 'b-', linewidth=2, label='Boundary Layer')
plt.plot(x_range, turbulenceGrid_kde(x_range), 'r-', linewidth=2, label='Turbulence Grid')
plt.plot(x_range, gaussian_pdf, 'k--', linewidth=2, label='Gaussian (S=0, K=3)')

plt.xlabel(r'$u\' / \sigma$')
plt.ylabel('Probability Density')
plt.title('PDF of Normalized Velocity Fluctuations')
plt.legend()
plt.grid(True)
plt.xlim(-4, 4) 
plt.show()

plt.figure(figsize=(10, 6))

plt.semilogx(x_range, boundryLayer_kde(x_range), 'b-', linewidth=2, label='Boundary Layer')
plt.semilogx(x_range, turbulenceGrid_kde(x_range), 'r-', linewidth=2, label='Turbulence Grid')
plt.semilogx(x_range, gaussian_pdf, 'k--', linewidth=2, label='Gaussian (S=0, K=3)')

plt.xlabel(r'$u\' / \sigma$')
plt.ylabel('Probability Density')
plt.title('PDF of Normalized Velocity Fluctuations in semilog plot')
plt.legend()
plt.grid(True) 
plt.show()