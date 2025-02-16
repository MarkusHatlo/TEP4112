import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

data = scipy.io.loadmat('windturbine_shear_PIVdata.mat')

x1 = data['X1'].flatten()  
x3 = data['X3'].flatten()  
U1 = data['U1']  
U3 = data['U3'] 
u1_prime = data['u1_tseries']  
u3_prime = data['u3_tseries']  
U0 = 9.91  
D = 0.21   


x1_D = x1 / D
x3_D = x3 / D

print(u1_prime)

# (a) Plot U3/U0 Field
plt.figure(figsize=(8, 6))
plt.contourf(x1_D, x3_D, U3/U0, levels=50, cmap='coolwarm')
plt.colorbar(label='U3/U0')
plt.xlabel('x1/D')
plt.ylabel('x3/D')
plt.title('Normalized U3/U0 Field')
plt.show()

# (b) Calculate & Plot Local Turbulence Intensity
turb_intensity_1 = np.abs(u1_prime) / U1
turb_intensity_3 = np.abs(u3_prime) / U1

plt.figure(figsize=(8, 6))
plt.contourf(x1_D, x3_D, turb_intensity_1, levels=50, cmap='plasma')
plt.colorbar(label='Turbulence Intensity u1/U1')
plt.xlabel('x1/D')
plt.ylabel('x3/D')
plt.title('Local Turbulence Intensity (u1)')
plt.show()

plt.figure(figsize=(8, 6))
plt.contourf(x1_D, x3_D, turb_intensity_3, levels=50, cmap='plasma')
plt.colorbar(label='Turbulence Intensity u3/U1')
plt.xlabel('x1/D')
plt.ylabel('x3/D')
plt.title('Local Turbulence Intensity (u3)')
plt.show()