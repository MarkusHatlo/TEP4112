import numpy as np
import scipy.io as sp
import matplotlib.pyplot as plt

data = sp.loadmat("Xwire_ui")
u = np.float64(data['u'].squeeze())
v = np.float64(data['v'].squeeze())
sampling_frequency = 62500

#a)

u_prime = u - np.mean(u)
v_prime = v - np.mean(v)

u_rms = np.std(u_prime)
v_rms = np.std(v_prime)
ratio = u_rms / v_rms
print("Ratio u'_1/u'_2 =", ratio)

#b)

#When there is isotropy in all directions all the rms values should equal each other. Meaning I can just use the 1/2u_rms 3 times
k_full = 1.5 * (u_rms**2)


#When there is only isotropy between x2 and x3 v_rms adn w_rms should equal each other. Therefore i use 1/2*v_rms 2 times
k_plane = 0.5 * (u_rms**2 + 2 * (v_rms**2))

print("TKE (full isotropy) =", k_full)
print("TKE (x2-x3 isotropy) =", k_plane)

percentage_difference = np.abs(k_full - k_plane) / k_full * 100
print("Percentage difference between TKE estimates: {:.2f}%".format(percentage_difference))


#c)

u_norm = u_prime / np.std(u_prime)
v_norm = v_prime / np.std(v_prime)

H, xedges, yedges = np.histogram2d(u_norm, v_norm, 50, density=True)

xcenters = (xedges[:-1] + xedges[1:]) / 2
ycenters = (yedges[:-1] + yedges[1:]) / 2
X, Y = np.meshgrid(xcenters, ycenters)

plt.figure(figsize=(8, 8))

contour = plt.contourf(X, Y, H.T, levels=20, cmap='viridis')
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)


plt.xlabel("Normalized $u'_1$", fontsize=14)
plt.ylabel("Normalized $u'_2$", fontsize=14)
plt.title("Contour Plot of Normalized Velocity Fluctuations", fontsize=16)
plt.colorbar(contour, label='Probability Density')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlim(-3, 3)
plt.ylim(-3, 3)

plt.show()

corr_coef = np.corrcoef(u_prime, v_prime)[0, 1]
print("Correlation coefficient between u'_1 and u'_2 =", corr_coef)

#d)

nu = 1.556e-5

U_mean = np.mean(u)

# Time step
dt = 1.0 / sampling_frequency

du_dt = np.gradient(u_prime, dt)  
dv_dt = np.gradient(v_prime, dt)

du_dx = du_dt / U_mean
dv_dx = dv_dt / U_mean

mean_du_dx_sq = np.mean(du_dx**2)
mean_dv_dx_sq = np.mean(dv_dx**2)


epsilon_1 = 15.0 * nu * mean_du_dx_sq


epsilon_2 = 3.0 * nu * (mean_du_dx_sq + 2.0 * mean_dv_dx_sq)


percent_diff = (abs(epsilon_1 - epsilon_2) / epsilon_1) * 100.0

print(f"Epsilon (full 3D isotropy): {epsilon_1:.6e} m^2/s^3")
print(f"Epsilon (x2-x3 isotropy) : {epsilon_2:.6e} m^2/s^3")
print(f"Percentage difference    : {percent_diff:.2f}%")
