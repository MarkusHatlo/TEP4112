import scipy.io
import numpy as np
import scipy.signal as sg
import matplotlib.pyplot as plt
import statsmodels.api as sm

# Load the .mat file
file_path = "XW_Sample_Sq39_farfield(1).mat"
dataValues = scipy.io.loadmat(file_path)

# Extract required variables
u = np.float64(dataValues['u'].squeeze())
nu = (dataValues['nu']).item()
fs = (dataValues['fs']).item()

#a

fc = fs/2.5
fc_new = 0
LPF_order = 7
k = 0
difference = 1


def calcDissipationRate(u,u_prime,U_mean):
    dt = 1.0 / fs

    du_dt = np.gradient(u_prime, dt)  

    du_dx = du_dt / U_mean

    mean_du_dx_sq = np.mean(du_dx**2)

    return 15.0 * nu * mean_du_dx_sq

#print(fc)
while difference > 0.01 or k == 9:
    b, a = sg.butter(LPF_order, fc*2/fs, btype='lowpass', output='ba')
    u_filt = sg.filtfilt(b,a,u)
    U_filt = np.mean(u_filt)
    up_filt = u_filt - U_filt
    dissipation = calcDissipationRate(u_filt,up_filt,U_filt)
    eta = (nu**3/dissipation)**(1/4)
    kolmogorov_frequency = U_filt/(2*np.pi*eta)
    fc_new = 1.1*kolmogorov_frequency
    difference = abs(fc-fc_new)
    fc = fc_new
    k += 1

#print(round(fc))

#b
BIN = 14

Freq, E_f = sg.welch(up_filt, fs, nperseg=2**(BIN+1))


E_k = E_f*U_filt/(2*np.pi)
k_space = 2*np.pi*fc



y_aksen = E_k*eta/nu**2
x_aksen = k_space*Freq


# Part 4d - Calculate integral length scale L11
maxLags = 40000
acf_result = sm.tsa.acf(up_filt, nlags=maxLags, fft=True)
Buu = acf_result
lags = np.arange(len(Buu))
rs = lags * (U_filt/fs)

# First, make a regular plot to identify the zero crossing
plt.figure(figsize=(10, 6))
plt.plot(rs, Buu)
plt.ylabel('Buu')
plt.xlabel('rs [m]')
plt.axhline(y=0, color='r', linestyle='--')
plt.title('Autocorrelation of u′ (linear scale)')
plt.grid(True)
plt.show()

# Now make a log scale plot similar to Figure 1 in the assignment
# (avoiding log(0) by starting from index 1)
plt.figure(figsize=(10, 6))
plt.semilogx(rs[1:], Buu[1:])  # Start from index 1 to avoid log(0)
plt.ylabel('Buu')
plt.xlabel('rs [m]')
plt.axhline(y=0, color='r', linestyle='--')
plt.title('Autocorrelation of u′ (log scale on x-axis)')
plt.grid(True)
plt.show()

# Find the first zero-crossing
zero_crossing_indices = np.where(Buu <= 0)[0]
if len(zero_crossing_indices) > 0:
    first_zero_idx = zero_crossing_indices[0]
    print(f"First zero crossing at index {first_zero_idx}, rs = {rs[first_zero_idx]} m")
else:
    first_zero_idx = len(Buu) - 1
    print("No zero crossing found, using entire range")

# Integrate from rs = 0 up to the first zero-crossing
# Use the actual rs values for integration
L11 = np.trapz(Buu[:first_zero_idx+1], rs[:first_zero_idx+1])
print(f"Integral length scale L11 = {L11:.6f} m")

# Calculate the number of decades between L11 and eta
decades = np.log10(L11) - np.log10(eta)
print(f"Number of decades between L11 and η: {decades:.2f}")
print(f"L11 = {L11:.6f} m")
print(f"η = {eta:.10f} m")