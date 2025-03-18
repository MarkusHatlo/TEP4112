import scipy.io
import numpy as np
import sympy as sp
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
n = 0
difference = 1


def calcDissipationRate(u,u_prime,U_mean):
    dt = 1.0 / fs

    du_dt = np.gradient(u_prime, dt)  

    du_dx = du_dt / U_mean

    mean_du_dx_sq = np.mean(du_dx**2)

    return 15.0 * nu * mean_du_dx_sq

#print(fc)
while difference > 0.01 or n == 9:
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
    n += 1

print(round(fc),eta)

#b
BIN = 14

Freq, E_f = sg.welch(up_filt, fs, nperseg=2**(BIN+1))

E_k = E_f*U_filt/(2*np.pi)
k_space = 2*np.pi*Freq/U_filt

y_aksen = E_k*eta/nu**2
x_aksen = k_space*eta


plt.figure(figsize=(10, 5))
plt.loglog(x_aksen,y_aksen)
plt.loglog(x_aksen,x_aksen**(-5/3))
plt.ylabel("$E_k$η/ν\u00b2")
plt.xlabel('kη')
plt.ylim(1e-6)
plt.axvline(x=3e-3, color='r', linestyle='--')
plt.axvline(x=1e-1, color='r', linestyle='--')
plt.title('$E_k$η/ν\u00b2 vs. kη')
plt.legend(['The normalized k-space spectra',r'$k^{-5/3}$'])
plt.show()

#d

maxLags = 40000
acf_result = sm.tsa.acf(up_filt, nlags=maxLags, fft=True)
Buu = acf_result
lags = np.arange(len(Buu))

rs = lags * (U_filt/fs)

#Plot this to help visulise the values
plt.figure(figsize=(10, 6))
plt.plot(rs, Buu)
plt.ylabel(r'$B_{uu}$')
plt.xlabel('rs [m]')
plt.axhline(y=0, color='r', linestyle='--')
plt.title('Autocorrelation of u′ (linear scale)')
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 6))
plt.semilogx(rs[1:], Buu[1:])
plt.ylabel(r'$B_{uu}$')
plt.xlabel('rs [m]')
plt.axhline(y=0, color='r', linestyle='--')
plt.title('Autocorrelation of u′ (log scale on x-axis)')
plt.show()

zero_crossing_indices = np.where(Buu <= 0)[0]
first_zero_idx = zero_crossing_indices[0]
L11 = np.trapz(Buu[:first_zero_idx+1], rs[:first_zero_idx+1])

decades = np.log10(L11/eta)
print(f"L11 {L11}")
print(f"Eta {eta}")
print(f"Decades {decades}")