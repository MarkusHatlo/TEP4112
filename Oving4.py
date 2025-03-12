import scipy.io
import numpy as np
import scipy.signal as sg
import matplotlib.pyplot as plt

# Load the .mat file
file_path = "XW_Sample_Sq39_farfield(1).mat"
dataValues = scipy.io.loadmat(file_path)

# Extract required variables
u = np.float64(dataValues['u'].squeeze())
nu = (dataValues['nu']).item()
fs = (dataValues['fs']).item()

#1

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

#2
BIN = 14

Freq, E_f = sg.welch(up_filt, fs, nperseg=2**(BIN+1))

#3
E_k = E_f*U_filt/(2*np.pi)
k_space = 2*np.pi*fc


#4
y_aksen = E_k*eta/nu**2
x_aksen = k_space*Freq


plt.figure(figsize=(10, 5))
plt.loglog(x_aksen,y_aksen)
plt.loglog(x_aksen,1e15*x_aksen**(-5/3))
plt.ylabel("E_k*eta/nu**2")
plt.xlabel('K_eta')
plt.ylim(1e-6)
plt.title('E_k*eta/nu**2 vs k_space*f')
plt.show()


# maxLags = 40000
# Buu, lags,mode = np.correlate(up_filt,maxLags,'full')
# rs = lags * (U_filt/fs)

# plt.loglog(rs, Buu)
# plt.xlabel('Buu')
# plt.ylabel('rs [m]')
# plt.show()