import scipy.io
import numpy as np
import scipy.signal as sg

# Load the .mat file
file_path = "XW_Sample_Sq39_farfield(1).mat"
dataValues = scipy.io.loadmat(file_path)
print(dataValues.keys())

# Extract required variables
u = np.float(dataValues['u'].squeeze)
nu = np.float(dataValues['nu'])
fs = np.float(dataValues['fs'])

fc = fs/2.5
LPF_order = 7
k = 0
difference = 1


def calcDissipationRate(u,u_prime,U_mean):
    dt = 1.0 / fc

    du_dt = np.gradient(u_prime, dt)  

    du_dx = du_dt / U_mean

    mean_du_dx_sq = np.mean(du_dx**2)

    return 15.0 * nu * mean_du_dx_sq


while difference < 0.01 or k == 9:
    b, a = sg.butter(LPF_order, fc*2/fs, btype='lowpass', output='ba')
    u_filt = sg.filtfilt(b,a,u)
    U_filt = np.mean(u)
    up_filt = u - U_filt
    dissipation = calcDissipationRate(u_filt,up_filt,U_filt)
    f_eta = (nu**3/dissipation)**(1/4)
    kolmogorov_frequency = U_filt/(2*np.pi*f_eta)
    fc_new = 1.1*f_eta
    difference = abs(fc-fc_new)
    fc = fc_new
    k += 1
    print(k)

print(fc_new)