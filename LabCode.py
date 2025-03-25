import scipy.io
import numpy as np
import sympy as sp
import scipy.signal as sg
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats as sps
from scipy.stats import norm
from scipy.optimize import curve_fit

plt.rcParams['font.size'] = 22

def calcDissipationRate(u,u_prime,U_mean):
    dt = 1.0 / fs

    du_dt = np.gradient(u_prime, dt)  

    du_dx = du_dt / U_mean

    mean_du_dx_sq = np.mean(du_dx**2)

    return 15.0 * nu * mean_du_dx_sq

y_aksen_list = []

def plot_spectra():
    BIN = 14

    Freq, E_f = sg.welch(up_filt, fs, nperseg=2**(BIN+1))


    E_k = E_f*U_filt/(2*np.pi)
    k_space = 2*np.pi*Freq/U_filt

    y_aksen = E_k*eta/nu**2
    x_aksen = k_space*eta

    y_aksen_list.append(E_k*eta/nu**2)

    if len(y_aksen_list) == 3:

        plt.figure(figsize=(10, 5))
        plt.loglog(x_aksen,y_aksen_list[0])
        plt.loglog(x_aksen,y_aksen_list[1])
        plt.loglog(x_aksen,y_aksen_list[2])
        plt.loglog(x_aksen,2*x_aksen**(-5/3))
        plt.ylabel("$E_k$η/ν\u00b2")
        plt.xlabel('kη')
        plt.ylim(1e-6)
        #plt.axvline(x=3e6, color='r', linestyle='--')
        #plt.axvline(x=1e8, color='r', linestyle='--')
        plt.title(f'$E_k$η/ν\u00b2 vs. kη')
        plt.legend(['25M','52M','80M',r'$k^{-5/3}$'])
        plt.show()

def plot_autocorrelation():
    #Plot this to help visulise the values
    # plt.figure(figsize=(10, 6))
    # plt.plot(rs, Buu)
    # plt.ylabel(r'$B_{uu}$')
    # plt.xlabel('rs [m]')
    # plt.axhline(y=0, color='r', linestyle='--')
    # plt.title(f'Autocorrelation of u′ (linear scale) - File: {i}M')
    # plt.grid(True)
    # plt.show()

    plt.figure(figsize=(10, 6))
    plt.semilogx(rs[1:], Buu[1:])
    plt.ylabel(r'$B_{uu}$')
    plt.xlabel('rs [m]')
    plt.axhline(y=0, color='r', linestyle='--')
    plt.title(f'Autocorrelation of u′ (log scale on x-axis)- File: {i}M')
    plt.show()

def calc_decades():
    zero_crossing_indices = np.where(Buu <= 0)[0]
    first_zero_idx = zero_crossing_indices[0]
    L11 = np.trapz(Buu[:first_zero_idx+1], rs[:first_zero_idx+1])

    decades = np.log10(L11/eta)
    print(f"L11 {L11}")
    print(f"Eta {eta}")
    print(f"Decades {decades}")

decay = []
decay_position = []
powerlaw = []
u_filt_kde_list = []

def plot_statistics():
    x_range = np.linspace(-4, 4, 1000)
    gaussian_pdf = sps.norm.pdf(x_range)

    u_filt_normalized = up_filt/u_filt_std
    u_filt_kde = sps.gaussian_kde(u_filt_normalized)

    print(f"This is the standard deviation for flow {i:.3g}M : {u_filt_std}\n")

    u_skewness = sps.skew(u_filt)
    print(f"This is the skewness for flow {i:.3g}M : {u_skewness}\n")

    u_kurtosis= sps.kurtosis(u_filt, fisher=False)    
    print(f"This is the kurtosis for flow {i:.3g}M: {u_kurtosis}\n")

    u_filt_kde_list.append(u_filt_kde(x_range))

    if len(u_filt_kde_list) == 3:

        plt.figure(figsize=(10, 6))
        plt.plot(x_range, u_filt_kde_list[0], 'r-', linewidth=2, label='Turbulence Grid', color = 'red')
        plt.plot(x_range, u_filt_kde_list[1], 'r-', linewidth=2, label='Turbulence Grid', color = 'blue')
        plt.plot(x_range, u_filt_kde_list[2], 'r-', linewidth=2, label='Turbulence Grid', color = 'green')
        plt.plot(x_range, gaussian_pdf, 'k--', linewidth=2, label='Gaussian (S=0, K=3)')
        plt.xlabel(r'$u\' / \sigma$')
        plt.ylabel('Probability Density')
        plt.title('PDF of Normalized Velocity Fluctuations')
        plt.legend()
        plt.grid(True)
        plt.xlim(-4, 4)
        plt.show()


# Load the .mat file
#for i in [25,29,34,38,43,47,52,57,61,66,70,75,80]:
for i in [25,52,80]:
#for i in [25]:
    file_path = f"{i}M_processed.mat"
    dataValues = scipy.io.loadmat(fr"Group 8\processed\{file_path}")

    print("-------------------")
    print(f"This is file {i}M")
    print("-------------------\n")

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

    maxLags = 40000
    acf_result = sm.tsa.acf(up_filt, nlags=maxLags, fft=True)
    Buu = acf_result
    lags = np.arange(len(Buu))
    rs = lags * (U_filt/fs)

    u_filt_std = np.std(up_filt)

    # plot_spectra()

    # plot_autocorrelation()

    # calc_decades()

    # plot_statistics()



    # k_full = 1.5 * (u_filt_std**2)
    # print("TKE (full isotropy) =", k_full)

    decay.append(round(float(np.mean(up_filt**2)),4))
    decay_position.append(i)


    print('Calculations completed\n')

#------------------------------------------------------------------------------
#Kode for decay
first_position = decay_position[0]
first_value = decay[0]

scaling_coefficient = first_value / (first_position**(-1.14))

powerlaw = [scaling_coefficient * (i**(-1.14)) for i in decay_position]


def power_law(x, A, n):
    return A * x**(-n)

params, params_covariance = curve_fit(power_law, decay_position, decay)
A_fit, n_fit = params
powerlaw_fit = power_law(np.array(decay_position), A_fit, n_fit)

print(f"Fitted parameters: A = {A_fit:.4f}, n = {n_fit:.4f}")
print(f'This is the decay numbers A = {scaling_coefficient:.4g} deacy = {decay} and position = {decay_position}')

std_curve_fit = np.sqrt(np.diag(params_covariance))
print(f"Parameter errors:  Δn = {std_curve_fit[1]:.4g}")


plt.figure(figsize=(10, 6))
plt.plot(decay_position,decay,'o', linestyle='')
plt.plot(decay_position,powerlaw_fit, color = 'green')
plt.ylabel(r'$\bar{u_1^2}$')
plt.xlabel('Position in the streamwise direction')
plt.title(r'The decay of $\bar{u_1^2}$ in the streamwise direction')
plt.show()
#------------------------------------------------------------------------------










