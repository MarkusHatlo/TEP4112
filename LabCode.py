import scipy.io
import numpy as np
import sympy as sp
import scipy.signal as sg
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats as sps
from scipy.stats import norm
from scipy.optimize import curve_fit
from matplotlib.ticker import LogLocator
import matplotlib.ticker as plti

plt.rcParams.update({
    'text.usetex': False,           # Use LaTeX to render text
    'font.family': 'serif',        # Use serif font
    'font.size': 22                # Keep your desired font size
})

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

L11_list = []
decades_list = []

def calc_decades():
    print(f"Min Buu: {np.min(Buu)}, Max Buu: {np.max(Buu)}")
    if i == 57:
        zero_crossing_indices = np.where(Buu <= 0.006)[0]
        first_zero_idx = zero_crossing_indices[0]
    elif i == 61:
        zero_crossing_indices = np.where(Buu <= 0.0129)[0]
        first_zero_idx = zero_crossing_indices[0]
    elif i == 70:
        zero_crossing_indices = np.where(Buu <= 0.009)[0]
        first_zero_idx = zero_crossing_indices[0]
    elif i == 75:
        zero_crossing_indices = np.where(Buu <= 0.005)[0]
        first_zero_idx = zero_crossing_indices[0]
    else:
        zero_crossing_indices = np.where(Buu <= 0)[0]
        first_zero_idx = zero_crossing_indices[0]
    L11 = np.trapz(Buu[:first_zero_idx+1], rs[:first_zero_idx+1])

    decades = np.log10(L11/eta)
    print(f"L11 {L11}")
    print(f"Eta {eta}")
    print(f"Decades {decades}")

    L11_list.append(float(L11))
    decades_list.append(float(decades))

    if len(L11_list) == 13:
            print(f"L11 {L11_list}")
            print(f"Eta {eta_list}")
            print(f"Decades {decades_list}")
    #     plt.figure(figsize=(10, 5))
    #     plt.plot([25,29,34,38,43,47,52,57,61,66,70,75,80],L11_list)
    #     plt.xlabel('Streamwise position from the turbulence grid [M]')
    #     plt.ylabel('Integral length scale [m]')
    #     plt.title('Integral length scale at each streamwise position')
    #     plt.show()

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
        plt.semilogy(x_range, u_filt_kde_list[0], linewidth=2, label=r'$25M$', color = 'red')
        plt.semilogy(x_range, u_filt_kde_list[1], linewidth=2, label=r'$52M$', color = 'blue')
        plt.semilogy(x_range, u_filt_kde_list[2], linewidth=2, label=r'$80M$', color = 'green')
        plt.semilogy(x_range, gaussian_pdf, 'k--', linewidth=2, label=r'$\text{Gaussian}~(S=0,~K=3)$')
        plt.xlabel(r'$u\' / \sigma$')
        plt.ylabel('Probability Density')
        plt.title('PDF of Normalized Velocity Fluctuations',)
        plt.legend()
        plt.grid(True)
        plt.xlim(-4, 4)
        plt.show()

eta_list = []

# Load the .mat file
for i in [25,29,34,38,43,47,52,57,61,66,70,75,80]:
# for i in [25,29,52,80]:
# for i in [70]:
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

    eta_list.append(float(eta))
    if i == 25:
        U0 = U_filt

    # maxLags = 40000
    # acf_result = sm.tsa.acf(up_filt, nlags=maxLags, fft=True)
    # Buu = acf_result
    # lags = np.arange(len(Buu))
    # rs = lags * (U_filt/fs)

    # u_filt_std = np.std(up_filt)

    #plot_spectra()

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

decay_ref = decay
decay = [val / (U0**2) for val in decay_ref]

# powerlaw = [scaling_coefficient * (i**(-1.14)) for i in decay_position]


def power_law(x, A, n, x0):
    return A * (x-x0)**(-n)


params, params_covariance = curve_fit(power_law, decay_position, decay)
A_fit, n_fit, x0_fit = params
powerlaw_fit = power_law(np.array(decay_position), A_fit, n_fit, x0_fit)

print(f"Fitted parameters: A = {A_fit:.4f}, n = {n_fit:.4f} x0 = {x0_fit:.4g}")
print(f'This is the decay numbers A = {scaling_coefficient:.4g} deacy = {decay} and position = {decay_position}')

std_curve_fit = np.sqrt(np.diag(params_covariance))
print(f"Parameter errors:  Δn = {std_curve_fit[1]:.4g}")

decay_position_riktig = decay_position[2:]
decay_riktig = decay[2:]
# powerlaw_riktig = powerlaw[1:]

params_riktig, params_covariance_riktig = curve_fit(power_law, decay_position_riktig, decay_riktig)
A_fit_riktig, n_fit_riktig, x0_fit_riktig = params_riktig
powerlaw_fit_riktig = power_law(np.array(decay_position_riktig), A_fit_riktig, n_fit_riktig, x0_fit_riktig)

print(f"Fitted parameters: A = {A_fit_riktig:.4f}, n = {n_fit_riktig:.4f} x0 = {x0_fit_riktig:.4g}")
print(f'This is the decay numbers A = {scaling_coefficient:.4g} deacy = {decay_riktig} and position = {decay_position_riktig}')

std_curve_fit_riktig = np.sqrt(np.diag(params_covariance_riktig))
print(f"Parameter errors:  Δn = {std_curve_fit_riktig[1]:.4g}")

# n1 = power_law(np.array(decay_position), A_fit, 1, x0_fit)
# n12 = power_law(np.array(decay_position), A_fit, 1.2, x0_fit)
# n14 = power_law(np.array(decay_position), A_fit, 1.4, x0_fit)

first_x = decay_position[0]
first_y = decay[0]

# Create a function for theoretical lines that pass through the first point
def theoretical_line(x, n):
    # Calculate coefficient A so the line passes through (first_x, first_y)
    A = first_y * (first_x**n)
    return A * (x**(-n))

# Calculate the three theoretical lines
n1_line = theoretical_line(np.array(decay_position), 1.0)
n12_line = theoretical_line(np.array(decay_position), 1.2)
n14_line = theoretical_line(np.array(decay_position), 1.4)


plt.figure(figsize=(10, 6))
plt.loglog(decay_position,decay,'o', linestyle='')
plt.loglog(decay_position,powerlaw_fit, color = 'red', label=r'Fitting all points $(n=1.45, x_0 = -5.2)$')
plt.loglog(decay_position_riktig,powerlaw_fit_riktig, color = 'green', label=r'Fitting without 1. point $(n=1.19, x_0 = 3.4)$')
plt.loglog(decay_position,n1_line, color='k', linestyle='--', label=r'Complete self-preservation $(n=1.0)$')
plt.loglog(decay_position,n12_line, color='g', linestyle='--', label=r'Saffman turbulence $(n=1.2)$')
plt.loglog(decay_position,n14_line, color='b', linestyle='--', label=r'Batchelor turbulence $(n=10/7)$')
plt.xticks([30, 40, 60])
plt.xlim(20, None)
# plt.yticks([0.5, 0.7, 1.0, 1.6])
plt.ylabel(r'$\overline{u^{\prime 2}}/U_0$ [-]')
plt.xlabel(r'$x_1/M$ - $x_0/M$ [-]')
#plt.title(r'The decay of $\overline{u^{\prime 2}}$ in the streamwise direction')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.show()
#------------------------------------------------------------------------------


# # Test plot to verify fonts
# x = np.linspace(0, 10, 100)
# y = np.sin(x)
# plt.figure(figsize=(8, 5))
# plt.plot(x, y)
# plt.title(r"Test Plot with Latin Modern: $\sin(x)$")
# plt.xlabel(r"$x$ values")
# plt.ylabel(r"$\sin(x)$")
# plt.tight_layout()
# plt.show()


n_list = []
positions = []

for i in range(8):
    decay_position_loop = decay_position[i:]
    decay_loop = decay[i:]

    if len(decay_position_loop) <= 3:
        print(f"Stopping at iteration {i}: not enough data points left")
        break

    try:
        # Increase maxfev and provide better initial guess based on previous fit
        if i == 0:
            params_loop, params_covariance_loop = curve_fit(power_law, decay_position_loop, decay_loop, 
                                                          maxfev=100)
        else:
            # Use previous parameters as initial guess
            params_loop, params_covariance_loop = curve_fit(power_law, decay_position_loop, decay_loop, 
                                                          p0=params_prev, 
                                                          maxfev=1000)
        
        A_fit_loop, n_fit_loop, x0_fit_loop = params_loop
        powerlaw_fit_loop = power_law(np.array(decay_position), A_fit_loop, n_fit_loop, x0_fit_loop)
        n_list.append(n_fit_loop)
        positions.append(decay_position[i])  # Store corresponding position
        
        # Save for next iteration
        params_prev = params_loop
        
        print(f"Iteration {i}: n = {n_fit_loop:.4f}")
        
    except RuntimeError as e:
        print(f"Skipping iteration {i}: {str(e)}")
        continue

# Only plot if we have data
if n_list:
    plt.figure(figsize=(10,10))
    plt.plot(positions, n_list, color='red', marker='o', label='n exponent')
    plt.axhline(y=1.0, color='k', linestyle='--', label=r'Complete Self-Preservation $(n=1.0)$')
    plt.axhline(y=1.2, color='g', linestyle='--', label=r'Saffman Turbulence $(n=1.2)$')
    plt.axhline(y=10/7, color='b', linestyle='--', label=r'Batchelor Turbulence $(n=10/7)$')
    
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plti.ScalarFormatter())
    ax.xaxis.set_minor_formatter(plti.ScalarFormatter())
    
    plt.ylabel('Value of n')
    plt.xlabel(r'$x_1/M$')
    #plt.title('Variation of n exponent with reduced data points')
    plt.legend()
    plt.show()
else:
    print("No successful fits to plot")







# %%
