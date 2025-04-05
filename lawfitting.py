import scipy.io
import numpy as np
import sympy as sp
import scipy.signal as sg
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats as sps
from scipy.stats import norm
from scipy.optimize import curve_fit

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
eta_list = []
dissipation_list = []

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

    if len(u_filt_kde_list) == 1:

        plt.figure(figsize=(10, 6))
        plt.semilogy(x_range, u_filt_kde_list[0], linewidth=2, label=r'$25M$', color = 'red')
        # plt.semilogy(x_range, u_filt_kde_list[1], linewidth=2, label=r'$52M$', color = 'blue')
        # plt.semilogy(x_range, u_filt_kde_list[2], linewidth=2, label=r'$80M$', color = 'green')
        plt.semilogy(x_range, gaussian_pdf, 'k--', linewidth=2, label=r'$\text{Gaussian}~(S=0,~K=3)$')
        plt.xlabel(r'$u\' / \sigma$')
        plt.ylabel('Probability Density')
        plt.title('PDF of Normalized Velocity Fluctuations',)
        plt.legend()
        plt.grid(True)
        plt.xlim(-4, 4)
        plt.show()

# Function for power law with virtual origin
def power_law(x, A, n, x0):
    return A * ((x - x0) ** (-n))

# Function to find the optimal virtual origin based on Lavoie et al. (2007) method
def find_optimal_virtual_origin(positions, values, x0_range):
    """
    Implements the method from Lavoie et al. (2007) to find the optimal virtual origin
    that gives the most consistent value of n over a range of fitting windows.
    """
    mq_values_by_x0 = {}
    plateau_widths = {}
    
    # For each potential virtual origin
    for x0 in x0_range:
        mq_values = []
        start_positions = []
        
        # Try different starting positions for the fit
        for start_idx in range(len(positions) - 3):  # Ensure at least 4 points for each fit
            # Get positions from this starting point to the end
            pos_window = positions[start_idx:]
            val_window = values[start_idx:]
            
            # Skip if any adjusted position is non-positive
            if any((p - x0) <= 0 for p in pos_window):
                continue
            
            # Perform a power law fit with fixed x0
            try:
                def fixed_x0_power_law(x, A, n):
                    return A * ((x - x0) ** (-n))
                
                params, _ = curve_fit(fixed_x0_power_law, pos_window, val_window)
                A_fit, n_fit = params
                
                mq_values.append(n_fit)
                start_positions.append(positions[start_idx])
            except:
                continue
        
        # If we couldn't get enough fits, skip this x0
        if len(mq_values) < 3:
            continue
        
        # Find the longest plateau (range where n is approximately constant)
        mq_values_by_x0[x0] = (start_positions, mq_values)
        
        # Measure the width of the plateau
        plateau_width = 0
        for i in range(len(mq_values) - 1):
            if abs(mq_values[i+1] - mq_values[i]) < 0.02:  # Tolerance for what's considered constant
                current_width = 2
                for j in range(i+2, len(mq_values)):
                    if abs(mq_values[j] - mq_values[i]) < 0.05:  # Slightly larger tolerance for the whole plateau
                        current_width += 1
                    else:
                        break
                
                plateau_width = max(plateau_width, current_width)
        
        plateau_widths[x0] = plateau_width
    
    # Find the x0 with the widest plateau
    if not plateau_widths:
        return None, None, {}
    
    best_x0 = max(plateau_widths, key=plateau_widths.get)
    
    # Find the average n value within the plateau for the best x0
    start_pos, mq_vals = mq_values_by_x0[best_x0]
    
    # Get the start and end index of the longest plateau
    best_plateau_start = 0
    best_plateau_width = 0
    for i in range(len(mq_vals) - 1):
        if abs(mq_vals[i+1] - mq_vals[i]) < 0.02:
            current_width = 2
            for j in range(i+2, len(mq_vals)):
                if abs(mq_vals[j] - mq_vals[i]) < 0.05:
                    current_width += 1
                else:
                    break
            
            if current_width > best_plateau_width:
                best_plateau_width = current_width
                best_plateau_start = i
    
    # Calculate the average n within the plateau
    best_n = np.mean(mq_vals[best_plateau_start:best_plateau_start+best_plateau_width])
    
    return best_x0, best_n, mq_values_by_x0

# Function to analyze decay according to power law
def analyze_decay():
    if len(decay) < 5:
        print("Not enough data points for decay analysis")
        return
    
    print("\n--- Power Law Decay Analysis ---")
    print(f"Decay positions: {decay_position}")
    print(f"Decay values: {decay}")
    
    # Only use positions after x/M = 30 where turbulence is likely fully developed
    min_position = 30
    valid_indices = [i for i, pos in enumerate(decay_position) if pos >= min_position]
    
    if len(valid_indices) < 5:
        print("Not enough data points in the fully developed region")
        return
    
    valid_positions = [decay_position[i] for i in valid_indices]
    valid_decay = [decay[i] for i in valid_indices]
    
    print(f"Using positions from {min_position}M onwards: {valid_positions}")
    
    # Try a more focused range of virtual origins based on physics
    # Typically x₀ is between 0 and 15 mesh lengths upstream of the grid
    x0_range = np.linspace(-5, 15, 41)
    
    # Try logarithmic fitting to reduce sensitivity to noise
    # Transform to log space for more stable fitting
    log_positions = np.log(valid_positions)
    log_decay = np.log(valid_decay)
    
    # Simple linear fit in log space to get initial estimate
    slope, intercept = np.polyfit(log_positions, log_decay, 1)
    initial_n = -slope
    print(f"Initial estimate from log-log fit: n ≈ {initial_n:.4f}")
    
    # Use the improved method for finding x₀
    best_x0, best_n, mq_data = find_optimal_virtual_origin(valid_positions, valid_decay, x0_range)
    
    
    if best_x0 is None:
        print("Could not find a stable virtual origin with Lavoie's method")
    else:
        print(f"\nLavoie method results:")
        print(f"Optimal virtual origin (x0): {best_x0:.4f}")
        print(f"Corresponding decay exponent (n): {best_n:.4f}")
        
        # Fit with the optimal virtual origin
        def fixed_x0_power_law(x, A, n):
            return A * ((x - best_x0) ** (-n))
        
        params, params_covariance = curve_fit(fixed_x0_power_law, decay_position, decay)
        A_fit, n_fit = params
        
        std_errors = np.sqrt(np.diag(params_covariance))
        print(f"Final fit with fixed x0: A = {A_fit:.4g}, n = {n_fit:.4f} ± {std_errors[1]:.4f}")
        
        # Calculate the power-law prediction
        x_range = np.linspace(min(decay_position), max(decay_position), 100)
        powerlaw_fit = fixed_x0_power_law(x_range, A_fit, n_fit)
        
        # Plot the results
        plt.figure(figsize=(12, 8))
        
        # Plot 1: Original data with fitted curve
        plt.subplot(2, 1, 1)
        plt.loglog(np.array(decay_position) - best_x0, decay, 'o', markersize=8)
        
        # Create a smoother curve for the fit line
        x_smooth = np.logspace(
            np.log10(min(np.array(decay_position) - best_x0)), 
            np.log10(max(np.array(decay_position) - best_x0)), 
            100
        )
        y_smooth = A_fit * (x_smooth ** (-n_fit))
        
        plt.loglog(x_smooth, y_smooth, 'r-', linewidth=2)
        plt.xlabel(r'$(x - x_0)/M$')
        plt.ylabel(r'$\overline{u^{\prime 2}}$')
        plt.title(f'Power Law Fit: $\\overline{{u^{{\\prime 2}}}} = {A_fit:.4g} \\cdot ((x - {best_x0:.4g})/M)^{{-{n_fit:.4f}}}$')
        plt.grid(True, which="both", ls="-", alpha=0.2)
        
        # Plot 2: Linear scale with fit
        plt.subplot(2, 1, 2)
        plt.plot(decay_position, decay, 'o', markersize=8)
        plt.plot(x_range, fixed_x0_power_law(x_range, A_fit, n_fit), 'r-', linewidth=2)
        plt.xlabel('Position (x/M)')
        plt.ylabel(r'$\overline{u^{\prime 2}}$')
        plt.title('Power Law Fit (Linear Scale)')
        plt.grid(True)
        
        plt.tight_layout()
        plt.show()
        
        # Plot the variation of n with starting position for the best x0
        if best_x0 in mq_data:
            start_positions, mq_values = mq_data[best_x0]
            
            plt.figure(figsize=(10, 6))
            plt.plot(start_positions, mq_values, 'o-', markersize=8)
            plt.axhline(y=best_n, color='r', linestyle='--')
            plt.xlabel(r'Starting Position $(x/M)$')
            plt.ylabel(r'Decay Exponent $(n)$')
            #plt.title(f'Stability of Decay Exponent for $x_0 = {best_x0:.4f}$')
            plt.grid(True)
            plt.show()
            
            # Also compare to theoretical values
            plt.figure(figsize=(10, 6))
            plt.axhline(y=1.0, color='k', linestyle='--', label='Complete Self-Preservation (n=1.0)')
            plt.axhline(y=1.2, color='g', linestyle='--', label='Saffman Turbulence (n=1.2)')
            plt.axhline(y=10/7, color='b', linestyle='--', label='Batchelor Turbulence (n=10/7)')
            plt.plot(start_positions, mq_values, 'ro-', markersize=8, label=f'Measured (n={best_n:.4f})')
            plt.xlabel('Starting Position (x/M)')
            plt.ylabel('Decay Exponent (n)')
            plt.title('Comparison with Theoretical Decay Exponents')
            plt.legend()
            plt.grid(True)
            plt.show()
            
    # Method 2: Try also the traditional curve_fit with all parameters free
    try:
        params, params_covariance = curve_fit(power_law, decay_position, decay)
        A_fit, n_fit, x0_fit = params
        powerlaw_fit = power_law(np.array(decay_position), A_fit, n_fit, x0_fit)
        
        std_errors = np.sqrt(np.diag(params_covariance))
        print(f"\nTraditional curve_fit results:")
        print(f"A = {A_fit:.4g}, n = {n_fit:.4f} ± {std_errors[1]:.4f}, x0 = {x0_fit:.4g} ± {std_errors[2]:.4f}")
        
        # Plot the results
        plt.figure(figsize=(10, 6))
        plt.loglog(np.array(decay_position) - x0_fit, decay, 'o', markersize=8)
        
        # Create a smoother curve for the fit line
        x_smooth = np.logspace(
            np.log10(min(np.array(decay_position) - x0_fit)), 
            np.log10(max(np.array(decay_position) - x0_fit)), 
            100
        )
        y_smooth = A_fit * (x_smooth ** (-n_fit))
        
        plt.loglog(x_smooth, y_smooth, 'r-', linewidth=2)
        plt.xlabel(r'$(x - x_0)/M$')
        plt.ylabel(r'$\overline{u^{\prime 2}}$')
        plt.title(f'Traditional Curve Fit: $\\overline{{u^{{\\prime 2}}}} = {A_fit:.4g} \\cdot ((x - {x0_fit:.4g})/M)^{{-{n_fit:.4f}}}$')
        plt.grid(True, which="both", ls="-", alpha=0.2)
        plt.show()
        
    except Exception as e:
        print(f"Traditional curve_fit failed: {e}")
     
    # Method 3: Check if the Taylor microscale grows linearly with (x-x0)
    if 'dissipation_list' in globals() and len(dissipation_list) == len(decay):
        # Calculate Taylor microscale according to λ² = 15ν⟨u²⟩/ε
        lambda_squared = []
        for i in range(len(decay)):
            lambda_sq = 15 * nu * decay[i] / dissipation_list[i]
            lambda_squared.append(lambda_sq)
        
        # If we found a good x0 with Lavoie's method, check if λ² vs (x-x0) is linear
        if best_x0 is not None:
            adjusted_positions = np.array(decay_position) - best_x0
            
            # Only use positions where (x-x0) > 0
            valid_indices = adjusted_positions > 0
            if np.sum(valid_indices) > 2:
                valid_positions = adjusted_positions[valid_indices]
                valid_lambda_sq = np.array(lambda_squared)[valid_indices]
                
                # Linear fit
                slope, intercept = np.polyfit(valid_positions, valid_lambda_sq, 1)
                
                # Calculate mλ from slope (should equal -n if power law is valid)
                if slope > 0:  # Ensure positive slope as λ² should grow with x
                    m_lambda = -10 * nu / slope
                    
                    print(f"\nTaylor microscale analysis:")
                    print(f"λ² grows linearly with (x-x0) with slope = {slope:.4g}")
                    print(f"This gives mλ = {m_lambda:.4f}")
                    print(f"Comparison: mq = {best_n:.4f}, mλ = {m_lambda:.4f}")
                    
                    # Plot λ² vs (x-x0)
                    plt.figure(figsize=(10, 6))
                    plt.scatter(valid_positions, valid_lambda_sq, s=80)
                    
                    x_fit = np.linspace(min(valid_positions), max(valid_positions), 100)
                    y_fit = slope * x_fit + intercept
                    
                    plt.plot(x_fit, y_fit, 'r-', linewidth=2)
                    plt.xlabel(r'$(x - x_0)/M$')
                    plt.ylabel(r'$\lambda^2$ [m$^2$]')
                    plt.title(f'Taylor Microscale Growth: $\\lambda^2 = {slope:.4g} \\cdot (x-x_0)/M + {intercept:.4g}$')
                    plt.grid(True)
                    plt.show()

# Load the .mat file
u_filt_kde_list = []

# Uncomment this line for all positions
# positions = [25, 29, 34, 38, 43, 47, 52, 57, 61, 66, 70, 75, 80]

# For testing with fewer positions
positions = [25, 29, 34, 38, 43, 47, 52, 57, 61, 66, 70, 75, 80]

for i in positions:
    file_path = f"{i}M_processed.mat"
    dataValues = scipy.io.loadmat(fr"Group 8\processed\{file_path}")

    print("-------------------")
    print(f"This is file {i}M")
    print("-------------------\n")

    # Extract required variables
    u = np.float64(dataValues['u'].squeeze())
    nu = (dataValues['nu']).item()
    fs = (dataValues['fs']).item()

    fc = fs/2.5
    fc_new = 0
    LPF_order = 7
    n = 0
    difference = 1

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
    dissipation_list.append(float(dissipation))

    # maxLags = 40000
    # acf_result = sm.tsa.acf(up_filt, nlags=maxLags, fft=True)
    # Buu = acf_result
    # lags = np.arange(len(Buu))
    # rs = lags * (U_filt/fs)

    # u_filt_std = np.std(up_filt)
    
    # Store the velocity variance (turbulent kinetic energy assuming isotropy)
    decay.append(float(np.mean(up_filt**2)))
    decay_position.append(i)

    print('Calculations completed\n')

# Run the power law analysis once all data is collected
analyze_decay()