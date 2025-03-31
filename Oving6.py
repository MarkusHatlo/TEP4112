import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

dataValues = scipy.io.loadmat('E:/Skole/TEP4112/data/Instantaneous.mat')
nearValues = scipy.io.loadmat('E:/Skole/TEP4112/data/NearWake.mat')
farValues = scipy.io.loadmat('E:/Skole/TEP4112/data/FarWake.mat')

instant_D = dataValues['D']
instant_U_0 = dataValues['U_0']
instant_u_1 = dataValues['u_1']
instant_x_1 = dataValues['x_1']
instant_x_2 = dataValues['x_2']

near_D = nearValues['D']
near_U_0 = nearValues['U_0']
near_U_1 = nearValues['U_1']
near_up_1up_1 = nearValues['up_1up_1']
near_x_1 = nearValues['x_1']
near_x_2 = nearValues['x_2']

far_D = farValues['D']
far_U_0 = farValues['U_0']
far_U_1 = farValues['U_1']
far_up_1up_1 = farValues['up_1up_1']
far_x_1 = farValues['x_1']
far_x_2 = farValues['x_2']

#Task a

instant_u_1_mean = np.mean(instant_u_1)
instant_u_1_prime = instant_u_1 - instant_u_1_mean

instant_x_1_normalized = instant_x_1/instant_D
instant_x_2_normalized = instant_x_2/instant_D
instant_u_1_normalized = instant_u_1/instant_U_0


X1, X2 = np.meshgrid(instant_x_1_normalized, instant_x_2_normalized)

colors = [(0, 0, 0),      # black
          (0.2, 0, 0.4),  # dark blue/purple
          (0.8, 0, 0.4),  # red/magenta
          (1, 0.5, 0),    # orange
          (1, 1, 0)]      # yellow


cm = LinearSegmentedColormap.from_list('velocity_cmap', colors, N=256)

# Create the figure and axis
fig, ax = plt.subplots(figsize=(12, 5))

# Create the heatmap plot
# If u1_u0 needs to be transposed to match the meshgrid, use u1_u0.T
im = ax.pcolormesh(X1, X2, instant_u_1_normalized, cmap=cm, shading='gouraud')

# Add a colorbar
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('$U_1/U_0$', rotation=0, labelpad=15, y=0.5)

# # Add vertical dashed lines
# ax.axvline(x=2, color='black', linestyle='--')
# ax.axvline(x=4, color='black', linestyle='--')
# ax.axvline(x=16, color='black', linestyle='--')
# ax.axvline(x=18, color='black', linestyle='--')

# # Add arrows between specific vertical lines
# arrow_y = 2
# ax.annotate('', xy=(4, arrow_y), xytext=(2, arrow_y),
#             arrowprops=dict(arrowstyle='<->', color='black'))
# ax.annotate('', xy=(18, arrow_y), xytext=(16, arrow_y),
#             arrowprops=dict(arrowstyle='<->', color='black'))

# Set labels
ax.set_xlabel('$x_1/D$')
ax.set_ylabel('$x_2/D$')

# Set the axis limits
ax.set_xlim(np.min(instant_x_1_normalized), np.max(instant_x_1_normalized))
ax.set_ylim(np.min(instant_x_2_normalized), np.max(instant_x_2_normalized))

# Show the plot
plt.tight_layout()
plt.show()

#Task b

near_U_1_normalized = near_U_1/near_U_0
far_U_1_normalized = far_U_1/far_U_0

near_total_columns = near_U_1_normalized.shape[1]
far_total_columns = far_U_1_normalized.shape[1]

columns_per_profile_near = near_total_columns // 9
columns_per_profile_far = far_total_columns // 9

near_wake_profiles = []
near_wake_positions = []

for i in range(9):
    start_idx = i * columns_per_profile_near
    end_idx = min((i + 1) * columns_per_profile_near, near_total_columns)
    
    avg_profile = np.mean(near_U_1_normalized[:, start_idx:end_idx], axis=1)
    near_wake_profiles.append(avg_profile)
    
    pos = np.mean(near_x_1[0, start_idx:end_idx]) / near_D[0, 0]
    near_wake_positions.append(pos)

far_wake_profiles = []
far_wake_positions = []

for i in range(9):
    start_idx = i * columns_per_profile_far
    end_idx = min((i + 1) * columns_per_profile_far, far_total_columns)
    
    # Average the velocity profiles over this streamwise extent
    avg_profile = np.mean(far_U_1_normalized[:, start_idx:end_idx], axis=1)
    far_wake_profiles.append(avg_profile)
    
    # Calculate approximate x1/D position (center of each segment)
    pos = np.mean(far_x_1[0, start_idx:end_idx]) / far_D[0, 0]
    far_wake_positions.append(pos)

near_x_2_normalized = near_x_2 / near_D
far_x_2_normalized = far_x_2 / far_D

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), sharey=True)

# Plot near wake profiles
for i, profile in enumerate(near_wake_profiles):
    ax1.plot(profile, near_x_2_normalized[:, 0], label=f'x₁/D ≈ {near_wake_positions[i]:.1f}')

# Plot far wake profiles
for i, profile in enumerate(far_wake_profiles):
    ax2.plot(profile, far_x_2_normalized[:, 0], label=f'x₁/D ≈ {far_wake_positions[i]:.1f}')

# Set labels and legends
ax1.set_xlabel('$U_1/U_0$')
ax1.set_ylabel('$x_2/D$')
ax1.set_title('Near Wake Region')
ax1.legend()
ax1.grid(True)

ax2.set_xlabel('$U_1/U_0$')
ax2.set_title('Far Wake Region')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()

#c

near_turbulence_intensity = np.sqrt(near_up_1up_1) / near_U_0 * 100
far_turbulence_intensity = np.sqrt(far_up_1up_1) / far_U_0 * 100


near_total_columns = near_turbulence_intensity.shape[1]
far_total_columns = far_turbulence_intensity.shape[1]

num_profiles = 5

columns_per_profile_near = near_total_columns // num_profiles
columns_per_profile_far = far_total_columns // num_profiles

near_wake_ti_profiles = []
near_wake_positions = []

for i in range(num_profiles):
    start_idx = i * columns_per_profile_near
    end_idx = min((i + 1) * columns_per_profile_near, near_total_columns)
    
    avg_profile = np.nanmean(near_turbulence_intensity[:, start_idx:end_idx], axis=1)

        
    near_wake_ti_profiles.append(avg_profile)
    
    pos = np.mean(near_x_1[0, start_idx:end_idx]) / near_D[0, 0]
    near_wake_positions.append(pos)

far_wake_ti_profiles = []
far_wake_positions = []

for i in range(num_profiles):
    start_idx = i * columns_per_profile_far
    end_idx = min((i + 1) * columns_per_profile_far, far_total_columns)

    avg_profile = np.nanmean(far_turbulence_intensity[:, start_idx:end_idx], axis=1)
        
    far_wake_ti_profiles.append(avg_profile)
    
    pos = np.mean(far_x_1[0, start_idx:end_idx]) / far_D[0, 0]
    far_wake_positions.append(pos)


near_x_2_normalized = near_x_2 / near_D
far_x_2_normalized = far_x_2 / far_D

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), sharey=True)

for i, profile in enumerate(near_wake_ti_profiles):
    ax1.plot(profile, near_x_2_normalized[:, 0], label=f'x₁/D ≈ {near_wake_positions[i]:.1f}')

for i, profile in enumerate(far_wake_ti_profiles):
    ax2.plot(profile, far_x_2_normalized[:, 0], label=f'x₁/D ≈ {far_wake_positions[i]:.1f}')

ax1.set_xlabel('Turbulence Intensity (%)')
ax1.set_ylabel('$x_2/D$')
ax1.set_title('Near Wake Region - Turbulence Intensity')
ax1.legend()
ax1.grid(True)

ax2.set_xlabel('Turbulence Intensity (%)')
ax2.set_title('Far Wake Region - Turbulence Intensity')
ax2.legend()
ax2.grid(True)


max_ti = 30 
ax1.set_xlim(0, max_ti)
ax2.set_xlim(0, max_ti)

plt.tight_layout()
plt.show()