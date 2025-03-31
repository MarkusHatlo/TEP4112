import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Your existing code for loading data
dataValues = scipy.io.loadmat('E:/Skole/TEP4112/data/Instantaneous.mat')
nearValues = scipy.io.loadmat('E:/Skole/TEP4112/data/NearWake.mat')
farValues = scipy.io.loadmat('E:/Skole/TEP4112/data/FarWake.mat')

# Print more details about the data structure
print("Near wake data:")
for key in nearValues.keys():
    if key not in ['__header__', '__version__', '__globals__']:
        print(f"{key} shape: {nearValues[key].shape}, type: {type(nearValues[key])}")
        if np.isnan(nearValues[key]).any():
            print(f"  WARNING: {key} contains NaN values")
        if np.size(nearValues[key]) == 0:
            print(f"  WARNING: {key} is empty")

print("\nFar wake data:")
for key in farValues.keys():
    if key not in ['__header__', '__version__', '__globals__']:
        print(f"{key} shape: {farValues[key].shape}, type: {type(farValues[key])}")
        if np.isnan(farValues[key]).any():
            print(f"  WARNING: {key} contains NaN values")
        if np.size(farValues[key]) == 0:
            print(f"  WARNING: {key} is empty")

# Now let's proceed with caution for the rest of the code
near_D = nearValues['D']
near_U_0 = nearValues['U_0']
near_U_1 = nearValues['U_1']
near_x_1 = nearValues['x_1']
near_x_2 = nearValues['x_2']

far_D = farValues['D']
far_U_0 = farValues['U_0']
far_U_1 = farValues['U_1']
far_x_1 = farValues['x_1']
far_x_2 = farValues['x_2']

# Normalize data - check for zeros to avoid division issues
near_U_1_normalized = np.divide(near_U_1, near_U_0, out=np.zeros_like(near_U_1), where=near_U_0!=0)
far_U_1_normalized = np.divide(far_U_1, far_U_0, out=np.zeros_like(far_U_1), where=far_U_0!=0)

# Calculate the number of columns for each profile
near_total_columns = near_U_1_normalized.shape[1] if len(near_U_1_normalized.shape) > 1 else 1
far_total_columns = far_U_1_normalized.shape[1] if len(far_U_1_normalized.shape) > 1 else 1

# Make sure we have enough data for 9 profiles
columns_per_profile_near = max(1, near_total_columns // 9)
columns_per_profile_far = max(1, far_total_columns // 9)

print(f"\nNear wake U_1 shape: {near_U_1_normalized.shape}")
print(f"Near wake columns per profile: {columns_per_profile_near}")
print(f"Far wake U_1 shape: {far_U_1_normalized.shape}")
print(f"Far wake columns per profile: {columns_per_profile_far}")

# Try a simpler approach first to debug the issue
# Let's see if we can get any profiles at all
# We'll adjust the number of profiles later if needed

# For near wake data
near_wake_profiles = []
near_wake_positions = []

# Only proceed if there are at least 9 profiles possible
if near_total_columns >= 9:
    for i in range(9):
        start_idx = i * columns_per_profile_near
        end_idx = min((i + 1) * columns_per_profile_near, near_total_columns)
        
        # Skip empty slices
        if start_idx >= end_idx:
            continue
            
        # Average the velocity profiles over this streamwise extent
        # Use nanmean to handle NaN values
        avg_profile = np.nanmean(near_U_1_normalized[:, start_idx:end_idx], axis=1)
        
        # Skip if the profile contains only NaN values
        if np.all(np.isnan(avg_profile)):
            continue
            
        near_wake_profiles.append(avg_profile)
        
        # Calculate approximate x1/D position (center of each segment)
        # Handle potential divide by zero
        if near_D.size > 0 and near_D[0, 0] != 0:
            pos = np.mean(near_x_1[0, start_idx:end_idx]) / near_D[0, 0]
        else:
            pos = i  # Fallback to index if D is zero or empty
        near_wake_positions.append(pos)

# Same for far wake data
far_wake_profiles = []
far_wake_positions = []

if far_total_columns >= 9:
    for i in range(9):
        start_idx = i * columns_per_profile_far
        end_idx = min((i + 1) * columns_per_profile_far, far_total_columns)
        
        # Skip empty slices
        if start_idx >= end_idx:
            continue
            
        # Average the velocity profiles over this streamwise extent
        # Use nanmean to handle NaN values
        avg_profile = np.nanmean(far_U_1_normalized[:, start_idx:end_idx], axis=1)
        
        # Skip if the profile contains only NaN values
        if np.all(np.isnan(avg_profile)):
            continue
            
        far_wake_profiles.append(avg_profile)
        
        # Calculate approximate x1/D position (center of each segment)
        # Handle potential divide by zero
        if far_D.size > 0 and far_D[0, 0] != 0:
            pos = np.mean(far_x_1[0, start_idx:end_idx]) / far_D[0, 0]
        else:
            pos = i  # Fallback to index if D is zero or empty
        far_wake_positions.append(pos)

# Check if we have any profiles to plot
print(f"\nNumber of near wake profiles generated: {len(near_wake_profiles)}")
print(f"Number of far wake profiles generated: {len(far_wake_profiles)}")

# Only plot if we have profiles
if len(near_wake_profiles) > 0 or len(far_wake_profiles) > 0:
    # Normalize x_2 coordinates
    near_x_2_normalized = near_x_2 / near_D if near_D.size > 0 and near_D[0, 0] != 0 else near_x_2
    far_x_2_normalized = far_x_2 / far_D if far_D.size > 0 and far_D[0, 0] != 0 else far_x_2
    
    # Plot the profiles
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
else:
    print("WARNING: No valid profiles were generated. Check your data for NaN values or empty arrays.")