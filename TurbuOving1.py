import numpy as np
import scipy.io as sp
import matplotlib.pyplot as plt

data = sp.loadmat("Xwire_ui")
u = np.float64(data['u'].squeeze())
v = np.float64(data['v'].squeeze())

#2 a)

sampling_frequency = 62500

def mean(time,speed):
    samples = time*sampling_frequency
    return np.mean(speed[:samples])

U1_15 = mean(15,u)
U1_20 = mean(20,u)
U1_25 = mean(25,u)
U1 = mean(len(u),u)

print()
print(f"This is for 15, 20 ,25 seconds {U1_15} {U1_20} {U1_25} This is all the datapoints {U1}")
# Yes the U1 is stationary

#b)

u_prime = u - np.mean(u)
v_prime = v - np.mean(v)

mean_u_prime = np.mean(u_prime)
mean_v_prime = np.mean(v_prime)

print()
print(f"This is the mean fluctuation of u and v: {mean_u_prime} and {mean_v_prime}")


#c)

var_u = np.var(u)
var_v = np.var(v)

print()
print(f"This is the variance of u and v: {var_u} and {var_v}")

#d)

turb_int = v_prime/U1

#e)

t_start = 3 * sampling_frequency
t_stop = 4 * sampling_frequency

t = np.linspace(3, 4, t_stop - t_start)
t_u = u[t_start:t_stop] 
t_u_prime = u_prime[t_start:t_stop]


plt.figure(figsize=(10, 5))
plt.plot(t, t_u, label="Velocity (u)")
plt.plot(t,t_u_prime, label='Fluctuating component (u`)', color="g")
plt.axhline(U1, color='r', linestyle='--', label="Mean Velocity U1")
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.title("Velocity over Time")
plt.legend()
plt.show()