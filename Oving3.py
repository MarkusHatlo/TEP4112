import numpy as np
import scipy.io as sp
import matplotlib.pyplot as plt

data = sp.loadmat("Xwire_ui")
u = np.float64(data['u'].squeeze())
v = np.float64(data['v'].squeeze())


sampling_frequency = 62500


u_prime = u - np.mean(u)
v_prime = v - np.mean(v)

mean_u_prime = np.mean(u_prime)
mean_v_prime = np.mean(v_prime)

isotropy = u_prime/v_prime

print(isotropy)

#b)
var_u = np.var(u)
var_v = np.var(v)

TKE = 1/2*(var_u+var_v)

print(TKE)

#c)

plt.figure(figsize=(8, 6))
plt.hist2d(u_prime, v_prime, bins=50, density=True, cmap='viridis')
plt.xlabel("u'_1 (m/s)")
plt.ylabel("u'_2 (m/s)")
plt.title("JPDF of u'_1 and u'_2")
plt.colorbar(label='Probability Density')
plt.show()

corr_coef = np.corrcoef(u_prime, v_prime)[0, 1]
print("Correlation coefficient between u'_1 and u'_2 =", corr_coef)