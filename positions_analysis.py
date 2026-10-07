#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
D = 0.2
x = np.array([11,22,44], dtype=float)
I = np.array([1.19E-02,
3.15E-03,
8.00E-04,
])

n, b = np.polyfit(np.log(x), np.log(I), 1)
n = -n

print(f"Exponent = {n:.3f}")
print(f"A = {np.exp(b):.3e}")

plt.plot(x, I, 'o', label='Data')
plt.plot(x, np.exp(b) * x**(-n), label=f'Fit: A*x^(-{n:.3f})')
plt.xlabel('x')
plt.ylabel('I')
plt.legend()
plt.show()

#%%
D = 0.5
x = np.array([3,5,11], dtype=float)
I = np.array([3.18E-01,
1.30E-01,
2.97E-02,
])

n, b = np.polyfit(np.log(x), np.log(I), 1)
n = -n

print(f"Exponent = {n:.3f}")
print(f"A = {np.exp(b):.3e}")

plt.figure(figsize=(8, 6))
plt.plot(x, I, 'o', label='Data')
plt.plot(x, np.exp(b) * x**(-n), label=f'Fit: A*x^(-{n:.3f})')
plt.xlabel('x')
plt.ylabel('I')
plt.legend()
plt.show()