#%%

import numpy as np

# ------------------------
# PARAMETERS
# ------------------------
L = 100
steps = 100000

D = 1.0

alpha = 1.0      # drift
beta = 0.5       # Zeeman spin bias
gamma = 0.8      # EMCHA coupling
delta = 0.6      # eMChA nonlinear term

# ------------------------
# SYSTEM STATE
# ------------------------
x = L // 2
s = np.random.choice([-1, 1])

# ------------------------
# FIELDS
# ------------------------
F = 0.2
B = 1.0
chi = 1

#%%
def rates(F, s, B, chi):
    forward_exponent = (
        alpha * F
        + beta * s * B
        + gamma * chi * s * B * F
        + delta * chi * B * F**2
    )

    backward_exponent = (
        -alpha * F
        - beta * s * B
        - gamma * chi * s * B * F
        - delta * chi * B * F**2
    )

    Wf = D * np.exp(forward_exponent)
    Wb = D * np.exp(backward_exponent)

    return Wf, Wb

def step(x, s):
    Wf, Wb = rates(F, s, B, chi)

    P_total = Wf + Wb
    r = np.random.rand()

    if r < Wf / P_total:
        x += 1
    else:
        x -= 1

    return x, s

positions = []

for t in range(steps):
    x, s = step(x, s)
    positions.append(x)

positions = np.array(positions)
current = (positions[-1] - positions[0]) / steps

#%%
#F -> -F
#B -> -B

# Delta_I = I(F, B) - I(-F, B)

#%%
# I(F,B) ~ a F + b * chi * B * F**2

# V = np.array([...])
# I = np.array([...])
# coeffs = np.polyfit(V, I, 2)
# quadratic term ~ eMChA
