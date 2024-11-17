import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.integrate import quad

matplotlib.use("TkAgg")

# Dados

sigma = 5.67e-8
C1 = 3.742e8
C2 = 1.439e4


# Funções

def f(lambdaT):
    return C1/(sigma*lambdaT**5 * (np.exp(C2/lambdaT) - 1))

def F(lambdaT):
    if lambdaT <= 1e-10:
        return 0
    return quad(f, 1e-10, lambdaT)[0]  # Integral numérica

def emissividade(T):
    return 0.36*F(2*T) + 0.2*(F(4*T) - F(2*T))


# Plot

e = []
T = np.linspace(500, 3000, 100)
for t in T:
    e.append(emissividade(t))

plt.style.use("seaborn-v0_8")

plt.plot(T, e, color="red")
plt.title("Emissividade em função da temperatura")
plt.xlabel("Temperatura (K)")
plt.ylabel("Emissividade")
plt.show()
