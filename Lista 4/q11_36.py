import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

matplotlib.use("TkAgg")


# Dados

D = 0.02
L = 2
N = 100  # número de tubos
Th = 373
Tc_i = 290
h_h = 10000
A = N*math.pi*D*2*L
cp_c = 4180
k = 0.64
mu = 577e-6
Pr = 3.77


# Funções

def Re(m):
    return 4*(m/N)/(math.pi*D*mu)

def h_c(m):
    return 0.023 * Re(m)**0.8 * Pr**0.4 * (k/D)

def U(m, Rf):
    return 1/(1/h_h + 1/h_c(m) + Rf)

def NUT(m, Rf):
    return U(m, Rf)*A/(m*cp_c)

def epsilon(m, Rf):
    return 1 - math.exp(-NUT(m, Rf))

def q(m, Rf):
    return epsilon(m, Rf)*(m*cp_c)*(Th - Tc_i)

def tc_o(m, Rf):
    return Tc_i + q(m, Rf)/(m*cp_c)


# Plots

m_c = np.linspace(5, 20, 100)
Tc_o = [[], [], []]
Rf = [0, 2e-4, 5e-4]
for i in range(len(Rf)):
    for m in m_c:
        Tc_o[i].append(tc_o(m, Rf[i]))

colors = ["red", "green", "blue"]

plt.style.use("seaborn-v0_8")

for i in range(len(Tc_o)):
    plt.plot(m_c, Tc_o[i], color=f"{colors[i]}", label=f"Rf = {Rf[i]} m²K/W")

plt.xlabel("Vazão de água (kg/s)")
plt.ylabel("Temperatura de saída da água (K)")
plt.title("Temperatura da água em função da vazão e da incrustação")
plt.legend()
plt.show()
