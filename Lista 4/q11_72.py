import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

matplotlib.use("TkAgg")


# Dados (propriedades a T = 350 K)

D = 0.07
L = 2
N = 1200
NT = 30
ST = 0.14
A = math.pi*D*L*N
rho = 0.995
cp = 1009
nu = 20.92e-6
k = 0.03
Pr = 0.7
Pr_s = 0.69
h_fg = 2183e3  # Ts = 400 K
Ts = 400
Tc_i = 300


# Funções

def V(m):
    return m/(rho*NT*ST*L)

def Re_max(m):
    return 2*V(m)*D/nu

def h_o(m):
    return (k/D) * 0.27 * Re_max(m)**0.63 * Pr**0.36 * (Pr/Pr_s)**(1/4)

def NUT(m):
    return h_o(m)*A/(m*cp)

def epsilon(m):
    return 1 - math.exp(-NUT(m))

def Tc_o(m):
    return Tc_i + epsilon(m)*(Ts - Tc_i)

def q(m):
    return epsilon(m)*m*cp*(Ts - Tc_i)

def taxa_cond(m):
    return q(m)/h_fg


# Plots

tc_o = []
Q = []
m_cd = []
m_c = np.linspace(10, 50, 100)

for m in m_c:
    tc_o.append(Tc_o(m))
    Q.append(q(m))
    m_cd.append(taxa_cond(m))

plt.style.use("seaborn-v0_8")

# Temperatura de saída
plt.plot(m_c, tc_o, color="red")
plt.xlabel("Vazão de ar (kg/s)")
plt.ylabel("Temperatura de saída do ar (K)")
plt.title("Temperatura do ar em função da vazão de ar")
plt.show()

# Taxa de calor
plt.plot(m_c, Q, color="green")
plt.xlabel("Vazão de ar (kg/s)")
plt.ylabel("Taxa de calor (W)")
plt.title("Taxa de calor em função da vazão de ar")
plt.show()

# Taxa de condensação
plt.plot(m_c, m_cd, color="blue")
plt.xlabel("Vazão de ar (kg/s)")
plt.ylabel("Taxa de condensação (kg/s)")
plt.title("Taxa de condensação em função da vazão de ar")
plt.show()
