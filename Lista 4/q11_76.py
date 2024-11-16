import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

matplotlib.use("TkAgg")

# Dados

D = 0.01
L = 8.11
N = 100
m_h = 1.5
cp_c = 4178
h_fg = 2304e3
k = 0.628
mu = 700e-6
Pr = 4.6
h_e = 5000
Rf = 3e-4
Th = 355
A = math.pi*D*L*N


# Funções

def Re(m):
    return 4*(m/N)/(math.pi*D*mu)

def h_i(m):
    return (k/D) * 0.023 * Re(m)**0.8 * Pr**0.4

def U(m):
    return 1/(1/h_e + 1/h_i(m) + Rf)

def NUT(m):
    return U(m)*A/(m*cp_c)

def epsilon(m):
    return 1 - math.exp(-NUT(m))

def q(m, Tc_e):
    return epsilon(m)*m*cp_c*(Th - Tc_e)

def taxa_cond(m, Tc_e):
    return q(m, Tc_e)/h_fg


# Plot

Tc_e = [275, 280, 285]  # Temperaturas de entrada (coloque quantas quiser)
m_c = np.linspace(15, 25)
m_cd = []
for i in range(len(Tc_e)):
    m_cd.append([])

for i in range(len(Tc_e)):
    for m in m_c:
        m_cd[i].append(taxa_cond(m, Tc_e[i]))

# colors = ["red", "green", "blue"] -> lista de cores

plt.style.use("seaborn-v0_8")
for i in range(len(m_cd)):
    plt.plot(m_c, m_cd[i], label=f"Tc_e = {Tc_e[i]} K")  # Adicionar color=f"{colors[i]}" se tiver a lista de cores
plt.xlabel("Vazão de água (kg/s)")
plt.ylabel("Taxa de condensação (kg/s)")
plt.title("Influência da temperatura e da vazão de água na condensação")
plt.legend()
plt.show()

