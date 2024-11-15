import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

matplotlib.use("TkAgg")

# Dados do problema

D = 0.025
L = 3
mu = 467e-6
mu_s = 855e-6
k = 0.653
h_fg = 244000
cp = 4185
rho = 770
V = 0.186
Pr = 2.99

Tinf = 27.4
Tmi = 60


# Funções

def Re(m):
    return 4*m/(math.pi*mu*D)


def Nu(m):
    return 0.027*Re(m)**(4/5)*Pr**(1/3)*(mu/mu_s)**0.14


def h(m):
    return Nu(m)*k/D


def Tmo(m):
    return Tinf - (Tinf - Tmi)*math.exp(-math.pi*D*L*h(m)/(m*cp))


def q(m):
    return m*cp*(Tmi - Tmo(m))


def t(m):
    return rho*V*h_fg/q(m)


m = np.linspace(0.1, 0.5, 100)
Q = []
T = []

for i in m:
    Q.append(q(i))
    T.append(Tmo(i))

# print(f"Para m = 0.5 kg/s, t = {t(0.5)/3600} h")

plt.style.use("seaborn-v0_8")  # Estilo do gráfico (outros em https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html)
plt.plot(m, Q, color="blue")
plt.xlabel("Vazão mássica (kg/s)")
plt.ylabel("Taxa de calor (W)")
plt.title("Comportamento da taxa de calor em função da vazão")
plt.show()

plt.plot(m, T, color="red")
plt.xlabel("Vazão mássica (kg/s)")
plt.ylabel("Temperatura de saída (ºC)")
plt.title("Comportamento da temperatura de saída em função da vazão")
plt.show()
