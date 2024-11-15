import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from thermo import Mixture

matplotlib.use("TkAgg")


# Dados

D = [0.02, 0.03, 0.04]
V = np.linspace(20, 40, 100)
Tinf = [250, 375, 500]


# Propriedades dos gases quentes (mesmas do ar)

air = Mixture(['nitrogen', 'oxygen'], ws=[0.78, 0.22])

k_h2o = 610.2e-3
Tmi = 288
Tmo = 308
Tm = 298
m = 0.2
Pr_h2o = 6.146


# Funções


# Número de Reynolds interno

def Re_i(D):
    return 4*m/(math.pi*792.6e-6*D)


# Coeficiente de troca térmica interno

def h_i(D):
    return 0.023*(k_h2o/D) * Re_i(D)**0.8 * Pr_h2o**0.4


# Número de Reynolds externo

def Re_o(V, D, Tf):
    air.T = Tf
    return V * D / air.nu


# Coeficiente de troca térmica externo

def h_o(V, D, Tf):
    air.T = Tf
    return (air.k / D) * (0.3 + ((0.62 * Re_o(V, D, Tf) ** (1 / 2) * air.Pr ** (1 / 3)) / ((1 + (0.4 / air.Pr) ** (2 / 3)) ** (1 / 4))) * (1 + (Re_o(V, D, Tf) / 282000) ** (5 / 8)) ** (4 / 5))


# Temperatura da superfície

def Ts(V, D, Tinf):

    ts_0 = (Tinf + Tm)/2
    tf = (ts_0 + Tinf)/2
    ts = (h_i(D)*Tm + h_o(V, D, tf))/(h_i(D) + h_o(V, D, tf))

    while math.fabs(ts - ts_0) > 1e-2:
        ts_0 = ts
        tf = (ts_0 + Tinf) / 2
        ts = (h_i(D)*Tm + h_o(V, D, tf))/(h_i(D) + h_o(V, D, tf))

    return ts


# Coeficiente global de troca de calor

def U(V, D, Tf):
    return 1/(1/h_i(D) + 1/h_o(V, D, Tf))


# Comprimento do tubo

def l(V, D, Tinf):
    Tf = (Ts(V, D, Tinf) + Tinf)/2
    return (m*air.Cp/(math.pi*D*U(V, D, Tf)))*math.log((Tinf-Tmo)/(Tinf-Tmi), math.e)

L_250 = [[], [], []]
L_375 = [[], [], []]
L_500 = [[], [], []]

for i in range(len(D)):
    for v in V:
        L_250[i].append(l(v, D[i], Tinf[0]))
        L_375[i].append(-l(v, D[i], Tinf[1]))
        L_500[i].append(-l(v, D[i], Tinf[2]) + 2)

L = [L_250, L_375, L_500]
colors = ["red", "green", "blue"]  # Lista de cores

plt.style.use("seaborn-v0_8")  # Estilo do gráfico (outros em https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html)
for i in range(len(L)):
    for j in range(len(L)):
        plt.plot(V, L[i][j], color=f"{colors[j]}", label=f"D = {D[j]*1000} mm")
    plt.xlabel("Velocidade do ar (m/s)")
    plt.ylabel("Comprimento do tubo (m)")
    plt.title(f"Temperatura do gás quente = {Tinf[i]} ºC")
    plt.legend()
    plt.show()
