from scipy.optimize import fsolve
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

matplotlib.use("TkAgg")


# Dados

nu = 14.82e-6
alfa = 20.92e-6
k_ar = 0.0253
k_iso = 0.03
eps = 0.6
A = 0.65
Tsi = 5
Tinf = 25
sigma = 5.67e-8
g = 9.81
Pr = 0.71
beta = 0.00347


# Funções

# Número de Rayleigh
def Ra(T):
    return g*beta*(Tinf - T)/(alfa*nu)

# Coeficiente de troca térmica
def h(T):
    return k_ar * (0.825 + 0.387*Ra(T)**(1/6) / (1 + (0.497/Pr)**(9/16))**(8/27))**2


L = np.linspace(0, 0.025, 100)
Tse = []
for l in L:
    if l == 0:
        Tse.append(Tsi)
    else:
        def temp(T):
            return h(T)*A*(Tinf - T) + eps*sigma*A*((Tinf+273)**4 - (T+273)**4) - (k_iso*A/l)*(T - Tsi)
        Tse.append(float(fsolve(temp, [1, 1])[0]))


Q = []
for i in range(len(L)):
    if L[i] == 0:
        Q.append(h(Tsi)*A*(Tinf - Tsi) + eps*sigma*A*((Tinf+273)**4 - (Tsi+273)**4))
    else:
        Q.append(k_iso*A*(Tse[i]-Tsi)/L[i])


# Plot

plt.style.use("seaborn-v0_8")  # Estilo do gráfico (outros em https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html)
plt.plot(L, Tse, color="red")
plt.xlabel("Espessura do isolante (m)")
plt.ylabel("Temperatura externa (ºC)")
plt.title("Temperatura da superfície externa")
plt.show()

plt.plot(L, Q, color="blue")
plt.xlabel("Espessura do isolante (m)")
plt.ylabel("Taxa de calor (W)")
plt.title("Ganho de calor devido ao isolante")
plt.show()
