import math
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

matplotlib.use("TkAgg")

# Dados (p = 0.1 bar)

D = 0.008
Tsat = 320
rho_v = 0.072
h_fg = 2390e3
g = 9.81


# Propriedades para 280 <= Ts <= 300

rho_l = {
    280: 1/(1e-3),
    285: 1/(1e-3),
    290: 1/(1.001e-3),
    295: 1/(1.002e-3),
    300: 1/(1.003e-3)
}

cp_l = {
    280: 4198,
    285: 4189,
    290: 4184,
    295: 4181,
    300: 4179
}

mu_l = {
    280: 1422e-6,
    285: 1225e-6,
    290: 1080e-6,
    295: 959e-6,
    300: 855e-6
}

k_l = {
    280: 582e-3,
    285: 590e-3,
    290: 598e-3,
    295: 606e-3,
    300: 613e-3
}


# Interpolação das propriedades

def interpolacao(prop, T):
    temperaturas = list(prop.keys())
    for i in range(len(temperaturas)):
        if T <= temperaturas[i]:
            T_inf = temperaturas[i-1]
            T_sup = temperaturas[i]
            return prop[T_inf] + ((prop[T_sup]-prop[T_inf])/5) * (T - T_inf)
        else:
            continue


def Rho_l(T):
    return interpolacao(rho_l, T)


def K_l(T):
    return interpolacao(k_l, T)


def Cp_l(T):
    return interpolacao(cp_l, T)


def Mu_l(T):
    return interpolacao(mu_l, T)


def h_fg_linha(T):
    return h_fg + 0.68*Cp_l(T)*(Tsat - T)


# Cálculo da taxa de condensação

def h(N, T):
    return 0.729 * (g*Rho_l(T)*(Rho_l(T) - rho_v)*K_l(T)**3 * h_fg_linha(T)/(N*Mu_l(T)*(Tsat-T)*D))**(1/4)


def m(N, T):
    return h(N, T)*math.pi*D*(Tsat-T)/h_fg_linha(T)


Ts = np.linspace(280, 300, 100)
taxas = [[], [], []]
N = [2, 5, 10]

for i in range(len(N)):
    for t in Ts:
        taxas[i].append(m(N[i], t))


# Plot

colors = ["red", "green", "blue"]  # Lista de cores

plt.style.use("seaborn-v0_8")  # Estilo do gráfico (outros em https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html)
for i in range(len(taxas)):
    plt.plot(Ts, taxas[i], color=f"{colors[i]}", label=f"N = {N[i]}")
plt.xlabel("Temperatura da superfície (K)")
plt.ylabel("Taxa de condensação por comprimento (kg/s.m)")
plt.title("Influência do arranjo dos tubos na condensação")
plt.legend()
plt.show()
