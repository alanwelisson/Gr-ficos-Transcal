import math
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

matplotlib.use("TkAgg")

# Dados

T1 = 500
T2 = T1
D = 0.4
L = 0.2
e1 = 0.6
ed2 = 0.8
Tinf = 300
Tviz = Tinf
sigma = 5.67e-8
g = 9.81
Ad = math.pi*(D**2)/4

# Fatores de forma
F13 = 0.382
F12 = 0.618
F31 = F13
F32 = F12

J2 = sigma*T2**4


# Propriedades (Tabela do Incropera)

nu = {
    300 : 15.89e-6,
    350 : 20.92e-6,
    400 : 26.41e-6,
    450 : 32.39e-6,
    500 : 38.79e-6
}

k = {
    300 : 26.3e-3,
    350 : 30.0e-3,
    400 : 33.8e-3,
    450 : 37.3e-3,
    500 : 40.7e-3
}

rho = {
    300 : 1.1614,
    350 : 0.9950,
    400 : 0.8711,
    450 : 0.7740,
    500 : 0.6964
}

cp = {
    300 : 1007,
    350 : 1009,
    400 : 1014,
    450 : 1021,
    500 : 1030
}


# Funções

def interpolacao(prop, T):
    temperaturas = list(prop.keys())
    for i in range(len(temperaturas)):
        if T <= temperaturas[i]:
            T_sup = temperaturas[i]
            T_inf = temperaturas[i-1]
            return prop[T_inf] + (prop[T_sup] - prop[T_inf])/50 * (T - T_inf)

def q_rad(T):
    return ed2*sigma*Ad*(T**4 - Tviz**4)

def RaL(T):
    Tf = (T + Tinf)/2
    nu_ar = interpolacao(nu, Tf)
    alfa = interpolacao(k, Tf)/(interpolacao(rho, Tf)*interpolacao(cp, Tf))
    beta = 1/Tf
    return (g*beta*(T-Tinf)*0.1**3)/(nu_ar*alfa)

def h(T):
    Tf = (T + Tinf)/2
    k_ar = interpolacao(k, Tf)
    return (k_ar/0.1) * 0.54 * RaL(T)**(1/4)

def q_conv(T):
    return h(T)*Ad*(T - Tinf)

def Jd(T):
    return sigma*T**4

def J1(T):
    return ((e1*sigma*T1**4)/(1-e1) + F12*J2 + F13*Jd(T))/(e1/(1-e1) + 1)

def q_d(T):
    return Ad * (F31*(Jd(T) - J1(T)) + F32*(Jd(T) - J2))

def q_liq(T):
    return -(q_conv(T) + q_rad(T) + q_d(T))


# Estado estacionário (solução para q_lid(T) = 0)

T_ee = fsolve(q_liq, 400)[0]
print(f"Temperatura do disco no estado estacionário: {round(T_ee, 1)} K")


# Plot

T = np.linspace(300, 500, 200)
q = []
for t in T:
    q.append(q_liq(t))

plt.style.use("seaborn-v0_8")
plt.plot(T, q, color="red")
plt.xlabel("Temperatura (K)")
plt.ylabel("Taxa líquida de transferência de calor (W)")
plt.title("Variação da taxa líquida com a temperatura do disco")
plt.show()

