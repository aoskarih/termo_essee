import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np


plt.rcParams.update({'font.size': 16})

hbar = 5.304e-65
me = 4.581e-61
mp = 8.412e-58
G = 1.327e20
c = 2.998e8
chan_limit = 1.72

def func1A(m):
    return 3/20 *(9*np.pi/8)**(2/3)*hbar**2/me *(m/mp)**(5/3)

def funcB(m):
    return 3/5 * G * m**2


def func2A(m):
    return 3/8 * (9 * np.pi/8)**(1/3) * hbar*c*(m/mp)**(4/3)

def func2C(m):
    return 3/4 * 1/(9 * np.pi)**(1/3) * me**2 * c**3 /(hbar) * (m/mp)**(2/3)


def funcNewt(r):
    return 7.004e20 * r**(-3)

def funcReal(r):
    a = 2.325e-8
    b = 0.4817
    c = 7.277e-15
    d = 0.644
    return (a*r+b)/(np.exp(c*r**2)-d)



M = np.linspace(0, 2, 1000)
R = np.linspace(0.05e7, 4e7, 1000)

#s, v = curve_fit(func, t, r)
#print(s)

plt.figure(1)
plt.xlim(0, 2)
plt.ylim(0, 3)
plt.xlabel("$\mathrm{Massa} \; (\mathrm{M}_{\odot})$")
plt.ylabel("$\mathrm{Säde}\; (\mathrm{m} \cdot 10^7)$")

plt.plot(M, (2*func1A(M)/funcB(M))/1e7, label="$p \ll m c$")
plt.plot(M, (np.sqrt((func2A(M)-funcB(M))/(func2C(M))))/1e7, label="$p \gg m c$")
plt.plot([chan_limit, chan_limit], [0, 1e3], "--", label="$\mathrm{Chandrasekharin\; raja}$")
plt.legend()


plt.figure(2)
plt.xlim(0, 1.5)
plt.ylim(0, 2)
plt.plot(funcReal(R), R/1e7, "r", label="$\mathrm{TOV}$")
plt.plot(funcNewt(R), R/1e7, "--", label="$\mathrm{Klassinen}$")
plt.xlabel("$\mathrm{Massa} \; (\mathrm{M}_{\odot})$")
plt.ylabel("$\mathrm{Säde}\; (\mathrm{m} \cdot 10^7)$")

plt.legend()
plt.show()
