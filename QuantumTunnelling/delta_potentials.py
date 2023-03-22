from tqdm import tqdm
import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib

matplotlib.use("Agg")

c = 2.99792458E+8
mu0 = 4.0 * np.pi * 1.0E-7
epsilon0 = 1.0 / (4.0 * np.pi * c * c) * 1.0E+7
h = 6.6260896 * 1.0E-34
hbar = h / (2.0 * np.pi)
me = 9.10938215 * 1.0E-31
eV = 1.60217733 * 1.0E-19
I = 1j
Lz = 5000
Ld = 10
dz = 1.0E-11
N = 5000
delta_z = 0.5E-8
sigma = 2.0 * np.sqrt(2.0 * np.log(2.0)) / (delta_z)
dk = 30.0 / (delta_z * float(N + 1))
iz0 = Lz // 2
E0 = 10.0 * eV
k0 = np.sqrt(2.0 * me * E0 / (hbar ** 2))
omega0 = hbar / (2.0 * me) * (k0 ** 2)
V1 = 0
V2 = E0
ts = -150
te = 300
dt = 1.0 * 1.0E-16

z_dz_all = []
Psi_real_all = []
Psi_imag_all = []
Psi_ab_all = []


for tn in range(ts, te):
    t_real = dt * tn
    z_dz = []
    Psi_real = []
    Psi_imag = []
    Psi_ab = []
        
    for iz in range(0,Lz,Ld):
        z = dz * (iz - iz0)
        Psi = 0

        for jz in range(N):
            k = (k0 + dk * (jz - N / 2))
            omega = hbar / (2.0*me) * k ** 2
            E = hbar * omega
            k1 = cmath.sqrt(2.0*me*(E - V1)) / hbar
            if (E > V2):
                k2 = cmath.sqrt(2.0*me*(E - V2)) / hbar
            else:
                k2 = I * cmath.sqrt(2.0*me*(V2 - E)) / hbar
            r = (k1 - k2) / (k1 + k2)
            t = (2.0 * k1) / (k1 + k2)

            if (z < 0):
                Psi += (cmath.exp(I * k1 * z)+r* cmath.exp(-I * k1 * z))* cmath.exp(-I * omega * t_real) * cmath.exp(-1.0 / 2.0 * ((k - k0) / sigma ) ** 2)
            else:
                Psi += (t * cmath.exp(I * k2 * z))* cmath.exp(-I * omega * t_real) * cmath.exp(-1.0 / 2.0 * ((k - k0) / sigma ) ** 2)

        Psi = Psi / N
        z_dz.append(z / dz )
        Psi_real.append(Psi.real)
        Psi_imag.append(Psi.imag)
        Psi_ab.append(abs(Psi))

    z_dz_all.append(z_dz)
    Psi_real_all.append(Psi_real)
    Psi_imag_all.append(Psi_imag)
    Psi_ab_all.append(Psi_ab)        

fig = plt.figure(figsize = (10, 6))

def update(i, fig_title, A):
    if i != 0:
        plt.cla()                     

    plt.xlim(-2500, 2500)
    plt.ylim(-0.30, 0.30)
    plt.plot([0,0], [-0.30,0.30],color = "c")
    plt.scatter(z_dz_all[i],Psi_real_all[i],s=15,c="gold")
    plt.scatter(z_dz_all[i],Psi_imag_all[i],s=15,c="blueviolet")
    plt.scatter(z_dz_all[i],Psi_ab_all[i],s=15,c="green")
    plt.grid()
    plt.title("Collisions on delta-functional potentials")

ani = FuncAnimation(fig, update, fargs = ('Initial Animation! ', 2.0), \
    interval = 100, frames = 450)

ani.save('delta_potentials.mp4', writer="ffmpeg", dpi=300)

