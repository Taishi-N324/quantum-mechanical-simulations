import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib

matplotlib.use("Agg")

c = 2.99792458E+8
mu0 = 4.0 * np.pi * 1.0E-7
epsilon0 = 1.0 / (4.0 * np.pi * c ** 2) * 1.0E+7
h = 6.6260896 * 1.0E-34
hbar = h / (2.0 * np.pi)
me = 9.10938215 * 1.0E-31
eV = 1.60217733 * 1.0E-19
I = 1j
Lz, Ld = 10000, 20
dz = 1.0E-13
N = 5000
delta_z = 0.5E-8
sigma = 2.0 * np.sqrt(2.0 * np.log(2.0)) / (delta_z)
dk = 30.0 / (delta_z * (N + 1))
iz0 = Lz // 2 - Lz // 10
E0 = 10.0 * eV
k0 = np.sqrt(2.0 * me * E0 / hbar ** 2)
omega0 = hbar / (2.0 * me) * k0 ** 2
V1 = 0
V2 = 20.0 * eV
V3 = 0
d = 2000 * dz
ts, te = 0, 100
dt = 2.0 * np.pi / omega0 / (te - ts + 2)


class Matrix22Complex:
    def __init__(self):
        self.a = [[0+0j]*2 for _ in range(2)]

    def myset(self, i, j, k):
        self.a[i][j] = k

    def get(self, i, j):
        return self.a[i][j]

    def __mul__(self, c):
        temp = Matrix22Complex()
        for i in range(2):
            for j in range(2):
                temp.a[i][j] = 0
                for k in range(2):
                    temp.a[i][j] += self.a[i][k] * c.a[k][j]
        return temp

def ReflectionCoefficient(MM):
    return -(MM.get(1, 0) / MM.get(1, 1))

def TransmissionCoefficient(MM):
    return (MM.get(0, 0) * MM.get(1, 1) - MM.get(0, 1) * MM.get(1, 0)) / MM.get(1, 1)





z_dz_all = []
Psi_real_all = []
Psi_imag_all = []
Psi_ab_all = []

for tn in range(ts, te+1):
    t_real = dt * tn
    z_dz = []
    Psi_real = []
    Psi_imag = []
    Psi_ab = []
    
    for iz in range(0, Lz+1, Ld):
        z = dz * (iz - iz0)
        Psi = 0j
        k = k0
        omega = hbar / (2 * me) * k ** 2
        E = hbar * omega
        k1 = cmath.sqrt(k ** 2 - (2 * me * V1) / hbar ** 2)
        k2 = cmath.sqrt(k ** 2 - (2 * me * V2) / hbar ** 2)
        k3 = cmath.sqrt(k ** 2 - (2 * me * V3) / hbar ** 2)

        M21 = Matrix22Complex()
        T2 = Matrix22Complex()
        M32 = Matrix22Complex()
        M31 = Matrix22Complex()
        
        M21.myset(0,0,(1.0+0j + k1/k2)/2.0)
        M21.myset(0,1,(1.0+0j - k1/k2)/2.0)
        M21.myset(1,0,(1.0+0j - k1/k2)/2.0)
        M21.myset(1,1,(1.0+0j + k1/k2)/2.0)
        M32.myset(0,0,(1.0+0j + k2/k3)/2.0)
        M32.myset(0,1,(1.0+0j - k2/k3)/2.0)
        M32.myset(1,0,(1.0+0j - k2/k3)/2.0)
        M32.myset(1,1,(1.0+0j + k2/k3)/2.0)
        T2.myset(0,0,cmath.exp(I * k2 * d))
        T2.myset(0,1,0)
        T2.myset(1,0,0)
        T2.myset(1,1,cmath.exp(-I * k2 * d))
        M31 = M32 * T2 * M21

        r = ReflectionCoefficient(M31)
        t = TransmissionCoefficient(M31)
        
        if(z<0):
            Psi += (cmath.exp(I * k1 * z) + r * cmath.exp(-I * k1 * z))* cmath.exp(-I * omega * t_real)
        elif(z<d):
            Ap = M21.get(0, 0) * 1.0 + M21.get(0, 1) * r
            Am = M21.get(1, 0) * 1.0 + M21.get(1, 1) * r
            Psi += ( Ap * cmath.exp(I * k2 * z) + Am * cmath.exp(-I * k2 * z))* cmath.exp(-I * omega * t_real)
        else:
            Psi += (t * cmath.exp(I * k3 * (z - d)))* cmath.exp(-I * omega * t_real);

        z_dz.append(z/dz)
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

    plt.xlim(-4000, 6000)
    plt.ylim(-2.5, 2.5)
    plt.plot([0,0], [-2.3,2.3],color = "c")
    plt.plot([0,2000], [2.3,2.3],color = "c")
    plt.plot([2000,2000], [-2.3,2.3],color = "c")
    plt.plot([2000,2000], [-2.3,2.3],color = "c")
    plt.plot([-4000,0], [-2.3,-2.3],color = "c")
    plt.plot([2000,6000], [-2.3,-2.3],color = "c")
    plt.scatter(z_dz_all[i],Psi_real_all[i],s=15,c="y")
    plt.scatter(z_dz_all[i],Psi_imag_all[i],s=15,c="m")
    plt.scatter(z_dz_all[i],Psi_ab_all[i],s=15,c="green")
    plt.grid()
    plt.title("Tunnel effect by box-shaped potential barrier")

ani = FuncAnimation(fig, update, fargs = ('Initial Animation! ', 2.0), \
    interval = 100, frames = 100)

ani.save('box_potentials.mp4', writer="ffmpeg")


