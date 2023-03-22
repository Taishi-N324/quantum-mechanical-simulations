import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
import sys
import math

def hpn_calc(x, n):
    """
    Calculate the value of the Hermite polynomial of degree n at x.
    """
    if n == 0:
        return 1
    elif n == 1:
        return 2 * x
    else:
        memo = [0] * (n+1)
        memo[0] = 1
        memo[1] = 2 * x

        for i in range(2, n+1):
            memo[i] = 2 * x * memo[i-1] - 2 * (i-1) * memo[i-2]
        return memo[n]
    
def hpn(x,n):
    return hpn_calc(x,n)* cmath.exp(- x ** 2 / 2) / cmath.sqrt(2 ** n * math.factorial(n) * cmath.sqrt(np.pi))


def psiOscillator(x,t,n):
    return (me * omega / hbar) ** 0.25 * hpn(cmath.sqrt(me * omega / hbar) * x, n) * cmath.exp(-1j * omega *(n+0.5) * t)




# matplotlib.use("Agg")

h = 6.6260896 * 1.0E-34
hbar = h / (2.0 * np.pi)
me = 9.10938215 * 1.0E-31
# omega = 1.0E-15
omega = 1.0
    
wave_functions = []
x_coordinates = []
for i in range(10):
    wave_function = []
    x_coordinate = []
    for x_coor in range(-1000, 1000):
        wave_function.append(abs(psiOscillator(x_coor/10000, 0, i))**2)
        x_coordinate.append(x_coor/10000)
    wave_functions.append(wave_function)
    x_coordinates.append(x_coordinate)
        
fig = plt.figure(figsize = (10, 6))

for i in range(10):
    plt.scatter(x_coordinates[i], wave_functions[i], s=10,label="n={}".format(i))
plt.grid()
plt.legend()
plt.savefig("existence.png", format="png", dpi=300)
plt.show()


