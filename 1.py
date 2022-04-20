#LIBRERIAS PARA FUNCIONES
from cmath import exp
from random import triangular
import numpy as np
from scipy import signal
import matplotlib. pyplot as plt

"""PARTE 1
#VARIABLES SINOIDAL
signalTime = np.arange(-5, 5, 0.1); 
signalAmplitude = np.sin(signalTime)

#VARIABLE CUADRADA, TRIANGULAR Y DIENTE DE SIERRA
a = -2
b = 2
n = 2000
fs = 5

t = np.linspace(a, b, n, endpoint=True)

fig1 = plt.figure("PARTE 1")
fig1.subplots_adjust(hspace=0.5, wspace=0.5)

#ORDEN DE LAS SEÑALES
ax = fig1.add_subplot(2,2,1)
ax.set_ylabel('Amplitud')
ax.set_xlabel('Tiempo')
ax.set_title('SINOIDAL')
ax.plot(signalTime,signalAmplitude)


ax = fig1.add_subplot(2,2,2)
ax.set_ylabel('Amplitud')
ax.set_xlabel('Tiempo')
ax.set_title('CUADRADA')
ax.plot(t, signal.square(2 * np.pi * fs * t))

ax = fig1.add_subplot(2,2,3)
ax.set_ylabel('Amplitud')
ax.set_xlabel('Tiempo')
ax.set_title('TRIANGULAR')
ax.plot(signal.sawtooth(2 * np.pi * fs * t, 0.5))

ax = fig1.add_subplot(2,2,4)
ax.set_ylabel('Amplitud')
ax.set_xlabel('Tiempo')
ax.set_title('DIENTE DE SIERRA')
ax.plot(t , signal.sawtooth(2 * np.pi * fs * t))

plt.show()
"""
"""PARTE 2
#PARTE2
#EXPONENCIAL DECRECIENTE
#VARIABLES
u = lambda t: np.piecewise(t,t>=0,[1,0])

a = -2
b = 2
dt = 0.1
t = np.arange(a, b, dt)

fig2 = plt.figure("PARTE2")
fig2.subplots_adjust(hspace=0.5, wspace=0.5)


# PROCEDIMIENTO

ExpDecreciente = np.exp(-t) * (u(t) - u(t - 1))

# SALIDA - GRAFICA
ax = fig2.add_subplot(2,2,1)
ax.set_ylabel('impulso')
ax.set_xlabel('t')
ax.set_title('EXPONENCIAL DECRACIENTE')
ax.plot(t,ExpDecreciente)

#EXPONENCIAL CRECIENTE
# PROCEDIMIENTO

ExpCreciente = np.exp(t) * (u(t) - u(t - 1))

# SALIDA - GRAFICA
ax = fig2.add_subplot(2,2,2)
ax.set_ylabel('impulso')
ax.set_xlabel('t')
ax.set_title('EXPONENCIAL CRECIENTE')
ax.plot(t,ExpCreciente)

#IMPULSO  
# PROCEDIMIENTO

impulso = u(t) - u(t - 1)

# SALIDA - GRAFICA
ax = fig2.add_subplot(2,2,3)
ax.set_ylabel('impulso')
ax.set_xlabel('t')
ax.set_title('IMPULSO')
ax.plot(t,impulso)

#ESCALON
# PROCEDIMIENTO

escalon = u(t)

# SALIDA - GRAFICA
ax = fig2.add_subplot(2,2,4)
ax.set_ylabel('impulso')
ax.set_xlabel('t')
ax.set_title('ESCALON')
ax.plot(t,escalon)
plt.show()
#Sinc(x)
# PROCEDIMIENTO

sinc = np.sin(t)/t

fig3 = plt.figure("PARTE 2")
fig3.subplots_adjust(hspace=0.5, wspace=0.5)

# SALIDA - GRAFICA
ax = fig3.add_subplot(2,2,1)
ax.set_ylabel('impulso')
ax.set_xlabel('t')
ax.set_title('SINC')
ax.plot(t,sinc)

plt.show()
"""
"""PARTE 3
#PARTE 3
# INGRESO
# tiempo [a,b] simétrico alrededor de 0
b = 2 ; a = -b; dt =0.01; fs = 5

u = lambda t: np.piecewise(t,t>=0,[1,0])

x = lambda t: signal.square(2 * np.pi * fs * t)
h = lambda t: np.sin(t)/t


# PROCEDIMIENTO
ti = np.arange(a,b+dt,dt)
xi = x(ti)
hi = h(ti)

# Integral de Convolucion x[t]*h[t]
# corrección de magnitud por dt para en integral
yi = np.convolve(xi,hi,'same')*dt

# SALIDA - GRAFICA
plt.figure('PARTE 3')
plt.suptitle('Convolucion')

plt.subplot(211)
plt.plot(ti,xi,'b', label='x(t)')
plt.plot(ti,hi,'r', label='h(t)')
plt.legend()
plt.grid()

plt.subplot(212)
plt.plot(ti,yi,'m', label='x(t)*h(t)')
plt.xlabel('t')
plt.legend()
plt.grid()

plt.show()
"""
"""PARTE 4
#PARTE 4

#EXPONENCIAL DECRECIENTE
#VARIABLES
u = lambda t: np.piecewise(t,t>=0,[1,0])

b = 2 ; a = -b; dt =0.01; fs = 5

# PROCEDIMIENTO
#aperiodica
ExpDecreciente = lambda t:(np.exp(-t) * (u(t) - u(t-1)))
ExpCreciente = lambda t: np.exp(t) * (u(t) - u(t - 1))
impulso = lambda t : u(t) - u(t - 1)
escalon = lambda t: u(t)
sinc = lambda t: np.sin(t)/t

#periodica
sinoidal = lambda t: np.sin(t)
cuadrada = lambda t: signal.square(2 * np.pi * fs * t)
triangular = lambda t: signal.sawtooth(2 * np.pi * fs * t, 0.5)
dienteSierra = lambda t: signal.sawtooth(2 * np.pi * fs * t)


aperiodicas = {'DECRECIENTE' : ExpDecreciente, 'CRECIENTE' : ExpCreciente, 'IMPULSO' : impulso, 'ESCALON' : escalon, 'SINC', sinc}
periodicas = {'SINOIDAL' : sinoidal, 'CUADRADA' : cuadrada, 'TRIANGULAR' : triangular, 'DIENTESIERRA' : dienteSierra}

i = 1
for x,y in aperiodicas.items():
    fig4 = plt.figure("PARTE2")
    fig4.subplots_adjust(hspace=0.5, wspace=0.5)

    g = y(t)
    ax = fig4.add_subplot(2,2,2)
    ax.set_ylabel('impulso')
    ax.set_xlabel('t')
    ax.set_title(x + ' sin modificacion')
    ax.plot(t,g)

    h = 2*y(t-0.5)

    ax = fig4.add_subplot(2,1,2)
    ax.set_ylabel('impulso')
    ax.set_xlabel('t')
    ax.set_title(x + ' amplificada por 2 y movida 0.5')
    ax.plot(t,h)

    plt.show()
"""

"""PARTE 5
#PARTE 5

#SEÑALES PERIODICAS

sinoidal = lambda t: np.sin(t)
cuadrada = lambda t: signal.square(2 * np.pi * fs * t)
triangular = lambda t: signal.sawtooth(2 * np.pi * fs * t, 0.5)
dienteSierra = lambda t: signal.sawtooth(2 * np.pi * fs * t)
"""