#LIBRERIAS PARA FUNCIONES
import numpy as np
from scipy import signal
import matplotlib. pyplot as plt

#VARIABLES SINOIDAL
signalTime = np.arange(5, 20, 0.1); 
  
signalAmplitude = np.sin(signalTime)

#VARIABLE CUADRADA, TRIANGULAR Y DIENTE DE SIERRA
a = -2
b = 2
n = 2000
fs = 5

#t = np.linspace(a, b, n, endpoint=True)
T = 0.01
t = np.arange(-2,2 +T, T)

fig1 = plt.figure("SEÑALES")
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
ax.plot(t, signal.square(np.pi*t))

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
"""Parte 2
#EXPONENCIAL DECRECIENTE
#VARIABLES
u = lambda t: np.piecewise(t,t>=0,[1,0])

a = -5
b = 5
dt = 0.1
t = np.arange(a, b, dt)

# PROCEDIMIENTO

ExpDecreciente = np.exp(-t) * (u(t) - u(t - 1))

# SALIDA - GRAFICA
plt.figure(2)
plt.plot(t, ExpDecreciente)

plt.title('EXPONENCIAL DECRECIENTE')
plt.xlabel('t')
plt.ylabel('impulso')
plt.margins(dt)
plt.grid()
plt.show()

#EXPONENCIAL CRECIENTE
# PROCEDIMIENTO

ExpCreciente = np.exp(t) * (u(t) - u(t - 1))

# SALIDA - GRAFICA
plt.figure(2)
plt.plot(t,ExpCreciente)

plt.title('EXPONENCIAL CRECIENTE')
plt.xlabel('t')
plt.ylabel('impulso')
plt.margins(dt)
plt.grid()
plt.show()

#IMPULSO  
# PROCEDIMIENTO

impulso = u(t) - u(t - 1)

# SALIDA - GRAFICA
plt.figure(2)
plt.plot(t,impulso)

plt.title('IMPULSO')
plt.xlabel('t')
plt.ylabel('impulso')
plt.margins(dt)
plt.grid()
plt.show()

#ESCALON
# PROCEDIMIENTO

escalon = u(t)

# SALIDA - GRAFICA
plt.figure(2)
plt.plot(t,escalon)

plt.title('ESCALON')
plt.xlabel('t')
plt.ylabel('impulso')
plt.margins(dt)
plt.grid()
plt.show()

#Sinc(x)
# PROCEDIMIENTO

sinc = np.sin(t)/t

# SALIDA - GRAFICA
plt.figure(2)
plt.plot(t,sinc)

plt.title('SINC')
plt.xlabel('t')
plt.ylabel('impulso')
plt.margins(dt)
plt.grid()
plt.show()
"""

"""
# INGRESO
# tiempo [a,b] simétrico alrededor de 0
b = 5 ; a = -b; dt =0.01

u = lambda t: np.piecewise(t,t>=0,[1,0])

x = lambda t: 10*np.exp(-3*t)*u(t)
h = lambda t: (2*np.exp(-2*t)-np.exp(-t))*u(t)

# PROCEDIMIENTO
a1 = 0
b2 = 1
n1 = 10
fs1 = 5

t = np.linspace(a1, b2,1000, endpoint=True)


# PROCEDIMIENTO
ti = np.arange(a, b,dt); 
xi = signal.square(2 * np.pi * fs1 * ti)
hi = np.sin(ti)/ti

# Integral de Convolucion x[t]*h[t]
# corrección de magnitud por dt para en integral
yi = np.convolve(xi,hi,'same')*dt

# SALIDA - GRAFICA
plt.figure(1)
plt.suptitle('Integral de Convolución x(t)*h(t)')

plt.subplot(211)
plt.plot(ti,xi,'b', label='x(t)')
plt.plot(ti,hi,'r', label='h(t)')
plt.margins(dt)
plt.legend()
plt.grid()

plt.subplot(212)
plt.plot(ti,yi,'m', label='x(t)*h(t)')
plt.xlabel('t')
plt.legend()
plt.grid()

plt.show()
"""