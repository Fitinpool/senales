#LIBRERIAS PARA FUNCIONES
import numpy as np
from scipy import signal
import matplotlib. pyplot as plt

""" Parte 1
#VARIABLES SINOIDAL
signalTime = np.arange(5, 20, 0.1); 
  
signalAmplitude = np.sin(signalTime)

#VARIABLE CUADRADA, TRIANGULAR Y DIENTE DE SIERRA
a = 0
b = 1
n = 2000
fs = 5

t = np.linspace(a, b, n, endpoint=True)

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

#comentario lol