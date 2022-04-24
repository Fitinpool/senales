#LIBRERIAS PARA FUNCIONES
from cmath import exp
from random import triangular
import numpy as np
from scipy import signal,integrate
import matplotlib. pyplot as plt
"""
#PARTE 1
#VARIABLES
u = lambda t: np.piecewise(t,t>=0,[1,0])

b = 2 ; a = -b; dt =0.01; fs = 5

t = np.arange(a, b+dt, dt)

#periodica
sinoidal = lambda t: np.sin(t)
cuadrada = lambda t: signal.square(2 * np.pi * fs * t)
triangular = lambda t: signal.sawtooth(2 * np.pi * fs * t, 0.5)
dienteSierra = lambda t: signal.sawtooth(2 * np.pi * fs * t)

periodicas = {'SINOIDAL' : sinoidal, 'CUADRADA' : cuadrada, 'TRIANGULAR' : triangular, 'DIENTESIERRA' : dienteSierra}

i = 1
for x,y in periodicas.items():
    fig1 = plt.figure("PARTE1")
    fig1.subplots_adjust(hspace=0.5, wspace=0.5)
    fig1.suptitle('SEÑALES PERIODICAS')

    g = y(t)
    ax = fig1.add_subplot(2,2,i)
    ax.set_ylabel('impulso')
    ax.set_xlabel('t')
    ax.set_title(x)
    ax.plot(t,g)

    i += 1

plt.show()

#PARTE2
#EXPONENCIAL DECRECIENTE
#VARIABLES
u = lambda t: np.piecewise(t,t>=0,[1,0])

b = 2 ; a = -b; dt =0.01; fs = 5

t = np.arange(a, b+dt, dt)

# PROCEDIMIENTO
#aperiodica
ExpDecreciente = lambda t:(np.exp(-t) * (u(t) - u(t-1)))
ExpCreciente = lambda t: np.exp(t) * (u(t) - u(t - 1))
impulso = lambda t : u(t) - u(t - 1)
escalon = lambda t: u(t)
sinc = lambda t: np.sin(t)/t

aperiodicas = {'DECRECIENTE' : ExpDecreciente, 'CRECIENTE' : ExpCreciente, 'IMPULSO' : impulso, 'ESCALON' : escalon}

# SALIDA - GRAFICA
i = 1
for x,y in aperiodicas.items():

    h = y(t)

    fig2 = plt.figure("PARTE2")
    fig2.subplots_adjust(hspace=0.5, wspace=0.5)
    fig2.suptitle('SEÑALES APERIODICAS')

    ax = fig2.add_subplot(2,2,i)
    ax.set_ylabel('impulso')
    ax.set_xlabel('t')
    ax.set_title(x)
    ax.plot(t,h)
    i += 1

plt.show()

h = sinc(t)

fig2 = plt.figure("PARTE2")
fig2.subplots_adjust(hspace=0.5, wspace=0.5)
fig2.suptitle('SEÑALES APERIODICAS')

ax = fig2.add_subplot(2,2,1)
ax.set_ylabel('impulso')
ax.set_xlabel('t')
ax.set_title(x)
ax.plot(t,h)

plt.show()

#PARTE 3
#EXPONENCIAL DECRECIENTE
#VARIABLES
u = lambda t: np.piecewise(t,t>=0,[1,0])

b = 2 ; a = -b; dt =0.01; fs = 5

t = np.arange(a, b+dt, dt)

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


aperiodicas = {'DECRECIENTE' : ExpDecreciente, 'CRECIENTE' : ExpCreciente, 'IMPULSO' : impulso, 'ESCALON' : escalon, 'SINC' : sinc}
periodicas = {'SINOIDAL' : sinoidal, 'CUADRADA' : cuadrada, 'TRIANGULAR' : triangular, 'DIENTESIERRA' : dienteSierra}

# Integral de Convolucion x[t]*h[t]
for x,x1 in aperiodicas.items():

    xi = x1(t)
    i = 1
    for y,y1 in periodicas.items():
        hi = y1(t)

        yi = np.convolve(xi,hi,'same')*dt
        fig3 = plt.figure("PARTE3")
        fig3.subplots_adjust(hspace=0.5, wspace=0.5)
        fig3.suptitle('Convolucion ' + x + ' con señales periodicas.')

        ax = fig3.add_subplot(2,2,i)
        ax.set_ylabel('impulso')
        ax.set_xlabel('t')
        ax.set_title(x + ' convolucionada con ' + y)
        ax.plot(t,yi)
        i += 1

    plt.show()

#PARTE 4

#EXPONENCIAL DECRECIENTE
#VARIABLES
u = lambda t: np.piecewise(t,t>=0,[1,0])

b = 2 ; a = -b; dt =0.01; fs = 5

t = np.arange(a, b, dt)

# PROCEDIMIENTO
#aperiodica
ExpDecreciente = lambda t:(np.exp(-t) * (u(t) - u(t-1)))
ExpCreciente = lambda t: np.exp(t) * (u(t) - u(t - 1))
impulso = lambda t : u(t) - u(t - 1)
escalon = lambda t: u(t)
sinc = lambda t: np.sin(t)/t


aperiodicas = {'DECRECIENTE' : ExpDecreciente, 'CRECIENTE' : ExpCreciente, 'IMPULSO' : impulso, 'ESCALON' : escalon, 'SINC' : sinc}

for x,y in aperiodicas.items():
    fig4 = plt.figure("PARTE4")
    fig4.subplots_adjust(hspace=0.5, wspace=0.5)

    g = y(t)
    ax = fig4.add_subplot(2,2,1)
    ax.set_ylabel('impulso')
    ax.set_xlabel('t')
    ax.set_title(x + ' sin modificacion')
    ax.plot(t,g)

    h = 2*y(t-0.5)

    ax = fig4.add_subplot(2,2,2)
    ax.set_ylabel('impulso')
    ax.set_xlabel('t')
    ax.set_title(x + ' amplificada por 2 y movida 0.5')
    ax.plot(t,h)

    plt.show()
"""
#PARTE 5

#VARIABLES
u = lambda t: np.piecewise(t,t>=0,[1,0])

b = 2 ; a = -b; dt =0.01; fs = 5

t = np.arange(a, b, dt)

# PROCEDIMIENTO
#aperiodica
ExpDecreciente = lambda t:(np.exp(-t) * (u(t) - u(t-1)))
ExpCreciente = lambda t: np.exp(t) * (u(t) - u(t - 1))
impulso = lambda t : u(t) - u(t - 1)
escalon = lambda t: u(t)
sinc = lambda t: np.sin(t)/t


aperiodicas = {'DECRECIENTE' : ExpDecreciente, 'CRECIENTE' : ExpCreciente, 'IMPULSO' : impulso, 'ESCALON' : escalon}

i = 1
for x,y in aperiodicas.items():
    xi = sinc(t)
    hi = y(t)

    yi = np.convolve(xi,hi,'same')*dt

    print(yi)
    v, err = integrate.simps(yi.astype(float)**2, t)
    print(v)

    fig5 = plt.figure("PARTE5")
    fig5.subplots_adjust(hspace=0.5, wspace=0.5)
    fig5.suptitle('Convoluciones entre señales aperiodicas.')

    ax = fig5.add_subplot(2,2,i)
    ax.set_ylabel('impulso')
    ax.set_xlabel('t')
    ax.set_title('Señal SINC convolucionada con ' + x)
    ax.plot(t,yi)
    i += 1

plt.show()

xi = sinc(t)
hi = sinc(t)

yi = np.convolve(xi,hi,'same')*dt
fig5 = plt.figure("PARTE5")
fig5.subplots_adjust(hspace=0.5, wspace=0.5)
fig5.suptitle('Convoluciones entre señales aperiodicas.')

ax = fig5.add_subplot(2,2,1)
ax.set_ylabel('impulso')
ax.set_xlabel('t')
ax.set_title('Señal SINC convolucionada con SINC')
ax.plot(t,yi)

plt.show()