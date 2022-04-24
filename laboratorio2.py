#LIBRERIAS PARA FUNCIONES
from cmath import exp
from random import triangular
import numpy as np
from scipy import signal,integrate
import matplotlib. pyplot as plt

#Variables generales
u = lambda t: np.piecewise(t,t>=0,[1,0])

b = 2 ; a = -b; dt =0.01; fs = 5

t = np.arange(a, b+dt, dt)

T  = 2*np.pi/1

#periodicaa
sinoidal = lambda t: np.sin(t)
cuadrada = lambda t: signal.square(2 * np.pi * fs * t)
triangular = lambda t: signal.sawtooth(2 * np.pi * fs * t, 0.5)
dienteSierra = lambda t: signal.sawtooth(2 * np.pi * fs * t)

#aperiodica
ExpDecreciente = lambda t:(np.exp(-t) * (u(t) - u(t-1)))
ExpCreciente = lambda t: np.exp(t) * (u(t) - u(t - 1))
impulso = lambda t : u(t) - u(t - 1)
escalon = lambda t: u(t)
sinc = lambda t: np.sin(t)/t

aperiodicas = {'DECRECIENTE' : ExpDecreciente, 'CRECIENTE' : ExpCreciente, 'IMPULSO' : impulso, 'ESCALON' : escalon, 'SINC' : sinc}
aperiodicas1 = {'DECRECIENTE' : ExpDecreciente, 'CRECIENTE' : ExpCreciente, 'IMPULSO' : impulso, 'ESCALON' : escalon}
periodicas = {'SINOIDAL' : sinoidal, 'CUADRADA' : cuadrada, 'TRIANGULAR' : triangular, 'DIENTESIERRA' : dienteSierra}

#PARTE 1

i = 1
print('----------------------------------')
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

    energia = integrate.simps(g**2,t)
    potencia = energia*(1/(b-a))

    print('Energia de ' + x + ': ', energia)
    print('Potencia de ' + x + ': ', potencia)

    i += 1
print('----------------------------------')
plt.show()

#PARTE2

i = 1
for x,y in aperiodicas1.items():

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

    energia = integrate.simps(h**2,t)
    potencia = energia*(1/(b-a))

    print('Energia de ' + x + ': ', energia)
    print('Potencia de ' + x + ': ', potencia)

plt.show()

h = sinc(t)

fig2 = plt.figure("PARTE2")
fig2.subplots_adjust(hspace=0.5, wspace=0.5)
fig2.suptitle('SEÑALES APERIODICAS')

ax = fig2.add_subplot(2,2,1)
ax.set_ylabel('impulso')
ax.set_xlabel('t')
ax.set_title('SINC')
ax.plot(t,h)

energia = integrate.simps(h**2,t)
potencia = energia*(1/(b-a))

print('Energia de SINC: ', energia)
print('Potencia de SINC: ', potencia)
plt.show()

print('----------------------------------')

#PARTE 3

for x,x1 in aperiodicas.items():

    xi = x1(t)
    i = 1
    for y,y1 in periodicas.items():
        hi = y1(t)

        yi = np.convolve(hi,xi,'same')*dt
        fig3 = plt.figure("PARTE3")
        fig3.subplots_adjust(hspace=0.5, wspace=0.5)
        fig3.suptitle('Convolucion ' + x + ' con señales periodicas.')

        ax = fig3.add_subplot(2,2,i)
        ax.set_ylabel('impulso')
        ax.set_xlabel('t')
        ax.set_title(x + ' convolucionada con ' + y)
        ax.plot(t,yi)
        i += 1
        energia = integrate.simps(hi**2,t)
        potencia = energia*(1/(b-a))

        print('Energia de entrada de ' + y + ': ', energia)
        print('Potencia de entrada de ' + y + ': ', potencia)

        energia = integrate.simps(yi**2,t)
        potencia = energia*(1/(b-a))

        print('Energia de salida de ' + x + ' convolucionada: ', energia)
        print('Potencia de salida de ' + x + ' convolucionada: ', potencia)
        print('++++++++++++++++++++++++++++++++++++++++++')

    plt.show()    
print('----------------------------------')

#PARTE 4

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

    energia = integrate.simps(h**2,t)
    potencia = energia*(1/(b-a))

    print('Energia de ' + x + ': ', energia)
    print('Potencia de ' + x + ': ', potencia)
    print('+++++++++++++++++++++++++++++++++++')


    plt.show()
print('----------------------------------')

#PARTE 5

i = 1
for x,y in aperiodicas1.items():
    xi = sinc(t)
    hi = y(t)

    yi = np.convolve(xi,hi,'same')*dt

    fig5 = plt.figure("PARTE5")
    fig5.subplots_adjust(hspace=0.5, wspace=0.5)
    fig5.suptitle('Convoluciones entre señales aperiodicas.')

    ax = fig5.add_subplot(2,2,i)
    ax.set_ylabel('impulso')
    ax.set_xlabel('t')
    ax.set_title('Señal SINC convolucionada con ' + x)
    ax.plot(t,yi)
    i += 1

    energia = integrate.simps(hi**2,t)
    potencia = energia*(1/(b-a))

    print('Energia de entrada de ' + x + ': ', energia)
    print('Potencia de entrada de ' + x + ': ', potencia)

    energia = integrate.simps(yi**2,t)
    potencia = energia*(1/(b-a))

    print('Energia de salida de SINC convolucionada: ', energia)
    print('Potencia de salida de SINC convolucionada: ', potencia)
    print('++++++++++++++++++++++++++++++++++++++++++')

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

energia = integrate.simps(hi**2,t)
potencia = energia*(1/(b-a))

print('Energia de entrada de SINC: ', energia)
print('Potencia de entrada de SINC: ', potencia)

energia = integrate.simps(yi**2,t)
potencia = energia*(1/(b-a))

print('Energia de salida de SINC convolucionada: ', energia)
print('Potencia de salida de SINC convolucionada: ', potencia)
print('++++++++++++++++++++++++++++++++++++++++++')

plt.show()