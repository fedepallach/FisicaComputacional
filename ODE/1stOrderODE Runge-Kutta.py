import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def graficar(x_values, y_values):
    plt.plot(x_values, y_values,label='Euler modificado')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("Solving ODE using Euler's Method")
    plt.legend()
    plt.grid(True)
    plt.show()

def euler_simp_p(f, xi, yi):
    p = f(xi, yi)
    return p

def euler_modif_p(f, x, y, h):
    ymit = y + euler_simp_p(f,x,y) * h/2
    pendienteMitad = f(x+(h/2), ymit)
    return pendienteMitad

def euler_modif_p_k2(f, x, y, h):
    ymit2 = y + euler_modif_p(f,x,y,h) * h/2  # igual es: y + k1  * h/2
    pendienteMitad2 = f(x+(h/2), ymit2)
    return pendienteMitad2

def para_k3(f,x,y,h):
    pend = f(x + h/2 , y + euler_modif_p_k2(f,x,y,h) * h)
    return pend

def euler_mej_p(f,x,y,h):  #NO se usa para RUNGE-KUTTA
    pend1 = f(x,y)
    ysig = y + pend1 * h
    pend2 = f(x+h,ysig)
    pendProm = (pend1+pend2)/2
    return pendProm

def runge_kutta(f,x,y,h):
    k0 = euler_simp_p(f,x,y)
    k1 = euler_modif_p(f,x,y,h)
    k2 = euler_modif_p_k2(f,x,y,h)
    k3 = para_k3(f,x,y,h)
    resultado = y + (h/6) * (k0 + 2*k1 + 2*k2 + k3)
    return resultado

def obtenerValores(f,x0,y0,h,n):
    x_values = np.zeros(n+1)
    y_values = np.zeros(n+1)
    x_values[0]=x0
    y_values[0]=y0
    
    for i in range(1,n+1):
        x_values[i] = x_values[i-1] + h 
        y_values[i] = runge_kutta(f,x_values[i-1],y_values[i-1],h)
    
    return x_values,y_values

def main ():
    x0 = 0
    y0 = 10
    k = .2   #func de superficie de contacto, litros en la alberca, prop del material, interac igual si viento otra cosa dif no?
    h   = .025
    n = 800

    #Funcion
    def f(x, y):
        return -k * y

    # Resolver
    x_values, y_values = obtenerValores(f, x0, y0, h, n)

    #  Plot the solution
    #label = 'ODE Euler Modificado'
    graficar(x_values, y_values)

    graficarElError(x_values, y_values,n)




#Para ver el:  ERROR

def valorDeFuncion(x,k):
    valor = 10*math.exp(-k*x)
    return valor

def graficarElError(x_values, y_values, n):
    #Graficar y Analisis del error
    meanAbsoluteErrorList = []   #solo correrlo al inicio
    k =.2
    listRealValues = []
    listErrors = []
    for i in range (0,n+1):
        #print(x_values[i])
        #print(y_values[i])
        valorReal = valorDeFuncion(x_values[i],k)
        listRealValues.append(valorReal)
        listErrors.append(valorReal-y_values[i])
    #print(listRealValues)
    arrayRealValues = np.array(listRealValues)
    #print(arrayRealValues)
    arrayErrorValues = np.array(listErrors)
    #print(arrayErrorValues)

    # Plot the solution
    plt.plot(x_values, arrayErrorValues, label='Euler Mod Method Error')
    plt.xlabel('x')
    plt.ylabel('Error')
    plt.title("Solving ODE using Euler's Method")
    plt.legend()
    plt.grid(True)
    plt.show()

    #Mean Absolute Error
    meanAbsoluteError = np.mean(np.abs(arrayErrorValues))
    meanAbsoluteErrorList.append(meanAbsoluteError)

    print(meanAbsoluteErrorList)
    graficar(x_values, y_values)
    graficar(x_values, arrayRealValues)



#Lo que corre cuando en terminal: python nombreArch.py
main()