from sympy import init_session,plot,Symbol,diff,cos,exp,sympify,N,lambdify
import math
from numpy import sign
import numpy as np
from numpy.lib.scimath import sqrt
from matplotlib import pyplot as plt
#Definicion de variable
x = Symbol('x')
#Definicion de funciones
f1= (5/2)*exp(3*(x-1)) - (15/4)*exp(2*(x-1)) - (1/4)*exp(x-1) + ((1/5)*x)
f2= 0.7 - 0.3*cos(x) + 0.25*x
f21= f2 - 1
#Derivar f1
df1 = f1.diff(x)
#Transformar a expresion numerica la expresion simbolica df1 y f2
f = lambdify(x,df1)
fn1= lambdify(x,f1)
fn2 = lambdify(x,f2)
fn21 = lambdify(x,f21)
#Metodo de Muller
def Muller(f,x0,x1,x2,tol):
    error = 1e3
    x3 = 0
    it = 1
    while error > tol:
        c = f(x2)
        b = ((x0-x2)**2*(f(x1)-f(x2))-(x1-x2)**2*(f(x0)-f(x2)))/((x0-x2)*(x1-x2)*(x0-x1))
        a = ((x1-x2)*(f(x0)-f(x2))-(x0-x2)*(f(x1)-f(x2)))/((x0-x2)*(x1-x2)*(x0-x1))
        x3 = x2-(2*c)/(b+sign(b)*sqrt(b**2-4*a*c))
        error = abs((x3-x2)/x3)
        x0 = x1; x1 = x2; x2 = x3
        it+=1
    return [it,x3,error]
#Metodo de Newton
def Newton(f,df,xi,tol):
    x = xi
    error = 1e3
    n = 1
    while error > tol:
        x = x-f(x)/df(x)
        error = abs(f(x))
        n += 1
    return[n,x,error]
###Columna Izquierda
##Valuacion del metodo de Muller
iteraciones,raiz,error = Muller(f,-1.25,-1.5,-2,(1e-4))
#raiz,error = Newton(f,lambdify(x,df1.diff(x)),-1,1e-4)
valorEnRaiz = fn2(raiz)
#Intervalo de la raiz teniendo en cuenta el error
a = raiz-error
b = raiz+error
#Intervalo de la imagen
c = fn2(a)
d = fn2(b)
#Region bidimensional donde se encuentra el extremo de la columna izquierda
region = [(a,c),(a,d),(b,c),(b,d)]
#Area de la region anterior
area = (b-a)*(d-c)
print("Con el metodo de Muller")
print("Los puntos que conforman la region bidimensional donde se encuentra el extremo de la columna izquierda:",region)
print("Con el area: ",area)
print("Con la cantidad de iteraciones: ",iteraciones)
###Valuacion del metodo de Newton
iteracionesN,raizN,errorN = Newton(f,lambdify(x,df1.diff(x)),-1.25,1e-4)
valorEnRaiz = fn2(raiz)
#Intervalo de la raiz teniendo en cuenta el error
aN = raizN-errorN
bN = raizN+errorN
#Intervalo de la imagen
cN = fn2(aN)
dN = fn2(bN)
#Region bidimensional donde se encuentra el extremo de la columna izquierda
regionN = [(aN,cN),(aN,dN),(bN,cN),(bN,dN)]
#Area de la region anterior
areaN = (bN-aN)*(dN-cN)
print("Con el metodo de Newton")
print("Los puntos que conforman la region bidimensional donde se encuentra el extremo de la columna izquierda:",regionN)
print("Con el area: ",areaN)
print("Con la cantidad de iteraciones: ",iteracionesN)

##Impresion del resultado
##intervalo
#inter= np.arange(-2,1.7,step=0.01)
##Funciones
#plt.plot(inter,fn1(inter))
#plt.plot(inter,fn2(inter))
##Puntos que conforman la region
#plt.plot(a,c, marker="o", color="black")
#plt.plot(a,d, marker="o", color="blue")
#plt.plot(b,c, marker="o", color="green")
#plt.plot(b,d, marker="o", color="red")
#plt.plot(raiz,valorEnRaiz, marker="o", color="yellow")
#plt.grid()
#plt.show()

###Para obtener posicion de la columna derecha:
iteracionesf2,raizf2,errorf2 = Muller(fn21,-1,-0.5,0,(1e-4))
#Distancia entre columnas
distanciaColumnas = sqrt((raiz-raizf2)**2+(fn2(raiz)-fn2(raizf2))**2)
print(distanciaColumnas)
##Impresion del resultado
inter= np.arange(-2,1.7,step=0.01)
plt.plot(inter,fn1(inter))
plt.plot(inter,fn2(inter))
plt.plot(raiz,fn2(raiz),marker="o", color="black")
plt.plot(raizf2,fn2(raizf2),marker="o", color="black")
plt.grid()
plt.show()




