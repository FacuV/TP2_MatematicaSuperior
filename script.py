from sympy import init_session,plot,Symbol,diff,cos,exp,sympify,N,lambdify
import math
from numpy import sign
import numpy as np
from numpy.lib.scimath import sqrt
from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
#Definicion de variable
x = Symbol('x')
#Definicion de funciones
f1= (5/2)*exp(3*(x-1)) - (15/4)*exp(2*(x-1)) - (1/4)*exp(x-1) + ((1/5)*x)
f2= 0.7 - 0.3*cos(x) + 0.25*x
#Derivar f1
df1 = f1.diff(x)
#Transformar a expresion numerica la expresion simbolica df1 y f2
f = lambdify(x,df1)
fn1= lambdify(x,f1)
fn2 = lambdify(x,f2)
#Metodo de Muller
def Muller(f,x0,x1,x2,tol):
    error = 1e3
    x3 = 0
    while error > tol:
        c = f(x2)
        b = ((x0-x2)**2*(f(x1)-f(x2))-(x1-x2)**2*(f(x0)-f(x2)))/((x0-x2)*(x1-x2)*(x0-x1))
        a = ((x1-x2)*(f(x0)-f(x2))-(x0-x2)*(f(x1)-f(x2)))/((x0-x2)*(x1-x2)*(x0-x1))
        x3 = x2-(2*c)/(b+sign(b)*sqrt(b**2-4*a*c))
        error = abs((x3-x2)/x3)
        x0 = x1; x1 = x2; x2 = x3
    return [x3,error]

#Valuacion del metodo de Muller
raiz,error = Muller(f,-1.25,-1.5,-2,1e-4)
#Intervalo de la raiz teniendo en cuenta el error
a = raiz-error
b = raiz+error
#Intervalo de la imagen
c = fn2(a)
d = fn2(b)
print(a,b)
print(c,d)
#Region bidimensional donde se encuentra el extremo de la columna izquierda
region = [(a,c),(a,d),(b,c),(b,d)]
#Area de la region anterior
area = (b-a)*(d-c)
print(region,area)
#Impresion del resultado
#intervalo
inter= np.arange(-2,1.7,step=0.01)
#Funciones
plt.plot(inter,fn1(inter))
plt.plot(inter,fn2(inter))
#Puntos que conforman la region
plt.plot(a,c, marker="o", color="black")
plt.plot(a,d, marker="o", color="black")
plt.plot(b,c, marker="o", color="black")
plt.plot(b,d, marker="o", color="black")
plt.grid()
plt.show()
