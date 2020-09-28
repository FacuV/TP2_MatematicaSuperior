import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
from numpy import exp
#Funciones F1 y F2:
def f1(x): return (5/2)*exp(3*(x-1)) - (15/4)*exp(2*(x-1)) - (1/4)*exp(x-1) + ((1/5)*x)
def df1(x): return sy.Derivative(f1(x),x)

def f2(x): return 0.7 - 0.3*np.cos(x) + 0.25*x
#Intervalo util para graficar:
x = np.arange(-2,1.7,0.01)
#x = np.linspace(-2,1.5,100)
#print(df1(x).getn)
plt.plot(x,f1(x),label="f1(x)")
plt.plot(x,df1(x).getn,label="f2(x)")
plt.legend()
plt.grid()
plt.xlabel("x")
plt.ylabel("y")
plt.title("f(x)")
plt.show() 


