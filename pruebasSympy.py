from sympy import init_session,plot,Symbol,diff,cos,exp,sympify,N,lambdify
import math
init_session(use_latex=True)
x = Symbol('x')
#f1 = Function('f1') ((5/2)*exp(3*(x-1)) - (15/4)*exp(2*(x-1)) - (1/4)*exp(x-1) + ((1/5)*x))
#f2 = Function('f2') (0.7 - 0.3*cos(x) + 0.25*x)
#df1 = Function('df1') (f1.diff(x))
#aca lo estoy cambiando de funcion a Add
f1= (5/2)*exp(3*(x-1)) - (15/4)*exp(2*(x-1)) - (1/4)*exp(x-1) + ((1/5)*x)
f2= 0.7 - 0.3*cos(x) + 0.25*x
df1 = (f1.diff(x))
ddf1 = (df1.diff(x))
#evaluar
#print(sympify(df1).subs(x,1))

#Metodo de Newton
x0 = -1
i = 0
while(i < 0):
    x1 = N(x0 - sympify(df1).subs(x,x0)/sympify(ddf1).subs(x,x0),n=100)
    print(x1)
    x0=x1
    i+=1
#
f = lambdify(x,f1)
#Metodo de Muller
def Muller(a, b, c): 
    MAX_ITERATIONS = 10000
    res = 0; 
    i = 0; 
    while (True): 
        f1 = f(a); f2 = f(b); f3 = f(c); 
        d1 = f1 - f3;  
        d2 = f2 - f3; 
        h1 = a - c;  
        h2 = b - c; 
        a0 = f3; 
        a1 = (((d2 * pow(h1, 2)) - 
               (d1 * pow(h2, 2))) / 
              ((h1 * h2) * (h1 - h2))); 
        a2 = (((d1 * h2) - (d2 * h1)) / 
              ((h1 * h2) * (h1 - h2))); 
        x = ((-2 * a0) / (a1 + 
             abs(math.sqrt(a1 * a1 - 4 * a0 * a2)))); 
        y = ((-2 * a0) / (a1 - 
            abs(math.sqrt(a1 * a1 - 4 * a0 * a2)))); 

        # Taking the root which is  
        # closer to x2 
        if (x >= y): 
            res = x + c; 
        else: 
            res = y + c; 
  
        # checking for resemblance of x3  
        # with x2 till two decimal places 
        m = res * 100; 
        n = c * 100; 
        m = math.floor(m); 
        n = math.floor(n); 
        if (m == n): 
            break; 
        a = b; 
        b = c; 
        c = res; 
        if (i > MAX_ITERATIONS): 
            return null
            break; 
        i += 1; 
    if (i <= MAX_ITERATIONS): 
        return round(res, 10)

a = 0; 
b = 1; 
c = 2; 
Muller(a,b,c)
#plot(f1,(x,-2,1.5))
#plot(df1,(x,-2,1.7))
#plot(f2,(x,-3,3))