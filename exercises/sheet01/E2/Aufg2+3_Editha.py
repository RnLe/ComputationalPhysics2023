import numpy as np


def trapezoid(int_f, lower_a, upper_b, width_h):    #Typ-Zuweisung: aufpassen
    num_n = int( (upper_b - lower_a) / width_h )    #als int definieren, da sonst Datentyp nicht klar
    x_k = np.zeros(num_n+1)                         #bildet numpy-array der Länge num_n, vollst. mit nullen (zeros) gefüllt
    for k in range(num_n+1):       
        x_k[k] = lower_a + k * width_h            
    y = np.zeros(num_n+1)  
    #y = np.array([int_f(x) for x in x_k])          #kurze Version der nachfolgenden for-schleife                              
    for i in range(len(x_k)):       
        y[i] = int_f(x_k[i]) 
    integral = width_h * ( 1/2 * ( y[0]  ) + np.sum(y[1:num_n-1]) + 1/2 * y[num_n])
    return integral


def mittelpunktsregel(f, a, b, h):
    n = int( (b - a ) / h )
    integral = 0
    for k in range(n):
        x_k = a + (k+1/2)*h
        integral += h * f(x_k)
    return integral

def simpsonregel(f, a, b, h):
    n = int( (b - a ) / h )
    integral = (h/3) * (f(a) + f(b)) 
    for k in range(n):
        if (k%2) == 0:               #gerade:
            x_k = a + k*h 
            multiplicator = 2/3  
            integral += h * (2/3) * f(x_k)
        else:
           x_k = a + k*h
           multiplicator = 4/3
           integral += h * (4/3) * f(x_k) 
    return integral
             

def expontential(x):
    return np.exp(-x)/x

def sinus(x):                   #numerisches Problem; explizite Grenzwertbildung
    if x == 0:                  # evtl bereich von Betrag x einschränken für kleine x
        return 0
    else:
        return x * np.sin(1/x)

def test_function(x,a):
    return x*a
 

def integral_berechnen(integration_rule, integrand, lower_a, upper_b):
    altes_ergebnis = 0
    neues_ergebnis = 1
    rel_abweichung = (altes_ergebnis - neues_ergebnis) / neues_ergebnis
    h = 1
    while np.abs(rel_abweichung) > 10**(-4):
        h = h/2
        altes_ergebnis = neues_ergebnis
        neues_ergebnis = integration_rule(integrand, lower_a, upper_b, h)
        rel_abweichung = (altes_ergebnis - neues_ergebnis) / neues_ergebnis
    else:
        print(rel_abweichung)
        print('h = ', h )
        print('Ergebnis: ', neues_ergebnis)

#integral_berechnen( trapezoid, expontential, 1, 100 )

def main():
   integration_rule=[trapezoid, mittelpunktsregel, simpsonregel]
   integrand=[expontential, sinus]
   for i in integration_rule:
       for j in integrand:
           if integrand == expontential:
               a = 1
               b = 100
           else:
               a = 0
               b = 1
           integral_berechnen(i, j, a, b )

# integral_berechnen( trapezoid, sinus, 0, 1 )
# integral_berechnen( mittelpunktsregel, sinus, 0, 1 )
# integral_berechnen( simpsonregel, sinus, 0, 1 )
# integral_berechnen( trapezoid, expontential, 1, 100 )
# integral_berechnen( mittelpunktsregel, expontential, 1, 100 )
# integral_berechnen( simpsonregel, expontential, 1, 100 )

main()