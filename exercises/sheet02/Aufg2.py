import numpy as np 

'''
- Das Integral muss definiert werden
    - k1 muss definiert werden
- y_n ist Startwert
- y_n1 ergibt sich aus Startwert += Integral
'''

#DGL-Integrationsschritt
def rungekutta(f, S, a, b, h, y,
               p_La_gamma, rho_ga2_gamma):           #Funktion f, untere/obere Grenten, Schrittweite h, Wert y
    
    # Schritte für k1, etc. müssen parallel ausgeführt werden
    if f is dPsi_dS:
        R, Z, Psi = y[0], y[1], y[2]
        
        k1 = f(*y, p_La_gamma, rho_ga2_gamma)
        k2 = f(R, Z, Psi + k1*h/2, p_La_gamma, rho_ga2_gamma)
        k3 = f(R, Z, Psi + k2*h/2, p_La_gamma, rho_ga2_gamma)
        k4 = f(R, Z, Psi + k3*h, p_La_gamma, rho_ga2_gamma)
    else:
        k1 = f(y)
        k2 = f(y + k1*h/2)
        k3 = f(y + k2*h/2)
        k4 = f(y + h*k3)
        
    
    integral = (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    
    
    return integral

# Defintion der DGLn

def dR_dS(Psi):
    return np.cos(Psi)

def dZ_dS(Psi):
    return np.sin(Psi)

def dPsi_dS(R, Z, Psi, p_La_gamma, rho_ga2_gamma):
    if (R == 0):
        dPsi_dS = (p_La_gamma - rho_ga2_gamma * Z) / 2
    else:
        dPsi_dS = p_La_gamma - rho_ga2_gamma * Z - np.sin(Psi) / R
    return dPsi_dS


def rk4(R_0, S_0, Z_0, Psi_0,
        s_max, h, n,
        p_La_gamma, rho_ga2_gamma):
    
    R_Z_Psi_S = [(R_0, Z_0, Psi_0, S_0)]
    
    '''
    Form von  R_Z_Psi   
    i    R   Z   Psi  S
    0    0   0   0    0
    1    1   2   3    0.1
    .
    .
    .
    '''
    
    # R, Z, Psi, S = R_0, Z_0, Psi_0, S_0
    # values = [R,Z,Psi]
    
    # for i in range(len(n)):
    #     for v in values:
    #         values[i] = rungekutta(dR_dS, n[i], n[i+1], h, v)   # Funktioniert noch nicht, da values[i] nur von i=0-2 iteriert werden soll
    
    R, Z, Psi, S = [R_0], [Z_0], [Psi_0], [S_0]
    
    for i in range(len(n) - 1):
        R.append(rungekutta(    dR_dS,    S[i], n[i], n[i+1], h,  Psi[i],                  p_La_gamma, rho_ga2_gamma))
        Z.append(rungekutta(    dZ_dS,    S[i], n[i], n[i+1], h,  Psi[i],                  p_La_gamma, rho_ga2_gamma))
        Psi.append(rungekutta(  dPsi_dS,  S[i], n[i], n[i+1], h,  (R[i], Z[i], Psi[i]),    p_La_gamma, rho_ga2_gamma))
        
        S.append(S[i] + h)
        
        R_Z_Psi_S.append( (R[i+1], Z[i+1], Psi[i+1], S[i+1]) )
      
    return R_Z_Psi_S


def main():
    # Konstanten
    # Eventuell variieren
    s_max = 10      # obere Grenze als Anfangswert, hier variieren durch probieren
    h = 0.1         # Schrittweite
    n = np.arange(0, s_max, h)
    
    # Durch Aufgabe gegeben
    R_0, S_0, Z_0, Psi_0 = 0, 0, 0, 0       # Anfangswerte
    p_La_gamma, rho_ga2_gamma = 2, 0.1
    
    results = rk4(R_0, S_0, Z_0, Psi_0,
              s_max, h, n,
              p_La_gamma, rho_ga2_gamma)
    
    for i in range(len(results)):
        print(results[i])
    
main()