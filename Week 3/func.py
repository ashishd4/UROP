def Cd_Sternin(Re): # Reynolds number
    if Re == 0: 
        return 0.42
    else:
        return 24/Re + 4.4/np.sqrt(Re) + 0.42

def kernel(S, arg):
    import numpy as np
    
    if (arg == 'c'):
        return S
    if (arg == 'p'):
        return S*np.array([0.25, 0.5, 0.25])
    return 0

def kernel2(xp, arg, u, dx):
    i = int(xp/dx)
    
    if (arg == 'c'):
        return u[i]
    if (arg == 'p'):
        return 0.25*u[i-1] + 0.5*u[i] + 0.25*u[i+1]
    return 0