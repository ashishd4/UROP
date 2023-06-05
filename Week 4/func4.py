def Cd_Sternin(Re): # Reynolds number
    ########## DESCRIPTION ###########################################
    # This function derives and outputs the drag coefficient of a particle.
    # The only parameter is the Reynolds number (Re).
    ##################################################################
    
    if Re == 0: 
        return 0.42
    else:
        return 24/Re + 4.4/np.sqrt(Re) + 0.42

def Cd_Loth(Re, ug, up, Tg, Tp, D, rhog):
    ########## DESCRIPTION ###########################################
    # This function derives and outputs the drag coefficient of a particle.
    # These are the parameters: 
    # Re   - Reynolds number
    # ug   - Gas velocity, in m/s
    # up   - Particle velocity, in m/s
    # Tg   - Gas temperature, in K
    # Tp   - Particle temperature, in K
    # D    - Particle diameter, in m
    # rhog - Gas density, in kg/m^3
    ##################################################################
    
    import math
    import numpy as np
    
    R = 287 # J/(kg*K)
    
    gamma = 1.4
    Ma = abs(ug - up)/np.sqrt(gamma * R * Tg) # Mach number
    if Re == 0:
        return 0
    
    if Re > 45:
        if Ma > 1.5:
            CM = 2.18 - 0.13*math.tanh(0.9*Ma - 2.7)
        else:
            CM = 1.65 + 0.65*math.tanh(4*Ma - 3.4)
        
        if Ma > 0.8:
            GM = 5 + 40/(Ma**3)
        else:
            GM = 166*Ma**3 + 3.29*Ma**2 - 10.9*Ma + 20
        
        if Ma > 1:
            HM = 0.93 + 1/(3.5 + Ma**5)
        else:
            HM = 0.0239*Ma**3 + 0.212*Ma**2 - 0.074*Ma + 1
    
        return 24*(1 + 0.15*Re**0.687)*HM/Re + (0.42*CM)/(1 + 42500*GM*Re**(-1.16))

    if Ma == 0: # speed ratio
        s = 1
    else:
        s = Ma*np.sqrt(gamma/2)
        
    mg = 4.65e-26 # kg, atomic mass of gas
    dg = 4.16e-10 # m, molecule diameter
    lamb = 1/np.sqrt(2) * 1/(np.pi*dg**2) * mg/rhog
    Kn = lamb/D

    C_apdfm = (1 + 2*s**2)*math.exp(-s**2)/(np.sqrt(np.pi)*s**3) + (4*s**4 + 4*s**2 - 1)*math.erf(s)/(2*s**4)
    C_dfm = C_apdfm + 2*np.sqrt(np.pi*Tp/Tg)/(3*s)
    C_dfmRe = C_dfm/(1 + np.sqrt(Re/45)*(C_apdfm/1.63 - 1))
    
    f_Kn = 1/(1 + Kn*(2.514 + 0.8*math.exp(-0.55/Kn)))
    
    C_dKnRe = 24/Re * (1 + 0.15*Re**0.687)*f_Kn
        
    return (C_dKnRe + Ma**4*C_dfmRe)/(1 + Ma**4)

def kernel(S, arg):
    ########## DESCRIPTION ###########################################
    # This kernel function derives and outputs the gas velocity distribution in a small region.
    # These are the parameters: 
    # S   - velocity of gas at a specific point
    # arg - argument for what to return
    #       c - return S
    #       p - return a numpy array that distributes the value of S as a 25-50-25 distribution
    #       anything else - return 0
    ##################################################################
    
    import numpy as np
    
    if (arg == 'c'):
        return S
    if (arg == 'p'):
        return S*np.array([0.25, 0.5, 0.25])
    return 0

def kernel2(xp, arg, u, dx):
    ########## DESCRIPTION ###########################################
    # This kernel function derives and outputs the gas velocity distribution in a small region.
    # These are the parameters: 
    # xp  - particle position
    # arg - argument for what to return
    #       c - return S
    #       p - return the sum of gas velocities around the particle as implemented below
    #       anything else - return 0
    # u   - array of gas velocities
    # dx  - change in distance
    ##################################################################
    
    i = int(xp/dx)
    
    if (arg == 'c'):
        return u[i]
    if (arg == 'p'):
        return 0.25*u[i-1] + 0.5*u[i] + 0.25*u[i+1]
    return 0