def advection(D):
    ########## DESCRIPTION ###########################################
    # This function implements the advection equation to move a single particle through a gas flow field.
    # The only parameter is the particle diameter (D).
    # These are the outputs:
    # t_arr  - array of time
    # up_arr - array of particle velocity
    # ug_arr - array of gas velocity
    # Fd_arr - array of particle drag force
    # xp_arr - array of particle position
    ##################################################################
    
    ##### INITIALIZATION #################################
    c     = 10 # m/s                          # initial gas velocity at spike
    e     = 1 # m/s                           # initial gas velocity elsewhere
    N     = 100                               # number of plot points   
    Niter = 100                               # number of iterations for time
    x     = np.linspace(0,100,N)              #
    dx    = x[1]-x[0]                         # change in x
    dt    = 0.5*dx/c                          # change in t
    rho_g = 1.6 # kg/m^3                      # gas density
    dia_p = D                                 # particle diameter
    A     = np.pi*0.25*dia_p**2               # particle projected area
    rho_p = 2500 # kg/m^3                     # particle density
    M     = rho_p * 1/6 * dia_p**3 * np.pi    # particle mass
    mu    = 2e-5 # kg/ms                      # dynamic viscosity

    # ARRAYS
    Fd_arr = []
    Fd_arr.append(0.0)
    up_arr = []
    up_arr.append(0.0)
    xp_arr = []
    xp_arr.append(50.0)
    ug_arr = []
    ug_arr.append(e)
    t_arr = []
    t_arr.append(0)
    ######################################################

    u = np.ones(N)
    u[0] = e
    u[10:20] = 5

    for j in range(Niter):
        t_arr.append(t_arr[-1] + dt)

        for i in range(N-1):
            uL = u[i]
            uR = u[i+1]
            u[i+1] = (uL*c*dt + uR*(dx - c*dt))/dx

        u[0] = u[N-1]

        ############# PARTICLE SOLVER #############################
        ug_temp = kernel2(xp_arr[-1] , 'p' , u , dx)
        ug_arr.append(ug_temp)

        Re = rho_g * abs(ug_arr[-1] - up_arr[-1]) * dia_p / mu # Reynolds number

        # Cd = Cd_Sternin(0) 
        Cd = Cd_Loth(Re, ug_temp, up_arr[-1], 300, 300, dia_p, rho_g) # Tg = 300 K, Tp = 300 K
        Fd_arr.append(Cd * rho_g * A * abs(ug_temp - up_arr[-1]) * (ug_temp - up_arr[-1]) / 2)

        up_arr.append(up_arr[-1] + Fd_arr[-1] * dt / M)

        xp_arr.append(xp_arr[-1] + up_arr[-1] * dt)
        ###########################################################
    
    return t_arr, up_arr, ug_arr, Fd_arr, xp_arr