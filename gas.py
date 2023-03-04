# %%
import math

# %%
def ppc(sg):
    '''
    Calculates the psuedocritical pressure (ppc) in psia using the correlation developed by Sutton.
    
    Args:
        sg (float): Specific gravity (dimensionless), range between 0.57 and 1.68.
        
    Returns:
        ppc (float): Pseudocritical pressure in psia.
        
    Raises:
        ValueError: If the input sg is outside the range of 0.57 to 1.68.
    '''
    
    if 0.57 <= sg <= 1.68:
        ppc_value = 756.8 - (131 * sg) - (3.6 * sg**2) 
        return ppc_value   
    else:
        raise ValueError("The ppc correlation only works with specific gravity values between 0.57 and 1.68.")


# %%
def ppr(sg, p):
    '''
    Calculates the pseudoreduced pressure (ppr) as a dimensionless value.
    
    Args:
        sg (float): Specific gravity (dimensionless), range between 0.57 and 1.68.
        p (float): Pressure in psia.
        
    Returns:
        ppr (float): Pseudoreduced pressure (dimensionless).
        
    Raises:
        ValueError: If the input sg is outside the range of 0.57 to 1.68.
        ZeroDivisionError: If the ppc function returns a value of zero.
    '''
    try:
        ppc_value = ppc(sg)
        if ppc_value == 0:
            raise ZeroDivisionError("The ppc function returned a value of zero. Cannot divide by zero.")
        else:
            ppr_value = p / ppc_value
            return ppr_value
    except ValueError as ve:
        raise ve
    except ZeroDivisionError as zde:
        raise zde


# %%
def tpc(sg):
    '''
    Calculates the pseudocritical temperature (tpc) in Rankine scale using the correlation developed by Sutton.
    
    Args:
        sg (float): Specific gravity (dimensionless), range between 0.57 and 1.68.
        
    Returns:
        tpc (float): Pseudocritical tempertarue in rankin.
        
    Raises:
        ValueError: If the input sg is outside the range of 0.57 to 1.68.
    '''

    if 0.57 <= sg <= 1.68:
        tpc_value = 169.2 + (349.5 * sg) - (74 * sg**2) 
        return tpc_value   
    else:
        return ValueError("The tpc correlation only work with specific gravity values between 0.57 and 1.69")
    

# %%
def tpr(sg, t):
    '''
    Calculates the pseudoreduced temperature (tpr) as a dimensionless value.

    Args:
        sg (float): Specifit gravity (dimensionless), range between 0.57 and 1.68.
        t (float): Temperature in Fahrenheit.

    Returns:
        tpr (float): Psudoreduce temperature (dimensionless).

    Raises:
        ValueError: If the input sg is outside the range of 0.57 to 1.68.
        ZeroDivisionError: If the ppc function returns a value of zero.
    '''
    try:    
        tpc_value = tpc(sg)
        if tpc_value == 0:
            raise ZeroDivisionError("The tpc function returned a value of zero. Cannot divide by zero.")
        else:
            tpr_value = ((459.67 + t) / tpc_value)
            return tpr_value
    except ValueError as ve:
        raise ve
    except ZeroDivisionError as zde:
        raise zde



# %%
def secant_method(f, x0, x1, tol=1e-4, max_iter=100):
    """
    Implements the iterative method of secant for finding a root of a function.

    Args:
        f: the function to find the root of.
        x0, x1: the initial guesses for the root.
        tol: the tolerance for the solution.
        max_iter: the maximum number of iterations.

    Returns:
        The approximate root of the function.
    """
    for i in range(max_iter):
        fx0 = f(x0)
        fx1 = f(x1)
        if abs(fx1) < tol:
            return x1
        x = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        if abs(x - x1) < tol:
            return x
        x0, x1 = x1, x
    raise ValueError("The method failed to converge.")


# %%
def z_factor(sg,t,p):
    ''' 
    Calculate the z factor of gas (z_factor), it is used the Dranchuk and Abou-Kassen
    equation of state,more especifically find the root of the equuation of state with 
    the secant method.
    
    Args:
        sg (float): Specific gravity (dimensionless), range between 0.57 and 1.68.
        p (float): Pressure in psia.
        t (float): Temperature in Fahrenheit.

    Returns:
        z_factor (float): Z factor of gas (dimensionless).
     '''
    
    '''
    The firts step calculate the psudoreduced pressure and temperature for 
    the z gas factor.
    '''

    ppr_z= ppr(sg,p)
    tpr_z= tpr(sg,t)

    '''
    The next step is write the equation of state, and is more easy did it by parts.
    '''
    
    def f(z):
        # Constans of the equation of state.
        A1 = 0.3265
        A2 = -1.0700
        A3 = -0.5339
        A4 = 0.01569
        A5 = -0.05165
        A6 = 0.5475
        A7 = -0.7361
        A8 = 0.1844
        A9 = 0.1056
        A10 = 0.6134
        A11 = 0.7210

        # Writting the subequations of equation of state.
        density_ro = 0.27 * (ppr_z / (z * tpr_z))
        c1_tpr = A1 + (A2 / tpr_z) + (A3 / tpr_z**3) + (A4 / tpr_z**4) + (A5 / tpr_z**5)
        c2_tpr = A6 + (A7 / tpr_z) + (A8 / tpr_z**2)
        c3_tpr = A9 * ((A7 / tpr_z) + (A8 / tpr_z**2))
        c4_tpr_ro = A10 * (1 + (A11 * density_ro**2)) * (density_ro**2 / tpr_z**3) * math.exp(-A11 * density_ro**2)

        # Writting the equation of state.
        funtion_fz = z - (1 + (c1_tpr * density_ro) + (c2_tpr * density_ro**2) - (c3_tpr * density_ro**5) + c4_tpr_ro)

        return funtion_fz
    
    # Finding the root with the method of secant, in this case I select a seeds of 0.5 and 1

    root = secant_method(f,0.5,1)

    return root

# %%
def newton_raphson(f, df, x0, tol=1e-4, max_iter=1000):
    """
    Finds a root of the function f using the Newton-Raphson method.

    Parameters:
    f (function): The function to find the root of.
    df (function): The derivative of the function f.
    x0 (float): The initial guess for the root.
    tol (float): The tolerance for the root. Defaults to 1e-4.
    max_iter (int): The maximum number of iterations. Defaults to 100.

    Returns:
    float: The root of the function f.
    """

    x = x0
    for i in range(max_iter):
        fx = f(x)
        if abs(fx) < tol:
            return x
        dfx = df(x)
        if dfx == 0:
            raise ValueError("Newton-Raphson method failed: derivative is zero.")
        x = x - fx / dfx
    raise ValueError("Newton-Raphson method failed: maximum number of iterations exceeded.")


# %%
def z_factor_newton_raphson(sg,t,p):
    ''' 
    Calculate the z factor of gas (z_factor), it is used the Dranchuk and Abou-Kassen
    equation of state, and resolve it with Newton-Raphson method.
    
    Args:
        sg (float): Specific gravity (dimensionless), range between 0.57 and 1.68.
        p (float): Pressure in psia.
        t (float): Temperature in Fahrenheit.

    Returns:
        z_factor (float): Z factor of gas (dimensionless).
     '''
    
    '''
    The firts step calculate the psudoreduced pressure and temperature for 
    the z gas factor.
    '''

    ppr_z= ppr(sg,p)
    tpr_z= tpr(sg,t)

    '''
    The next step is write the equation of state, and is more easy did it by parts.
    '''
    
    # Writing constants of the equation of state

    A1 = 0.3265
    A2 = -1.0700
    A3 = -0.5339
    A4 = 0.01569
    A5 = -0.05165
    A6 = 0.5475
    A7 = -0.7361
    A8 = 0.1844
    A9 = 0.1056
    A10 = 0.6134
    A11 = 0.7210
    
    # Writing a funtion of Dranchuck and Abou-Kassen equiation state

    def f(z):
        
        # Writting the subequations of equation of state.
        density_ro = 0.27 * (ppr_z / (z * tpr_z))
        c1_tpr = A1 + (A2 / tpr_z) + (A3 / tpr_z**3) + (A4 / tpr_z**4) + (A5 / tpr_z**5)
        c2_tpr = A6 + (A7 / tpr_z) + (A8 / tpr_z**2)
        c3_tpr = A9 * ((A7 / tpr_z) + (A8 / tpr_z**2))
        c4_tpr_ro = A10 * (1 + (A11 * density_ro**2)) * (density_ro**2 / tpr_z**3) * math.exp(-A11 * density_ro**2)

        # Writting the equation of state.
        funtion_fz = z - (1 + (c1_tpr * density_ro) + (c2_tpr * density_ro**2) - (c3_tpr * density_ro**5) + c4_tpr_ro)

        return funtion_fz
    
    # Writing a derivative of Dranchuck and Abou-Kassen equiation state
    
    def df(z):
        
        # Writting the subequations of equation of state.
        density_ro = 0.27 * (ppr_z / (z * tpr_z))
        c1_tpr = A1 + (A2 / tpr_z) + (A3 / tpr_z**3) + (A4 / tpr_z**4) + (A5 / tpr_z**5)
        c2_tpr = A6 + (A7 / tpr_z) + (A8 / tpr_z**2)
        c3_tpr = A9 * ((A7 / tpr_z) + (A8 / tpr_z**2))
        
        # Writing the terms of the derivative equation of state.

        term1 = (c1_tpr * density_ro) / z
        term2 = (2 * c2_tpr * density_ro**2) / z
        term3 = (5 * c3_tpr * density_ro**5) / z
        term4 = ((2 * A10 * density_ro**2) / (pow(tpr_z,3) * z)) * ((1 + A11 * density_ro**2) - pow(A11 * density_ro**2, 2) * math.exp(-A11 * density_ro**2))


        #Writting the derivative equation of state
        derivative_fz = 1 + term1 + term2 - term3 + term4

        return derivative_fz

    
    # Finding the root with the method of Newton-Raphson

    root = newton_raphson(f,df,1)

    return root


# %%
def Bg(sg,t,p):
    '''
    Calculate the gas volume factor in (cu ft/SCF), relate the volume of gas in the reseirvoir to
    the volume on the surface (14.7 psia and 60 Farenheit).

    Args:
        sg (float): Specific gravity (dimensionless), range between 0.57 and 1.68.
        p (float): Pressure in psia.
        t (float): Temperature in Fahrenheit.

    Returns:
        Bg (float): Gas volume factor (cu ft/SCF).
    '''

    # Defining the the variables 
    pressure_reseirvoir = p
    temperature_reseirvoir = 459.67 + t 
    z_reseirvoir_conditions = z_factor_newton_raphson(sg,t,p)
    
    # Calculating the gas volume factor

    Bg = 0.02829 * ((z_reseirvoir_conditions * temperature_reseirvoir) / pressure_reseirvoir)

    return Bg


# %%
def Bg_bbl(sg,t,p):
    '''
    Calculate the gas volume factor in (bbl/SCF), relate the volume of gas in the reseirvoir to
    the volume on the surface (14.7 psia and 60 Farenheit).

    Args:
        sg (float): Specific gravity (dimensionless), range between 0.57 and 1.68.
        p (float): Pressure in psia.
        t (float): Temperature in Fahrenheit.

    Returns:
        Bg__bbl (float): Gas volume factor (bbl/SCF).
    '''

    # Defining the the variables 
    pressure_reseirvoir = p
    temperature_reseirvoir = 459.67 + t 
    z_reseirvoir_conditions = z_factor_newton_raphson(sg,t,p)
    
    # Calculating the gas volume factor

    Bg = 0.00504 * ((z_reseirvoir_conditions * temperature_reseirvoir) / pressure_reseirvoir)

    return Bg

# %%
def reseirvoir_gas_density(sg,t,p):
    '''
    Calculate the density of a reseirvoir gas.

    Args:
        sg (float): Specific gravity (dimensionless), range between 0.57 and 1.68.
        p (float): Pressure in psia.
        t (float): Temperature in Fahrenheit.

    Returns:
        reseirvoir_gas_density (float): Density of a reseirvoir gas (lb/cu ft)
    '''

    density = (28.97 * sg * p) / (z_factor(sg,t,p) * 10.73 * (t + 459.67))

    return density

# %%
def gas_viscosity(sg,t,p):
    '''
    Calculate the gas viscocity with the correlation develop by Lee, Gonzales and Eakin.

    Args:
        sg (float): Specific gravity (dimensionless), range between 0.57 and 1.68.
        p (float): Pressure in psia.
        t (float): Temperature in Fahrenheit.

    Returns:
        gas_viscosity (float): Viscosity in cp.
    '''
    
    # Calculating the molecular weight.
    mw = 28.97 * sg
    
    # Writing the subequations of correlation.
    density = (0.0014935 * p * mw) / (z_factor_newton_raphson(sg,t,p) * (459.67 + t))
    k = ((9.4 + (0.02 * mw)) * ((459.67 + t)**1.5)) / (209 + (19 * mw) + (459.67 + t))
    x = 3.5 + (986 / (459.67 + t)) + (0.01 * mw)
    y = 2.4 - (0.2*x)

    mg = 0.0001 * k * math.exp(x * density**y)

    return mg



