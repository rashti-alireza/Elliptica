def print_SLy_pwp(rho0_th, K0, Gammas, name, rho0_i, rho0_f, points):
    """Generates (.txt) data table for piecewise-polytropic EOS 
       in format [pressure] [rest-mass density] [energy density] [enthalpy].
       Takes input:
            (int) n: number of pwp intervals, 
            (float) K0: initial polytropic coefficient,
            (float list) rho0_th: pwp intervals in terms of rest-mass density,
            (float list) Gammas: polytrope gamma for each interval,
            (string) name: name of output file.
            (float) rho0_i and rho0_f: specifying range of values to print,
            (int) points: number of points to print.
        Assumes rho0 is increasing. Outputs in uniform rho0 steps."""
            
    with open(name, "w") as output_file:
        rho0 = rho0_i
        drho0 = (rho0_f - rho0_i)/points
        Ks = generate_Ks(rho0_th, K0, Gammas)
        
        for pt in range(points):
            data_pt = SLy_pwp(rho0_th, Ks, Gammas, rho0) 
            output_file.write(str(data_pt[0]) + " " + str(data_pt[1]) + " " + str(data_pt[2]) + " " + str(data_pt[3]) + "\n")
            rho0 += drho0

    
            
def SLy_pwp(rho0_th, Ks, Gammas, rho0_val):
    """Takes input:
            (float list) rho0_th: pwp intervals in terms of rest-mass density,
            (float list) Ks: polytropic coefficient for each interval,
            (float list) Gammas: polytrope gamma for each interval,
            (float) rho0_val: rest-mass density to evaluate.
        Returns (float list) [pressure, rest-mass density, energy density, enthalpy]
        calculated from piecewise polytrope."""
    
    output = [pressure(rho0_th, Ks, Gammas, rho0_val),
              rho0_val,
              energy_density(rho0_th, Ks, Gammas, rho0_val),
              enthalpy(rho0_th, Ks, Gammas, rho0_val)]
    
    return output
    
def generate_Ks(rho0_th, K0, Gammas):
    """Takes input:
            (float list) rho0_th: pwp intervals in terms of rest-mass density,
            float) K0: initial polytropic coefficient,
            (float list) Gammas: polytrope gamma for each interval.
        Returns (float list) of polytropic coefficients for piecewise polytrope.
        Reference: eq 11. https://arxiv.org/pdf/0812.2163.pdf."""
        
    Ks = [K0]
    for j in range(1,len(Gammas)):
        K_next = pressure(rho0_th, Ks, Gammas, rho0_th[j-1]) / (rho0_th[j-1] ** Gammas[j])
        Ks.append(K_next)
            
    return Ks
        
def find_interval(rho0_val, rho0_th):
    """Takes input:
            (float) rho0_val: rest-mass density value,
            (float list) rho0_th: pwp intervals in terms of rest-mass density.
       Returns (int) the interval of rho0_th that rho0_val occupies."""
    
    if (rho0_val < rho0_th[0]): return 0
    elif (rho0_val > rho0_th[-1]): return len(rho0_th)
    else:
        ctr = 0
        while (rho0_val > rho0_th[ctr]):
            ctr += 1
        return ctr
        
def pressure(rho0_th, Ks, Gammas, rho0_val):
    interval = find_interval(rho0_val, rho0_th)
    
    return Ks[interval] * (rho0_val ** Gammas[interval])

def energy_density(rho0_th, Ks, Gammas, rho0_val):
    interval = find_interval(rho0_val, rho0_th)
    
    return rho0_val + pressure(rho0_th, Ks, Gammas, rho0_val) / (Gammas[interval] - 1)
    
def enthalpy(rho0_th, Ks, Gammas, rho0_val):
    return (pressure(rho0_th, Ks, Gammas, rho0_val) + energy_density(rho0_th, Ks, Gammas, rho0_val)) / rho0_val
    

if __name__ == "__main__":
    points = 1000
    rho0_th = [2.3674e-04, 8.1147e-04, 1.6191e-03]
    K0 = 8.9493e-02
    Gammas = [1.3569e+00, 3.0050e+00, 2.9880e+00, 2.8510e+00]
    name = "eos_test_table.txt"
    rho0_i = 1E-9
    rho0_f = 3e-3
    Ks = generate_Ks(rho0_th, K0, Gammas)
    
    print_SLy_pwp(rho0_th, K0, Gammas, name, rho0_i, rho0_f, points)




        
