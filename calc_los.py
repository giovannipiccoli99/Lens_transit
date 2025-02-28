import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

def calc_r(s, b):
    '''
    Calculates the distance to the center of the planet.
    Input values:
    s: path variable (-inf initially, 0 when passing the planet, infinity at the end)
    b: impact parameter
    '''
    return np.sqrt(s**2 + b**2)

def nu(s, b):
    '''
    Input from Thomas here
    '''
    r = calc_r(s, b)

    #dummy function
    return 10**(-3)

def calc_los(b):
    '''
    Calculates the line-of-sight integral
    Input values:
    b: Impact parameter
    '''
    return integrate.quad(lambda s: nu(s, b), -np.inf, np.inf)

def main():

    b = 1
    l = calc_los(b)
    print(l)

if __name__ == "__main__":
    main()

