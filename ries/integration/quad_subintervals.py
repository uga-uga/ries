import numpy as np
from scipy.integrate import quad

def quad_subintervals(f, x):
    integral = 0.
    uncertainty_estimate = 0.

    for i in range(len(x)-1):
        quad_result = quad(f, x[i], x[i+1])
        integral += quad_result[0]
        uncertainty_estimate += quad_result[1]*quad_result[1]

    return (integral, np.sqrt(uncertainty_estimate))