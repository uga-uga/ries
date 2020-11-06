import numpy as np
from scipy.integrate import nquad

def nquad_subintervals(f, x0, x1_to_xn):
    integral = 0.
    uncertainty_estimate = 0.

    for i in range(len(x0)-1):
        nquad_result = nquad(f, [[x0[i], x0[i+1]], *x1_to_xn])
        integral += nquad_result[0]
        uncertainty_estimate += nquad_result[1]*nquad_result[1]

    return (integral, np.sqrt(uncertainty_estimate))