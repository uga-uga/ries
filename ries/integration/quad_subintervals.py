from scipy.integrate import quad

def quad_subintervals(f, x):
    integral = 0.

    for i in range(len(x)-1):
        integral += quad(f, x[i], x[i+1])[0]

    return integral