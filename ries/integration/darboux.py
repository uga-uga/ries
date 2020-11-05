import numpy as np

def darboux(f, x):
    dx = (x - np.roll(x, 1))[1:]
    fx = f(x)
    lower_upper = np.zeros((2, len(fx)-1))
    lower_upper[0] = fx[:-1]
    lower_upper[1] = fx[1:]
    return (
        np.sum(np.min(lower_upper, axis=0)*dx),
        np.sum(np.max(lower_upper, axis=0)*dx),
    )