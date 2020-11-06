import numpy as np

def darboux(f_or_fx, x):
    dx = (x - np.roll(x, 1))[1:]
    if callable(f_or_fx):
        fx = f_or_fx(x)
    else:
        fx = f_or_fx
    lower_upper = np.zeros((2, len(fx)-1))
    lower_upper[0] = fx[:-1]
    lower_upper[1] = fx[1:]
    lower_sum = np.sum(np.min(lower_upper, axis=0)*dx)
    upper_sum = np.sum(np.max(lower_upper, axis=0)*dx)
    return (
        lower_sum,
        abs(upper_sum-lower_sum),
    )