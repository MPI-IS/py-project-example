import numpy


def sech(x):
    """
    Calculates hyperbolic secant of the input without overflowing.

    Parameters
    ----------
    x : array_like
        input array

    Returns
    -------
    sech(x) : array_like
        an element-wise hyperbolic secant
    """
    e = numpy.exp(-abs(x))
    return 2*e / (1 + e**2)
