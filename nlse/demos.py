import logging

import numpy
import matplotlib.pyplot as plt

from nlse.solver import gnlse
from nlse.helpers import sech


logging.basicConfig(level=logging.DEBUG)


def soliton_collision():
    """
    This simulation demonstrates collision of two solitons in a classic NLSE.
    """

    def sod(f):
        return -1/2 * f**2

    def gamma(f):
        return 1

    def kerr(t, x, u):
        return abs(u)**2 * u

    t = numpy.linspace(0, 20, 2**9)
    x = numpy.linspace(-20, +20, 2**10)

    u0 = (
        sech(x - 10) * numpy.exp(-1j * 3.0 * x) +
        sech(x + 10) * numpy.exp(+1j * 1.0 * x))
    result = gnlse(t, x, u0, sod, gamma, kerr)

    plt.figure(figsize=(4, 6))

    plt.pcolormesh(
        x, t, abs(result.u),
        cmap="magma")
    plt.xlabel(r"Coordinate $x$, a.u.")
    plt.ylabel(r"Time $t$, a.u.")
    plt.tight_layout()

    plt.show()
