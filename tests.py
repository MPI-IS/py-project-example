from unittest import TestCase

from matplotlib import pyplot as plt
import numpy

from nlse.helpers import sech
from nlse.plotter import plot_simulation_results
from nlse.solver import gnlse


class SolitonTestCase(TestCase):
    """
    In this test case we check our solver on various known problems of
    soliton propagation. Those are mostly the input/output equality
    tests.
    """

    def test_fundamental_nlse_soliton_propagation(self):
        """
        The most straightforward test case. Take a classic NLSE and
        launch a fundamental soliton. Verify that the output envelope
        is the same as the input envelope.
        """

        def sod(f):
            return - 1/2 * f**2

        def gamma(f):
            return 1

        def kerr(t, x, u):
            return abs(u)**2 * u

        t = numpy.linspace(0, 10, 1000)
        x = numpy.linspace(-20, +20, 1000)
        u0 = sech(x)

        result = gnlse(t, x, u0, sod, gamma, kerr)
        self.assertTrue(result.successful)

        plot_simulation_results(result, title="Fundamental soliton")
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        u = result.u
        self.assertTrue(all(u[0, :] == u0))
        self.assertLess((abs(u[0, :]) - abs(u[-1, :])).max(), 1E-3)

    def test_second_order_soliton_propagation(self):
        """
        This is a slightly more complicated test. We take a
        second-order soliton and let it propagate for 10 periods. We
        check that the soliton shape is restored after each period.
        """

        def sod(f):
            return - 1/2 * f**2

        def gamma(f):
            return 1

        def kerr(t, x, u):
            return abs(u)**2 * u

        period = numpy.pi / 2

        t = numpy.linspace(0, 10 * period, 1000)  # one period = 100 points
        x = numpy.linspace(-20, +20, 1000)
        u0 = 2 / numpy.cosh(x)

        result = gnlse(t, x, u0, sod, gamma, kerr)
        self.assertTrue(result.successful)

        plot_simulation_results(result, title="Second-order soliton")
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        # Pick the initial profile and check it after each period.
        u = result.u
        for n in range(10):
            un = u[n * 100, :]
            self.assertLess((abs(u0) - abs(un)).max(), 1E-3)

    def test_cherenkov_radiation_in_tod_setting(self):
        """
        In here we test that a simple bright soliton launched into a
        third-order dispersion medium develops a radiation tail. The
        maximum of radiation we expect to see near a certain
        frequency predicted by the resonance condition.
        """

        beta3 = 0.15

        def sod(f):
            return - 1/2 * f**2 + 1/6 * beta3 * f**3

        def gamma(f):
            return 1

        def kerr(t, x, u):
            return abs(u)**2 * u

        t = numpy.linspace(0, 10, 1000)
        x = numpy.linspace(-20, +20, 1000)
        u0 = 3 / numpy.cosh(3*x)

        result = gnlse(t, x, u0, sod, gamma, kerr)
        self.assertTrue(result.successful)

        plot_simulation_results(result, title="Cherenkov radiation")
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        f = result.k
        v = result.v

        # We take a look at the output spectrum
        vout = v[-1, :]

        # We extract the radiation part of the spectrum (it should be
        # somewhere right to f = 10)
        frad = f[f > 10]
        vrad = vout[f > 10]

        # We look for the maximum radiation frequency.
        fidx = numpy.argmax(abs(vrad))
        fmax = frad[fidx]

        # We expect the maximum spectrum to be within 10% from
        # 3/beta3.
        self.assertLess(abs(fmax - 3/beta3) / fmax, 0.1)

    def test_radiation_stays_suppressed_by_a_linear_absorber(self):
        """
        In here we test that Cherenkov radiation from a bright soliton
        stays completely trapped by an absorbing boundary layer and
        that the presence of the boundary layer does not affect
        soliton behaviour.

        The test itself is simple: we first launch the soliton without
        the absorber and the with the absorber.
        """

        beta3 = 0.15

        def sod(f):
            return - 1/2 * f**2 + 1/6 * beta3 * f**3

        def gamma(f):
            return 1

        def kerr(t, x, u):
            return abs(u)**2 * u

        t = numpy.linspace(0, 10, 1000)
        x = numpy.linspace(-20, +20, 1000)
        u0 = 3 * sech(3*x)

        # Additional non-reflective absorbing layer near the
        # boundaries of the computational domain
        absorbing_profile = 50 * sech(x + 18) + 50 * sech(x - 18)

        def absorption(t, x, u):
            return 1j * absorbing_profile * u

        result_without_absorber = gnlse(t, x, u0, sod, gamma, kerr)
        result_with_absorber = gnlse(t, x, u0, sod, gamma, kerr, absorption)
        self.assertTrue(result_without_absorber.successful)
        self.assertTrue(result_with_absorber.successful)

        plot_simulation_results(result_without_absorber, title="Without absorber")
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        plot_simulation_results(result_with_absorber, title="With absorber")
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        uwa = result_without_absorber.u[-1, :]
        uw = result_with_absorber.u[-1, :]

        inner = abs(x) < 10  # inside and relatively far away from the boundary layer
        outer = abs(x) > 18  # outside the boundary layer

        # The difference in the inner region should be relatively
        # small: it's mostly the wrap-around Cherenkov radiation.
        # Let's say that that difference should be less than 5% of the
        # soliton peak.
        self.assertLess(
            abs(uw[inner] - uwa[inner]).max(),
            0.05 * abs(uw[inner]).max())

        # The outer region in the case with absorbing layer should
        # contain approximately nothing. We define nothing as less
        # than 0.1% of the original Cherenkov radiation.
        self.assertLess(
            abs(uw[outer]).max(),
            1E-3 * abs(uwa[outer]).max())
