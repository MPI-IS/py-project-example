from datetime import timedelta
import logging
import time

from humanize import precisedelta
from scipy import fft
from scipy import integrate
import numpy


def _report_progress(n, nt, timer_start):
    """
    Report integration progress to the log.
    """

    progress = (n / nt)
    timer_current = time.time()
    timer_elapsed = timer_current - timer_start
    timer_left = timer_elapsed * (1 - progress) / progress

    # Format timers to human-readable strings (hence 'hr_').
    hr_elapsed = precisedelta(
        timedelta(seconds=timer_elapsed),
        minimum_unit="seconds",
        format="%02d")
    hr_left = precisedelta(
        timedelta(seconds=timer_left),
        minimum_unit="seconds",
        format="%02d")

    logging.info(
        "  %5.2f%%  %s elapsed  %s left",
        100 * progress, hr_elapsed, hr_left)


class IntegrationResult:
    """
    This is a container class for the integration results.

    Parameters
    ----------
    t : array_like
        an array of time points where the integrated function is
        evaluated
    x : array_like
        the coordinate grid
    k : array_like
        the frequency grid
    u : array_like of shape len(t)×len(x)
        integrated field in time domain evaluated on a grid
    v : array_like of shape len(t)×len(x)
        integrated field in frequency domain evaluated on a grid
    successful : bool
        a flag indicating that the integration was successful
    error_code : int
        return code of the solver if the integration failed

    Note
    ----
    `u` and `v` are output matrices of the integrated field in
    coordinate-domain and spectral-domain representations. Row number
    is position in the time array, column number is the position in
    coordinate/frequency array.
    """

    def __init__(self, t, x, k, u, v, error_code=None):
        self.t = t
        self.x = x
        self.k = fft.fftshift(k)
        self.u = u
        self.v = v
        self.error_code = error_code

    @property
    def successful(self):
        return self.error_code is None


def gnlse(t, x, u0, beta, gamma, nonlin, lin=None, dt=None, gpc=None):
    """
    Integrate a GNLSE using the integrating factor method.

    This function integrates a generalized version of nonlinear Schrödinger
    equation

        ∂ₜ ũ = i β(k) ũ(t, k)
             + i γ(k) F{ N(t, x, u(t, x)) }
             + i F{ L(t, x, u(t, x)) },

    where ũ(t, k) is the spectrum of the unknown field, β(k) is a
    dispersive operator, and γ(k) is a gain coefficient. u(t, x) is the
    coordinate-domain representation of the same unknown field, and
    nonolinear operator N(t, x, u(t, x)) is defined in terms of that
    coordinate representation. L(t, x, u(t, x)) is an auxiliary linear
    operator that can be used to introduce an absorbing boundary layer.
    It does not have to have a physical meaning.

    The integration is performed using the integrating factor method as
    described by J.M. Dudley & J.R. Taylor in Chapter 3 of Supercontinuum
    Generation in Optical Fibers, CUP 2010. Instead of integrating the
    original equation we resort to integrating a modified version

        ∂ₜ v = i γ(k) F{ N(t, x, u(t, x)) }
             + i F{ L(t, x, u(t, x)) },

    where v = v(t, k) is the modified spectrum that is defined as

        ũ(t, k) = exp(i β(k) t) v(t, k).

    The modifed equation is supposed to be non-stiff, which allows us to
    apply almost any third-party solver. We chose a scipy-provided wrapper
    of ZVODE solver from ODEPACK, which offers error control and an
    adaptive scheme for step-size selection.

    Parameters
    ----------
    t : array_like
        an array of time points where the integrated function should
        be evaluated
    x : array_like
        the coordinate grid
    u0 : array_like
        the initial condition at t[0]
    beta : callable with the signature of beta(f)
        the dispersive profile as a function of frequency
    gamma : callable with the signature of gamma(f)
        frequency-dependent gain part of the nonlinear operator
    nonlin : callable with the signature of nonlin(t, x, u)
        time-domain part of the nonlinear operator
    lin : callable with the signature of lin(t, x, u)
        time-domain linear operator
    gpc : callable with the signature of gcp(t, x, f, u, v)
        an optional callback executed after computing up to the next
        grid point

    Returns
    -------
    result : an instance of IntegrationResult
    """

    # Pre-allocate the output matrices.
    nt = len(t)
    nx = len(x)

    u = numpy.zeros((nt, nx), dtype=complex)
    v = numpy.zeros((nt, nx), dtype=complex)

    # Put the initial conditions in time and frequency domains into
    # the output matrices
    v0 = fft.ifft(u0)
    u[0, :] = u0
    v[0, :] = fft.fftshift(v0)

    # Prepare the frequency scale and evaluate beta on the scale
    k = 2 * numpy.pi * fft.fftfreq(nx, x[1] - x[0])
    D = beta(k)
    G = gamma(k)

    # Prepare the RHS we feed to the solver
    def rhs(t_, v_):
        # Scale the spectrum by the accumulated phase shift due to
        # the dispersion and transform to coordinate space
        exp = numpy.exp(1j * D * (t_ - t[0]))
        u_ = fft.fft(exp * v_)

        # Apply nonlinear operator N() and, maybe, linear operator L()
        # as well, transform back to the modified spectrum and return.
        ret = G * fft.ifft(nonlin(t, x, u_))
        if lin:
            ret += fft.ifft(lin(t, x, u_))

        return 1j / exp * ret

    # Configure the solver. We pick ZVODE from ODEPACK which provides
    # error control and adaptive stepping out of the box. Since we are
    # working on a modified spectrum, we can hope the problem to be
    # non-stiff, but it's a good idea *not* to impose the choice of
    # the integration method and error tolerances -- this will be
    # estimated by the solver itself.
    ode = integrate.ode(rhs)
    ode.set_integrator("zvode", rtol=1E-6)

    # Those are the internal loop variables. `t_` holds the current
    # integration time, `v_` is the current modified spectrum.
    t_ = t[0]
    v_ = v0

    timer_start = time.time()
    logging.info("Integrating:")

    for n in range(1, nt):
        _report_progress(n, nt, timer_start)

        # Pick a step size. THIS IS NOT THE ADAPTIVE STEP-SIZE CONTROL
        # --- it is done by the solver. But sometimes the time step
        # given by the time grid is too coarse for the given precision
        # and the solver exhaust maximum number of integration steps
        # while integrating the equation over the grid step. Providing
        # a finer time grid will result in larger output matrices, and
        # for long simulations those might not fit in the memory
        # anymore. Therefore we provide an additional work-around:
        # user can specify an intermediate time step that is used
        # internally by the loop but not saved anywhere.
        dt_ = dt or t[n] - t_

        while t_ < t[n]:
            # Set the initial value for the current integration step.
            ode.set_initial_value(v_, t_)

            # Pick the next step -- either dt_ or until the next grid
            # point.
            if t_ + dt_ > t[n]:
                t_ = t[n]
            else:
                t_ += dt_

            # Integrate for one step and check the solver status.
            v_ = ode.integrate(t_)
            if not ode.successful():
                return IntegrationResult(
                    t, x, k, u, v, error_code=ode.get_return_code())

        # Calculate the proper spectrum and then save the spectrum and
        # the coordinate representation into the output matrices.
        exp = numpy.exp(1j * D * (t_ - t[0]))
        u[n, :] = fft.fft(exp * v_)
        v[n, :] = fft.fftshift(exp * v_)

        if gpc:
            gpc(t, x, k, u, v)

    logging.info("Done!")
    return IntegrationResult(t, x, k, u, v)
