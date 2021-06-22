from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm


def plot_simulation_results(result, title=None):
    """
    Plot simulation results.

    This creates a two-panel view of the simulation results, with time domain
    and frequency domain plotted side-by-side.

    Parameters
    ----------
    result : an instance of solver.IntegrationResult
        simulation results as returned by the solver

    title : str, optional
        optional plot title
    """

    fig = plt.figure(figsize=(6, 4))

    if title:
        plt.suptitle(title)

    # Time domain plot
    plt.subplot(1, 2, 1)
    plt.pcolormesh(
        result.x, result.t, abs(result.u),
        cmap="magma",
        norm=LogNorm(vmin=1E-6))
    plt.xlim(result.x.min(), result.x.max())
    plt.ylim(result.t.min(), result.t.max())
    plt.xlabel(r"Coordinate $x$, a.u.")
    plt.ylabel(r"Time $t$, a.u.")

    # Frequency domain plot
    plt.subplot(1, 2, 2)
    plt.pcolormesh(
        result.k, result.t, abs(result.v),
        cmap="jet",
        norm=LogNorm(vmin=1E-6))
    plt.xlim(result.k.min(), result.k.max())
    plt.ylim(result.t.min(), result.t.max())
    plt.xlabel(r"Wavenumber $k$, a.u.")
    plt.ylabel(r"Time $t$, a.u.")

    plt.tight_layout()

    return fig
