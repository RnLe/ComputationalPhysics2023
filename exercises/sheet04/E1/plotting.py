import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
import numpy as np

rc('text', usetex=True)

# LaTeX in axis-labels
matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

def plot_potential(phi: np.ndarray, axs: plt.Axes, fig: plt.Figure, title: str) -> matplotlib.image.AxesImage:
    # Plot the potential
    im = axs.imshow(phi, cmap='viridis', extent=(0, 1, 0, 1))
    axs.set_xlabel(r"$\tilde{x} = x / \sigma$")
    axs.set_ylabel(r"$\tilde{y} = y / \sigma$")
    fig.colorbar(im, ax= axs, label=r"$\tilde{\phi} = \phi / \gamma$")
    axs.set_title(title)
    return im

def plot_field(E: np.ndarray, axs: plt.Axes, fig: plt.Figure, im: matplotlib.image.AxesImage, N: int, title: str) -> None:
    axs.quiver(E[1], -E[0])
    axs.set_title(title)
    axs.set_aspect("equal")
    fig.colorbar(im, ax=axs, label=r"$\vec{\tilde{E}} = \vec{E}\sigma / \gamma$")
    _set_ticks(axs, N)
    axs.invert_yaxis()

def _set_ticks(axs: plt.Axes, N: int) -> None:
    # Set the ticks (pain)
    x_ticks = np.arange(0, N + N/5, step=N/5)  # Choose 6 ticks
    x_labels = np.linspace(0, 1, num=len(x_ticks))  # Scale them between 0 and 1
    axs.set_xticks(x_ticks)
    axs.set_xticklabels(x_labels.round(2))
    axs.set_yticks(x_ticks)
    axs.set_yticklabels(x_labels.round(2))

    axs.set_xlabel(r"$\tilde{x} = x / \sigma$")
    axs.set_ylabel(r"$\tilde{y} = y / \sigma$")