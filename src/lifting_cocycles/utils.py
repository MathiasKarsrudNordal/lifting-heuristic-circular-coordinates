from matplotlib.axes import Axes
from persim import plot_diagrams


def plot_dgms(
    ax: Axes,
    diagrams,
) -> Axes:

    if isinstance(diagrams, tuple):
        diagrams = list(diagrams)
    plot_diagrams(diagrams, ax=ax, size=100, legend=False)

    for l in list(ax.lines):
        l.set_zorder(-1)
        l.set_color("gray")
        l.set_alpha(0.5)
        l.set_linewidth(2.5)

    ax.locator_params(axis="both", nbins=4)
    ax.legend(
        loc="lower right",
        fontsize=ax.xaxis.label.get_fontsize(),
    )
    return ax
