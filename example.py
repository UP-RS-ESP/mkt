"""
Examples to help illustrate the module BPL.PY
=============================================

Created: Mon Apr 17, 2017  01:19PM
Last modified: Mon Apr 17, 2017  04:46PM
Copyright: Bedartha Goswami <goswami@uni-potsdam.de>

"""

import numpy as np
import matplotlib.pyplot as pl
import mkt


def show_examples():
    """
    Returns the MK tset results for artificial data.
    """
    # create artificial time series with trend
    n = 100
    C = [0.01, 0.001, -0.001, -0.01]
    e = 1.00
    t = np.linspace(0., 500, n)

    # set up figure
    fig, axes = pl.subplots(nrows=2, ncols=2, figsize=[16.00, 9.00])

    # loop through various values of correlation
    ALPHA = 0.01
    for c, ax in zip(C, axes.flatten()):
        # estimate the measurements 'x'
        x = c * t +  e * np.random.randn(n)
        x = np.round(x, 2)

        # get the slope, intercept and pvalues from the mklt module
        MK, m, c, p = mkt.test(t, x, precx=1E-3, alpha=ALPHA, Ha="upordown")

        # plot results
        ax.plot(t, x, "k.-", label="Sampled time series")
        ax.plot(t, m * t + c, "r-", label="Linear fit")
        ax.set_title(MK.upper() + "\np=%.3f, alpha = %.2f" % (p, ALPHA),
                     fontweight="bold", fontsize=10)

        # prettify
        if ax.is_last_row():
            ax.legend(loc="upper right")
            ax.set_xlabel("Time")
        if ax.is_first_col():
            ax.set_ylabel(r"Measurements $x$")
        if ax.is_first_row():
            ax.legend(loc="upper left")

    # save/show plot
    pl.show(fig)
    return None

if __name__ == "__main__":
    print("running example...")
    show_examples()
    print("done.")

