#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
from scipy import interpolate, signal

def xkcd_line(x, y, xlim=None, ylim=None, mag=1.0, f1=30, f2=0.05, f3=15):
    """
        Mimic a hand-drawn line from (x, y) data


        Definition
        ----------
        def xkcd_line(x, y, xlim=None, ylim=None, mag=1.0, f1=30, f2=0.05, f3=15):


        Input
        -----
        x, y          array_like; arrays to be modified


        Optional Input
        --------------
        xlim, ylim    data range; the assumed plot range for the modification.
                      If not specified,  they will be guessed from the  data
        mag           float; the magnitude of the distortion (default: 1.0)
        f1, f2, f3    int, float, int; filtering parameters.
                      f1 gives the size of the window (default: 50)
                      f2 gives the high-frequency cutoff (default: 0.01)
                      f3 gives the size of the filter (default: 15)


        Output
        ------
        x, y          ndarrays; the modified lines


        References
        ----------
        See xkcd below.


        Examples
        --------
        for line in ax.lines:
            x, y         = line.get_data()
            x_int, y_int = xkcd_line(x, y, xlim, ylim, mag, f1, f2, f3)
            line.set_data(x_int, y_int)


        History
        -------
        Written,  MC, Mar 2013
    """
    # assure array
    x = np.asarray(x)
    y = np.asarray(y)

    # get limits for rescaling
    if xlim is None: xlim = (x.min(), x.max())
    if ylim is None: ylim = (y.min(), y.max())
    if xlim[1] == xlim[0]: xlim = ylim
    if ylim[1] == ylim[0]: ylim = xlim

    # scale the data
    x_scaled = (x - xlim[0]) * 1. / (xlim[1] - xlim[0])
    y_scaled = (y - ylim[0]) * 1. / (ylim[1] - ylim[0])

    # compute the total distance along the path
    dx       = x_scaled[1:] - x_scaled[:-1]
    dy       = y_scaled[1:] - y_scaled[:-1]
    dist_tot = np.sum(np.sqrt(dx*dx + dy*dy))

    # number of interpolated points is proportional to the distance
    Nu = int(200 * dist_tot)
    u  = np.arange(-1, Nu + 1) * 1. / (Nu - 1)

    # interpolate curve at sampled points
    # k            = min(3, len(x) - 1)
    k            = min(3, x.size - 1)
    res          = interpolate.splprep([x_scaled, y_scaled], s=0, k=k)
    x_int, y_int = interpolate.splev(u, res[0])

    # we perturb perpendicular to the drawn line
    dx   = x_int[2:] - x_int[:-2]
    dy   = y_int[2:] - y_int[:-2]
    # horizontal or vertical lines
    # np.sign(np.cumsum(np.random.random(dx.size)-0.5)) emulates something like a Brownian motion
    # i.e. auto-correlated random walks around 0; just the sign interests here.
    eps = np.maximum(np.abs(np.amax(x_scaled)-np.amin(x_scaled)), np.abs(np.amax(y_scaled)-np.amin(y_scaled)))/Nu
    if np.all(np.abs(dx) < eps):
        dx = np.sign(np.cumsum(np.random.random(dx.size)-0.5)) * eps
    if np.all(np.abs(dy) < eps):
        dy = np.sign(np.cumsum(np.random.random(dx.size)-0.5)) * eps
    # equal distances
    if np.all(np.sign(dx) == np.sign(dx[0])):
        dx *= np.sign(np.cumsum(np.random.random(dx.size)-0.5))
    if np.all(np.sign(dy) == np.sign(dy[0])):
        dy *= np.sign(np.cumsum(np.random.random(dx.size)-0.5))
    dist = np.sqrt(dx * dx + dy * dy)

    # create a filtered perturbation
    # coeffs       = mag * np.random.normal(0, 0.01, len(x_int) - 2)
    coeffs       = mag * np.random.normal(0, 0.01, x_int.size - 2)
    b            = signal.firwin(f1, f2*dist_tot, window=('kaiser', f3))
    response     = signal.lfilter(b, 1, coeffs)
    x_int[1:-1] += response * dy / dist
    y_int[1:-1] += response * dx / dist

    # un-scale data
    x_int = x_int[1:-1] * (xlim[1] - xlim[0]) + xlim[0]
    y_int = y_int[1:-1] * (ylim[1] - ylim[0]) + ylim[0]

    return x_int, y_int


def xkcd(ax,
         mag=1.0,
         f1=50, f2=0.01, f3=15,
         bgcolor='w',
         title_size=None,
         xaxis_loc=None, yaxis_loc=None,
         xaxis_arrow='+', yaxis_arrow='+',
         ax_extend=0.1,
         xlabel_inside=0., ylabel_inside=0.,
         ticks=False,
         xticks_inside=0., yticks_inside=0.,
         ):
    """
        Make axis look hand-drawn

        This adjusts all lines, text, legends, and axes in the figure to look
        like xkcd plots, a webcomic from Randall Munroe. Other plot elements are not modified.


        Definition
        ----------
        def xkcd(ax,
                 mag=1.0,
                 f1=50, f2=0.01, f3=15,
                 bgcolor='w',
                 title_size=None,
                 xaxis_loc=None, yaxis_loc=None,
                 xaxis_arrow='+', yaxis_arrow='+',
                 ax_extend=0.1,
                 xlabel_inside=0., ylabel_inside=0.,
                 ticks=False,
                 xticks_inside=0., yticks_inside=0.,
                 ):


        Input
        -----
        ax    Axes instance the axes instance to be modified.


        Optional Input
        --------------

        mag                            float; the magnitude of the distortion (default: 1.0)
        f1, f2, f3                     int, float, int; filtering parameters.
                                       f1 gives the size of the window (default: 50)
                                       f2 gives the high-frequency cutoff (default: 0.01)
                                       f3 gives the size of the filter (default: 15)
        bgcolor                        str; color around lines so that axis look brocken,
                                       i.e. lines are overdrawn on axis (default: 'w')
        titel_size                     float; poitn size of plot title. If None, same size as axis labels.
                                       (default: None)
        xaxis_loc, yaxis_log           float; The locations to draw the x and y axes in data coordinates.
                                       If not specified, they will be drawn from the bottom left of the plot.
                                       (default: None)
        xaxis_arrow, yaxis_arrow       str; where to draw arrows on the x/y axes
                                       Options are '+', '-', '+-', or '' (default: '+')
        ax_extend                      float; How far (fractionally) to extend the drawn axes beyond
                                       the original axes limits (default: 0.1)
        xlabel_inside, ylabel_inside   float: By how much the labels are shifted (default: 0.0)

        The last two options are not working how with mc_plot_template
        ticks                          True: change tick labels; False: no tick labels are drawn (default: False)
        xticks_inside, yticks_inside   float: By how much the ticks are shifted (default: 0.0)


        Output
        ------
        ax is basically empty and all former elements are redrawn on plot.


        Note
        ----
        For reproducible plots, seed the random number generator before each new plot.
        If a new line was added, the old lines will look the same. The legend will be different though.


        References
        ----------
        This is the modified XKCD plot generator of Jake Vanderplas
        http://nbviewer.ipython.org/url/jakevdp.github.com/downloads/notebooks/XKCD_plots.ipynb

        The idea for this comes from work by Damon McDougall
        http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg25499.html


        Examples
        --------
        import matplotlib.pylab as plt
        fig = plt.figure(1)
        ax = fig.add_axes([0.1,0.1,0.5,0.5])
        ax.plot(range(10), label='Line')
        ax.set_title('Title')
        ax.set_xlabel('x label')
        ax.set_ylabel('y label')
        ax.legend()
        xkcd(ax)


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2013 Matthias Cuntz - mc (at) macu (dot) de

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.


        History
        -------
        Written,  MC, Mar 2013
    """
    import matplotlib.pylab as plt
    import matplotlib.font_manager as fm

    # remember random state for later resetting
    random_state = np.random.get_state()

    # Get axes aspect
    ext = ax.get_window_extent().extents
    aspect = (ext[3] - ext[1]) / (ext[2] - ext[0])

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    xspan = xlim[1] - xlim[0]
    yspan = ylim[1] - xlim[0]

    xax_lim = (xlim[0] - ax_extend * xspan,
               xlim[1] + ax_extend * xspan)
    yax_lim = (ylim[0] - ax_extend * yspan,
               ylim[1] + ax_extend * yspan)

    if xaxis_loc is None: xaxis_loc = ylim[0]
    if yaxis_loc is None: yaxis_loc = xlim[0]

    # Draw axes
    acolor  = ax.get_xaxis().get_gridlines()[0].get_color()
    xaxis = plt.Line2D([xax_lim[0], xax_lim[1]], [xaxis_loc, xaxis_loc],
                      linestyle='-', color=acolor)
    yaxis = plt.Line2D([yaxis_loc, yaxis_loc], [yax_lim[0], yax_lim[1]],
                      linestyle='-', color=acolor)

    # adjust the axes
    if ticks:
        for x, xtext in zip(ax.get_xticks(), ax.get_xticklabels()):
            ax.text(x, xaxis_loc - 0.08 * yspan * (2 * xticks_inside - 1), xtext.get_text(),
                fontsize=xtext.get_size(), ha='center', va='bottom' if xticks_inside else 'top', rotation=0)
        for y, ytext in zip(ax.get_yticks(), ax.get_yticklabels()):
            ax.text(yaxis_loc + 0.02 * xspan * (2 * yticks_inside - 1), y, ytext.get_text(),
                fontsize=ytext.get_size(), ha='left' if yticks_inside else 'right', va='center', rotation=0)

    # Label axes
    siz = ax.get_xaxis().get_label().get_size()
    ax.text(xax_lim[1], xaxis_loc - 0.2 * yspan * (2 * xlabel_inside - 1), ax.get_xlabel(),
            fontsize=siz, ha='right', va='bottom' if xlabel_inside else 'top', rotation=0)
    ax.text(yaxis_loc + 0.04 * xspan * (2 * ylabel_inside - 1), yax_lim[1], ax.get_ylabel(),
            fontsize=siz, ha='right', va='bottom' if ylabel_inside else 'top', rotation=84)

    # Title - default: same size as axis labels
    if title_size is not None:
        siz2 = title_size
    else:
        siz2 = siz
    ax.text(0.5 * (xax_lim[1] + xax_lim[0]), yax_lim[1], ax.get_title(),
            ha='center', va='bottom', fontsize=siz2)

    # Draw arrow-heads at the end of axes lines
    arr1       = 0.04 * np.array([-1, 0, -1])
    arr2       = 0.03 * np.array([-1, 0, 1])
    arr1[::2] += np.random.normal(0, 0.005 / 2, 2)
    arr2[::2] += np.random.normal(0, 0.005 / 2, 2)
    x, y = xaxis.get_data()
    if '+' in str(xaxis_arrow):
        ax.plot(x[-1] + arr1 * xspan * aspect,
                y[-1] + arr2 * yspan,
                color=acolor, lw=2)
    if '-' in str(xaxis_arrow):
        ax.plot(x[0] - arr1 * xspan * aspect,
                y[0] - arr2 * yspan,
                color=acolor, lw=2)

    x, y = yaxis.get_data()
    if '+' in str(yaxis_arrow):
        ax.plot(x[-1] + arr2 * xspan * aspect**2,
                y[-1] + arr1 * yspan / aspect,
                color=acolor, lw=2)
    if '-' in str(yaxis_arrow):
        ax.plot(x[0] - arr2 * xspan * aspect**2,
                y[0] - arr1 * yspan / aspect,
                color=acolor, lw=2)

    # Set the axis limits
    ax.set_xlim(xax_lim[0] - 0.1 * xspan,
                xax_lim[1] + 0.1 * xspan)
    ax.set_ylim(yax_lim[0] - 0.1 * yspan,
                yax_lim[1] + 0.1 * yspan)

    # The lines
    Nlines = len(ax.lines)
    lines  = [xaxis, yaxis] + [ax.lines.pop(0) for i in range(Nlines)]
    for line in lines:
        x, y = line.get_data()
        ls   = line.get_linestyle()
        if ls != 'None':
            x_int, y_int = xkcd_line(x, y, xlim, ylim, mag, f1, f2, f3)
        else:
            x_int, y_int = x, y
        # create foreground and background line
        lw = line.get_linewidth()
        line.set_linewidth(2*lw)
        line.set_data(x_int, y_int)

        # White surrounding of line makes them look overplot on axis
        if (line is not xaxis) and (line is not yaxis) and ls != 'None':
            line_bg = plt.Line2D(x_int, y_int, color=bgcolor, linewidth=2*lw+4)
            ax.add_line(line_bg)
        ax.add_line(line)

    # Change all the fonts to humor-sans.
    # from jams.find_in_path import find_in_path
    # fhumor = find_in_path('Humor-Sans.ttf') # in jams_python
    import os
    fhumor = os.path.join(os.path.dirname(__file__), 'Humor-Sans.ttf') # in jams_python/jams
    for text in ax.texts:
        tsize  = text.get_size()
        prop = fm.FontProperties(fname=fhumor, size=tsize)
        text.set_fontproperties(prop)

    # modify legend
    leg = ax.get_legend()
    if leg is not None:
        np.random.set_state(random_state) # restate random number generator for reproducible results
        leg.set_frame_on(False)
        for child in leg.get_children():
            if isinstance(child, plt.Line2D):
                x, y = child.get_data()
                child.set_data(xkcd_line(x, y, mag=10.*mag, f1=2*f1, f2=f2/10.))
                child.set_linewidth(2*child.get_linewidth())
            if isinstance(child, plt.Text):
                tsize  = child.get_size()
                prop = fm.FontProperties(fname=fhumor, size=tsize)
                child.set_fontproperties(prop)

    # remove standard axis
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()

    return ax


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import numpy as np
    # import jams
    # from position import position
    # pdffile='test_xkcd.pdf'
    # usetex  = False
    # textsize    = 13          # standard text size
    # lwidth      = 1.5         # linewidth
    # alwidth     = 1.0         # axis line width

    # if (pdffile == ''):
    #     outtype = 'x'
    # else:
    #     outtype = 'pdf'

    # import matplotlib as mpl
    # if (outtype == 'pdf'):
    #   mpl.use('PDF') # set directly after import matplotlib
    #   import matplotlib.pyplot as plt
    #   from matplotlib.backends.backend_pdf import PdfPages
    #   # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    #   mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
    #   mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    #   if usetex:
    #     mpl.rc('text', usetex=True)
    #   else:
    #     #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #     mpl.rc('font',**{'family':'serif','serif':['times']})
    #   mpl.rc('text.latex', unicode=True)
    #   mpl.rc('font', size=textsize)
    # else:
    #   import matplotlib.pyplot as plt
    #   mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
    #   mpl.rc('font', size=textsize)
    # mpl.rc('lines', linewidth=lwidth, color='black')
    # mpl.rc('axes', linewidth=alwidth, labelcolor='black')
    # mpl.rc('path', simplify=False) # do not remove

    # if (outtype == 'pdf'):
    #     print('Plot PDF ', pdffile)
    #     pdf_pages = PdfPages(pdffile)
    # else:
    #     print('Plot X')
    # figsize = mpl.rcParams['figure.figsize']
    # # figsize = [6.616, 9.352]

    # ifig = 0
    # nrow = 2
    # ncol = 1

    # # ----------------------------------------------------------------------------------
    # # Example 1
    # np.random.seed(1)

    # iplot = 1
    # fig = plt.figure(ifig)
    # pos = position(nrow,ncol,iplot,golden=False,figsize=figsize,left=0.1)
    # ax = fig.add_axes(pos)

    # x = np.linspace(0, 10, 100)
    # ax.plot(x, np.sin(x) * np.exp(-0.1 * (x - 5) ** 2), 'b', lw=1, label='sine')
    # ax.plot(x, -np.cos(x) * np.exp(-0.1 * (x - 5) ** 2), 'r', lw=1, label='cosine')
    # ax.set_title('check it out!')
    # ax.set_xlabel('x label')
    # ax.set_ylabel('y label')
    # ax.legend(loc='upper left', bbox_to_anchor=(0.7,0.4), ncol=1, handlelength=0)

    # xkcd(ax, xaxis_loc=0.0, yaxis_loc=1.0,
    #      xaxis_arrow='+-', yaxis_arrow='+-', xlabel_inside=1., title_size=textsize+2)

    # if (outtype == 'pdf'):
    #     pdf_pages.savefig(fig)
    #     plt.close()


    # # ----------------------------------------------------------------------------------
    # # Example 1 with third line
    # np.random.seed(1)

    # iplot = 1
    # fig = plt.figure(ifig)
    # pos = position(nrow,ncol,iplot,golden=False,figsize=figsize,left=0.1)
    # ax = fig.add_axes(pos)

    # x = np.linspace(0, 10, 100)
    # ax.plot(x, np.sin(x) * np.exp(-0.1 * (x - 5) ** 2), 'b', lw=1, label='sine')
    # ax.plot(x, -np.cos(x) * np.exp(-0.1 * (x - 5) ** 2), 'r', lw=1, label='cosine')
    # ax.plot(x, -np.cos(x+1.0) * np.exp(-0.1 * (x - 5) ** 2), 'm', lw=1, label='shift')
    # ax.set_title('check it out!')
    # ax.set_xlabel('x label')
    # ax.set_ylabel('y label')
    # ax.legend(loc='upper left', bbox_to_anchor=(0.7,0.4), ncol=1, handlelength=0)

    # xkcd(ax, xaxis_loc=0.0, yaxis_loc=1.0,
    #      xaxis_arrow='+-', yaxis_arrow='+-', xlabel_inside=1., title_size=textsize+2)

    # if (outtype == 'pdf'):
    #     pdf_pages.savefig(fig)
    #     plt.close()

    # # ----------------------------------------------------------------------------------
    # # Example 2

    # # Some helper functions
    # def norm(x, x0, sigma):
    #     return np.exp(-0.5 * (x - x0) ** 2 / sigma ** 2)

    # def sigmoid(x, x0, alpha):
    #     return 1. / (1. + np.exp(- (x - x0) / alpha))

    # # define the curves
    # x = np.linspace(0, 1, 100)
    # y1 = np.sqrt(norm(x, 0.7, 0.05)) + 0.2 * (1.5 - sigmoid(x, 0.8, 0.05))

    # y2 = 0.2 * norm(x, 0.5, 0.2) + np.sqrt(norm(x, 0.6, 0.05)) + 0.1 * (1 - sigmoid(x, 0.75, 0.05))

    # y3 = 0.05 + 1.4 * norm(x, 0.85, 0.08)
    # y3[x > 0.85] = 0.05 + 1.4 * norm(x[x > 0.85], 0.85, 0.3)

    # ifig  += 1
    # iplot = 1
    # fig = plt.figure(ifig)
    # ax = fig.add_axes(position(nrow,ncol,iplot,golden=False,figsize=figsize,left=0.1))

    # # draw the curves
    # ax.plot(x, y1, c='gray')
    # ax.plot(x, y2, c='blue')
    # ax.plot(x, y3, c='red')

    # ax.text(0.3, -0.1, "Yard")
    # ax.text(0.5, -0.1, "Steps")
    # ax.text(0.7, -0.1, "Door")
    # ax.text(0.9, -0.1, "Inside")

    # ax.text(0.05, 1.1, "fear that\nthere's\nsomething\nbehind me")
    # ax.plot([0.15, 0.2], [1.0, 0.2], '-k', lw=0.5)

    # ax.text(0.25, 0.8, "forward\nspeed")
    # ax.plot([0.32, 0.35], [0.75, 0.35], '-k', lw=0.5)

    # ax.text(0.9, 0.4, "embarrassment")
    # ax.plot([0.8, 1.0], [1.05, 0.55], '-k', lw=0.5)

    # ax.set_title("Walking back to my\nfront door at night:")

    # ax.set_xlim(0, 1)
    # ax.set_ylim(0, 1.5)

    # # modify all the axes elements in-place
    # xkcd(ax)

    # if (outtype == 'pdf'):
    #     pdf_pages.savefig(fig)
    #     plt.close()


    # if (outtype == 'pdf'):
    #     pdf_pages.close()
    # else:
    #     plt.show()
