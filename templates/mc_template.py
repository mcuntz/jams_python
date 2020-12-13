#!/usr/bin/env python
"""
usage: mc_template.py [-h] [-p plotname] [-s] [-t outtype] [-u] [-w] [args [args ...]]

This is the Python template for any new program of Matthias Cuntz.

positional arguments:
  args                  No positional arguments.

optional arguments:
  -h, --help            show this help message and exit
  -p plotname, --plotname plotname
                        Name of plot output file for types pdf, html, d3, or
                        plotly, and name basis for type png (default:
                        /Users/cuntz/prog/github/jams_python/jams/mcplot).
  -s, --serif           Use serif font; default sans serif.
  -t outtype, --type outtype
                        Output type is pdf, png, html, d3, or plotly
                        (default: open screen windows).
  -u, --usetex          Use LaTeX to render text in pdf, png and html.
  -w, --white           White lines on transparent or black background;
                        default: black lines on transparent or white background.


History
-------
Written,  Matthias Cuntz, Jul 2012
Modified, Matthias Cuntz, Jul 2013 - optparse->argparse
          Matthias Cuntz, Jul 2013 - extended to be lookup and gallery
          Matthias Cuntz, Dec 2013 - add png support
          Matthias Cuntz, Mar 2014 - split into individual templates
          Matthias Cuntz, Nov 2014 - script -> function
                       - pdf, png, html, or d3
          Matthias Cuntz, Sep 2015 - Serif and sans serif fonts
          Matthias Cuntz, Dec 2015 - white lines on black or
                                     transparent background
          Matthias Cuntz, Mar 2017 - bokeh, plotly
          Matthias Cuntz, Aug 2018 - use jams.plot snippets
          Matthias Cuntz, Nov 2018 - label=str(iplot) to each add_axes
                                     to suppress warning about future changes
          Matthias Cuntz, Apr 2019 - llcspace missing in calls to legend
          Matthias Cuntz, Dec 2020 - use class mcPlot
"""
from __future__ import division, absolute_import, print_function
import numpy as np
from jams import mcPlot


# -------------------------------------------------------------------------
# Class PlotIt based on Matthias' standard plotting class
#

class PlotIt(mcPlot):

    # -------------------------------------------------------------------------
    # init
    #
    def __init__(self, *args, **kwargs):
        """ initialisation """
        super().__init__(*args, **kwargs)
        # nrow, ncol, colours, etc.
        self.set_extra_layout()

    # -------------------------------------------------------------------------
    # special plot layout
    #
    def set_extra_layout(self):
        from jams.color import colours
        # layout and spaces
        self.nrow     = 3     # # of rows of subplots per figure
        self.ncol     = 2     # # of columns of subplots per figure
        self.hspace   = 0.09  # x-space between subplots
        self.vspace   = 0.04  # y-space between subplots
        self.right    = 0.9   # right space on page
        self.textsize = 11    # standard text size
        self.dxabc    = 0.02  # % of (max-min) shift to the right
                              # of left y-axis for a,b,c,... labels
        self.dyabc    = 0.90  # % of (max-min) shift up from lower x-axis
                                 # for a,b,c,... labels
        self.mcol1 = self.fgcolor     # obs line, obs 2000
        self.mcol2 = self.mcols[-3]   # model line
        self.mcol3 = colours('gray')  # obs 1970
        self.mcol4 = self.mcols[0]    # model 2000
        self.mcol5 = self.mcols[2]    # model 1970
        self.lcol1 = self.mcol1
        self.lcol2 = self.mcol2
        self.lcol3 = self.mcol3
        self.lcol4 = self.mcol4
        self.lcol5 = self.mcol5
        # legend
        self.llxbbox    = 1.0   # x-anchor legend bounding box
        self.llybbox    = 1.0   # y-anchor legend bounding box
        self.llhlength  = 1.5   # length of the legend handles
        # legend
        self.dxabc = 0.02  # % of (max-min) shift to the right
                           # of left y-axis for a,b,c,... labels
        self.dyabc = 0.02  # % of (max-min) shift up from lower x-axis
                           # for a,b,c,... labels

    # -------------------------------------------------------------------------
    # read data
    #
    def read_data(self):
        # do something
        nn = 100
        self.dat = np.arange(nn) / float(nn) * 4.*np.pi

    # -------------------------------------------------------------------------
    # plot fig
    #
    def plot_fig_sin(self):
        import matplotlib.pyplot as plt
        from jams import str2tex, position, abc2plot

        self.ifig += 1
        iplot  = 0
        print('    Plot - Fig ', self.ifig)
        fig = plt.figure(self.ifig)

        xlab = str2tex(r'4 $\pi$', usetex=self.usetex)
        ylab = str2tex('sine and cosine function', usetex=self.usetex)
        xlim = None
        ylim = None

        xx  = self.dat
        yy1 = np.sin(xx)
        yy2 = np.cos(xx)

        iplot += 1

        pos = position(self.nrow, self.ncol, iplot,
                       hspace=self.hspace, vspace=self.vspace)
        sub = fig.add_axes(pos, label=str(iplot))

        larr = []
        tarr = []

        tarr = ['sin']
        larr = sub.plot(xx, yy1)
        plt.setp(larr[-1], linestyle='-', linewidth=self.lwidth, marker=None,
                 color=self.lcols[0])

        tarr += ['cos']
        larr += sub.plot(xx, yy2)
        plt.setp(larr[-1], linestyle='-', linewidth=self.lwidth, marker=None,
                 color=self.lcols[-3])

        if xlab != '':
            plt.setp(sub, xlabel=xlab)
        if ylab != '':
            plt.setp(sub, ylabel=ylab)
        sub.grid(False)

        sub.spines['right'].set_color('none')
        sub.spines['top'].set_color('none')

        if xlim is not None:
            plt.setp(sub, xlim=xlim)
        if ylim is not None:
            plt.setp(sub, ylim=ylim)

        ll = sub.legend(larr, tarr, frameon=self.frameon, ncol=1,
                        labelspacing=self.llrspace,
                        handletextpad=self.llhtextpad,
                        handlelength=self.llhlength,
                        loc='upper left',
                        bbox_to_anchor=(self.llxbbox, self.llybbox),
                        scatterpoints=1, numpoints=1)
        plt.setp(ll.get_texts(), fontsize='small')

        abc2plot(sub, self.dxabc, self.dyabc, iplot, lower=True,
                 bold=True, usetex=self.usetex, mathrm=True)

        # import pdb
        # pdb.set_trace()

        self.plot_save(fig)


if __name__ == '__main__':

    import time as ptime
    t1 = ptime.time()

    desc   = "This is the Python template for any new program"
    desc  += " of Matthias Cuntz."
    argstr = "No positional arguments."
    iplot = PlotIt(desc=desc, argstr=argstr)

    # # Uncomment for xkcd-style
    # import matplotlib.pyplot as plt
    # plt.xkcd()
    iplot.read_data()
    iplot.plot_fig_sin()
    iplot.close()

    t2    = ptime.time()
    strin = ( '[m]: {:.1f}'.format((t2 - t1) / 60.)
              if (t2 - t1) > 60.
              else '[s]: {:d}'.format(int(t2 - t1)) )
    print('    Time elapsed', strin)
