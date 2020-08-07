#!/usr/bin/env python
################################################
#
# Plot Template:
#
# Description:
#    This is the python plot template of Stephan Thober, which builds on top of that by Matthias Cuntz.
#    It provides quick view on screen or PDF or PNG output in object oriented fashion that can be easily modified for a specific purpose.
#    It is designed for publication-ready PDF output so screen output is not optimal.
#
# License: see license statement below.
#
# author: Stephan Thober
# created: Mar 2019
#
################################################

import pandas as pd
import os
from os.path import isfile
import xarray as xr
import numpy as np
from ufz import position, get_brewer
import argparse
import textwrap
import ufz
import matplotlib.pyplot as plt
""" 
License
-------
This file is part of the JAMS Python package, distributed under the MIT License.

Copyright (c) 2016 Stephan Thober

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
"""


# object to hold data for plotting
class plotObject(object):

    def __init__(self, nrows=1, ncols=1):
        self.nrows = nrows
        self.ncols = ncols
        
        self.getArgs()

        self.readData('test.nc')

        self.setPlotVariables()
        self.plot()
        self.closePlot()


    def plot(self):
        # plot the data
        print('start plotting...')
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import matplotlib.colors as mcolors
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        # plot
        plotData = self.plotData.variables['plotdata']
        
        projection = ccrs.PlateCarree()

        cmap_blue = get_brewer('blues9')
        cmap_red = get_brewer('reds9')
        colors = np.vstack((cmap_blue.colors[::-1], cmap_red.colors))
        cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors, N=len(colors))

        levels = np.linspace(plotData.min(), plotData.max(), 10)
        norm = mpl.colors.BoundaryNorm(levels, levels.shape[0])

        fig = plt.figure(1)
        ax = fig.add_axes(position(2,1,1, right=0.8), projection=projection)
        # ax.set_extent((lon.min(), lon.max(), lat.min(), lat.max()))
        ax.set_title('Test title')
        ax.contourf(plotData, cmap=cmap, norm=norm)
        # self.addFeatures(ax)
                
        # draw a legend
        ax_legend = fig.add_axes(position(1, 1, 1, left=0.85))
        mpl.colorbar.ColorbarBase(ax_legend, cmap=cmap, norm=norm, orientation='vertical')

        # --------------------------------------------------------------------
        # save ---------------------------------------------------------------
        if (self.outtype == 'pdf'):
            self.pdf_pages.savefig(fig)
            plt.close(fig)
        elif (self.outtype == 'png'):
            fig.savefig(pngbase, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
            plt.close(fig)
        else:
            plt.show()     

            
    def readData(self, saveFile):

        # check whether saveFile exists
        # if so read plotData
        if isfile(saveFile):
            self.plotData = xr.open_dataset(saveFile)

            return

        # start reading data
        # this is very flexible and needs to be adjusted
        # for each application
        self.plotData = xr.Dataset({'plotdata': (['x', 'y'],  np.arange(100).reshape(10, 10))})
        #
        # store data to netcdf file
        self.plotData.to_netcdf(saveFile)


    def setPlotVariables(self):
        
        if (self.pdffile != ''):
            self.outtype = 'pdf'
        elif (self.pngbase != ''):
            self.outtype = 'png'
        else:
            self.outtype = 'x'

        # Main plot
        nrow        = 3           # # of rows of subplots per figure
        ncol        = 2           # # of columns of subplots per figure
        hspace      = 0.10        # x-space between subplots
        vspace      = 0.05        # y-space between subplots
        textsize    = 13          # standard text size
        dxabc       = 0.90        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
        dyabc       = 0.05        # % of (max-min) shift up from lower x-axis for a,b,c,... labels

        lwidth      = 1.5         # linewidth
        elwidth     = 1.0         # errorbar line width
        alwidth     = 1.0         # axis line width
        msize       = 1.0         # marker size
        mwidth      = 1.0         # marker edge width
        mcol1       = ufz.color.colours('red')        # primary marker colour
        mcol2       = '0.0'                     # secondary
        mcol3       = (202/255.,0/255.,32/255.) # third
        mcols       = ufz.color.colours(['blue','red','darkgray','orange','darkblue','black'])
        lcol1       = ufz.color.colours('blue')   # primary line colour
        lcol2       = '0.0'
        lcol3       = '0.0'
        lcols       = mcols

        # Map
        delon   = 60.        # spacing between longitude labels
        delat   = 60.        # spacing between latitude labels
        xsize   = textsize   # axis label size
        ncolor  = 10         # # of colors in plot
        cbsize  = textsize   # colorbar label size

        # Legend
        llxbbox     = -0.01       # y-anchor legend bounding box
        llybbox     = 0.04        # y-anchor legend bounding box
        llrspace    = 0.          # spacing between rows in legend
        llcspace    = 1.0         # spacing between columns in legend
        llhtextpad  = 0.4         # the pad between the legend handle and text
        llhlength   = 1.5         # the length of the legend handles
        frameon     = False       # if True, draw a frame around the legend. If None, use rc

        # PNG
        dpi         = 300
        transparent = False
        bbox_inches = 'tight'
        pad_inches  = 0

        import matplotlib as mpl
        if (self.outtype == 'pdf'):
            mpl.use('PDF') # set directly after import matplotlib
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages
            # Customize: http://matplotlib.sourceforge.net/users/customizing.html
            mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
            mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
            if self.usetex:
                mpl.rc('text', usetex=True)
            else:
                #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
                mpl.rc('font',**{'family':'serif','serif':['times']})
            mpl.rc('text.latex', unicode=True)
        elif (self.outtype == 'png'):
            mpl.use('Agg') # set directly after import matplotlib
            import matplotlib.pyplot as plt
            mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
            if self.usetex:
                mpl.rc('text', usetex=True)
            else:
                #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
                mpl.rc('font',**{'family':'serif','serif':['times']})
            mpl.rc('text.latex', unicode=True)
            mpl.rc('savefig', dpi=dpi, format='png')
        else:
            import matplotlib.pyplot as plt
            mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
        mpl.rc('font', size=textsize)
        mpl.rc('lines', linewidth=lwidth, color='black')
        mpl.rc('axes', linewidth=alwidth, labelcolor='black')
        mpl.rc('path', simplify=False) # do not remove

        if (self.outtype == 'pdf'):
            print('Plot PDF ', self.pdffile)
            self.pdf_pages = PdfPages(self.pdffile)
        elif (self.outtype == 'png'):
            print('Plot PNG ', self.pngbase)
        else:
            print('Plot X')
        # figsize = mpl.rcParams['figure.figsize']


    def getArgs(self):

        addargs = []
        pngbase = ''
        pdffile = ''
        usetex  = False
        parser  = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\
            description:
            This is the python plot template of Stephan Thober, which builds on top of that by Matthias Cuntz.
            It provides quick view on screen or PDF or PNG output in object oriented fashion that can be easily modified for a specific purpose.
            It is designed for publication-ready PDF output so screen output is not optimal.
            '''))
        parser.add_argument('-g', '--pngbase', action='store',
                                default=pngbase, dest='pngbase', metavar='pngbase',
                                help='Name basis for png output files (default: open screen window).')
        parser.add_argument('-p', '--pdffile', action='store',
                                default=pdffile, dest='pdffile', metavar='pdffile',
                                help='Name of pdf output file (default: open screen window).')
        parser.add_argument('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                                help="Use LaTeX to render text in pdf.")
        parser.add_argument('further', nargs='*', metavar='additional_args',
                            help='Any additional arguments such as filenames, etc.')
        
        args    = parser.parse_args()
        self.addargs = args.further
        self.pdffile = args.pdffile
        self.pngbase = args.pngbase
        self.usetex  = args.usetex

        if (pdffile != '') & (pngbase != ''):
            print('\nError: PDF and PNG are mutually exclusive. Only either -p or -g possible.\n')
            parser.print_usage()
            import sys
            sys.exit()

        del parser, args


    def closePlot(self):
        if (self.outtype == 'pdf'):
            self.pdf_pages.close()
        elif (self.outtype == 'png'):
            pass
        else:
            plt.show()

            
    def addFeatures(self, ax):
        states_provinces = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='50m',
                facecolor='none')

        ax.add_feature(cfeature.LAND)
        ax.add_feature(states_provinces, edgecolor='gray', lw=0.5)
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS)
        # ax.add_feature(cfeature.LAKES, alpha=0.5)


if __name__ == "__main__":
    print('Start')

    plot = plotObject()
    
    print('Done!')
