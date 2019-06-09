#!/usr/bin/env python
from __future__ import print_function
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

# -------------------------------------------------------------------------
# Command line arguments
#

pngbase   = ''
pdffile   = ''
usetex     = False

import argparse
import textwrap
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
          description:
            This is the python plot script for ploting streamflow for Aug 2002 over the EDgE domain and
            adding the location of gauges with size

          Example:
            python plot_Figure_-1.py

          Note:
          '''))

parser.add_argument('-g', '--pngbase', action='store',
                    default=pngbase, dest='pngbase', metavar='pngbase',
                    help='Name basis for png output files (default: open screen window)')
parser.add_argument('-p', '--pdffile', action='store',
                    default=pdffile, dest='pdffile', metavar='pdffile',
                    help='Name of pdf output file (default: open screen window).')
parser.add_argument('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf.")

# evaluate args
args = parser.parse_args()
pngbase   = args.pngbase
pdffile   = args.pdffile
usetex    = args.usetex

# consistency checks
exit = False
if (pdffile != '') & (pngbase != ''):
    print('PDF and PNG are mutually exclusive. Only either -p or -g possible.')
    parser.print_usage()
    exit = True
if exit:        
    import sys
    sys.exit()
del parser, args

# further import statements
import numpy as np

# -------------------------------------------------------------------------
# Customize plots
#
if (pdffile == ''):
    if (pngbase == ''):
        outtype = 'x'
    else:
        outtype = 'png'
else:
    outtype = 'pdf'

# Main plot
nrow       = 1
ncol       = 1
textsize   = 12          # standard text size
hspace     = 0.08        # x-space between subplots
wspace     = 0.04        # y-space between subplots
dxabc      = 0.05        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc      = 0.90        # % of (max-min) shift up from lower x-axis for a,b,c,... labels

lwidth     = 0.5         # linewidth
elwidth    = 1.0         # errorbar line width
alwidth    = 1.0         # axis line width
msize      = 1.0         # marker size
mwidth     = 0.5         # marker edge width
# color: 'b'|'g'|'r'|'c'|'m'|'y'|'k'|'w'
#        'blue'|'green'|'red'|'cyan'|'magenta'|'yellow'|'black'|'white'
#        hex string '#eeefff' | RGB tuple (1,0.5,1) | html names 'burlywod', 'chartreuse', ...
#        grayscale intensity, e.g. '0.7', 'k'='0.0'
myBlue          = '#00589C'
myLightBlue     = '#E5F4FF'
FigLightBlue    = '#6FA6B9'
otherBlue       = '#6FA6FF'     # rgb = 111,166,255
myRed           = '#CA373B'
myOrange        = '#FDAE61'
mcol1      = (202/255.,0/255.,32/255.)     # primary marker colour
mcol2      = '0.0'       # color of second markers
mcol3      = '0.0'       # color of third markers
lcol1      = (5/255.,113/255.,176/255.)       # primary line colour
lcol2      = '0.0'       # color of second lines
lcol3      = '0.0'       # color of third lines

# Legend
llxbbox    = -0.01       # y-anchor legend bounding box
llybbox    = 0.04        # y-anchor legend bounding box
llrspace   = 0.          # spacing between rows in legend
llcspace   = 1.0         # spacing between columns in legend
llhtextpad = 0.4         # the pad between the legend handle and text
llhlength  = 1.5         # the length of the legend handles
frameon    = False       # if True, draw a frame around the legend. If None, use rc
llxbbox2    = 0.60       # Tight bounding of symbol and text (w/o lines)
llhtextpad2 = 0.         #                   "
llhlength2  = 1.0        #                   "

# PNG
dpi         = 300
transparent = False
bbox_inches = 'tight'
pad_inches  = 0
import matplotlib as mpl
#
if (outtype == 'pdf'):
    mpl.use('PDF') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf    
    if usetex:
        mpl.rc('text', usetex=True)
        mpl.rc('text.latex', unicode=True)
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('font', size=textsize)
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    # mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    mpl.rc('figure', figsize=(11.69,8.27)) # a4 landscape
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
    mpl.rc('savefig', dpi=dpi, format='png')
else:
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(14,10))

mpl.rc('font', size=textsize)
mpl.rc('lines', linewidth=lwidth, color='black')
mpl.rc('axes', linewidth=alwidth, labelcolor='black')
mpl.rc('path', simplify=False) # do not remove
# mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
mpl.rc('figure', figsize=(11.69,8.27)) # a4 landscape
#mpl.rc('figure', figsize=(4./5.*12,4./5.*10.5))

if (outtype == 'pdf'):
    print('Plot PDF ', pdffile)
    pdf_pages = PdfPages(pdffile)
elif (outtype == 'png'):
    print('Plot PNG ', pngbase)
else:
    print('Plot X')

figsize = mpl.rcParams['figure.figsize']

ifig = 1
iplot = 0

print( 'Plot - Fig ', ifig)
fig = plt.figure(ifig)

# -----------------------------------------------------------------------------
# PLOT
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# geographical features
ocean         = cfeature.NaturalEarthFeature(category='physical', name='ocean',
                                             scale='50m', facecolor=cfeature.COLORS['water'])
state_borders = cfeature.NaturalEarthFeature(category='cultural', name='admin_0_countries',
                                             scale='50m',
                                             edgecolor='0.2',
                                             facecolor=cfeature.COLORS['land'])
# setup two projections
# projection = ccrs.LambertConformal(central_latitude=52, central_longitude=10)
projection = ccrs.PlateCarree()

ax = plt.axes(projection=projection)
lon, lat = np.meshgrid(np.arange(10), np.arange(50, 60))
ax.plot(lon, lat, marker='.', ms=10, mfc='k', mec='k', lw=0., transform=projection)

# some resources
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()])
ax.add_feature(ocean, zorder=-1)
ax.add_feature(state_borders)

# -----------------------------------------------------------------------------
# save
if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    fig.savefig(pngbase, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)
 
if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()
