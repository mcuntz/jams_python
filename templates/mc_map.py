#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
usage: mc_map.py [-h] [-g pngbase] [-p pdffile] [-t] [Input file]

This is the Python template for drawings maps of Matthias Cuntz: use with input file mc_map-foster_davy-snow_depth.nc

positional arguments:
  Input file            Mandatory input file.

optional arguments:
  -h, --help            show this help message and exit
  -g pngbase, --pngbase pngbase
                        Name basis for png output files (default: open screen
                        window).
  -p pdffile, --pdffile pdffile
                        Name of pdf output file (default: open screen window).
  -t, --usetex          Use LaTeX to render text in pdf.


License
-------
This file is part of the JAMS Python package, distributed under the MIT License.

Copyright (c) 2012-2019 Matthias Cuntz - mc (at) macu (dot) de

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
Written,  MC, Jul 2012
Modified, MC, Jul 2013 - optparse->argparse
          MC, Jul 2013 - extended to be lookup and gallery
          MC, Dec 2013 - add png support
          MC, Mar 2014 - split into individual templates
"""

# -------------------------------------------------------------------------
# Command line arguments
#

import argparse

pngbase = ''
pdffile = ''
usetex  = False
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''This is the Python template for drawings maps of Matthias Cuntz: use with input file mc_map-foster_davy-snow_depth.nc''')
parser.add_argument('-g', '--pngbase', action='store',
                    default=pngbase, dest='pngbase', metavar='pngbase',
                    help='Name basis for png output files (default: open screen window).')
parser.add_argument('-p', '--pdffile', action='store',
                    default=pdffile, dest='pdffile', metavar='pdffile',
                    help='Name of pdf output file (default: open screen window).')
parser.add_argument('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf.")
parser.add_argument('file', nargs='?', default=None, metavar='Input file',
                   help='Mandatory input file.')

args    = parser.parse_args()
pngbase = args.pngbase
pdffile = args.pdffile
usetex  = args.usetex
infile  = args.file

if (infile is None):
    print('\nError: Input file must be given.\n')
    parser.print_usage()
    import sys
    sys.exit()

if (pdffile != '') & (pngbase != ''):
    print('\nError: PDF and PNG are mutually exclusive. Only either -p or -g possible.\n')
    parser.print_usage()
    import sys
    sys.exit()

del parser, args

# import packages after help so that help with command line -h is fast
import numpy as np
import jams
import time
t1 = time.time()

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
nrow        = 1           # # of rows of subplots per figure
ncol        = 1           # # of columns of subplots per figure
hspace      = 0.10        # x-space between subplots
vspace      = 0.05        # y-space between subplots
right       = 0.9         # right space on page
textsize    = 13          # standard text size
dxabc       = 0.90        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc       = 0.05        # % of (max-min) shift up from lower x-axis for a,b,c,... labels

lwidth      = 1.5         # linewidth
elwidth     = 1.0         # errorbar line width
alwidth     = 1.0         # axis line width
msize       = 1.0         # marker size
mwidth      = 1.0         # marker edge width
mcol1       = jams.color.colours('red')        # primary marker colour
mcol2       = '0.0'                     # secondary
mcol3       = (202/255.,0/255.,32/255.) # third
mcols       = jams.color.colours(['blue','red','darkgray','orange','darkblue','black'])
lcol1       = jams.color.colours('blue')   # primary line colour
lcol2       = '0.0'
lcol3       = '0.0'
lcols       = mcols

# Map
delon   = 20.        # spacing between longitude labels
delat   = 10.        # spacing between latitude labels
xsize   = textsize   # axis label size
ncolor  = 10         # # of colors in plot
cbsize  = textsize   # colorbar label size

# Legend
llxbbox     = 0           # x-anchor legend bounding box
llybbox     = 1           # y-anchor legend bounding box
llrspace    = 0.          # spacing between rows in legend
llcspace    = 1.0         # spacing between columns in legend
llhtextpad  = 0.4         # the pad between the legend handle and text
llhlength   = 1.5         # the length of the legend handles
frameon     = False       # if True, draw a frame around the legend. If None, use rc

# PNG
dpi         = 300
transparent = False
bbox_inches = 'tight'
pad_inches  = 0.035

import matplotlib as mpl
if (outtype == 'pdf'):
    mpl.use('PDF') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
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


from mpl_toolkits.basemap import Basemap#, shiftgrid

# -------------------------------------------------------------------------
# Read file
#

undef = 9999.

print('Read input file ', infile)
foster_lat  = jams.readnetcdf(infile,'Y')
foster_lon  = jams.readnetcdf(infile,'X')
foster_snow = jams.readnetcdf(infile,'snow_depth') * 1e-3 # mm -> m

foster_lonh, foster_lath = jams.grid_mid2edge(foster_lon, foster_lat)
foster_snow_max = np.ma.amax(foster_snow, axis=0)

# -------------------------------------------------------------------------
# Plot
#

if (outtype == 'pdf'):
    print('Plot PDF ', pdffile)
    pdf_pages = PdfPages(pdffile)
elif (outtype == 'png'):
    print('Plot PNG ', pngbase)
else:
    print('Plot X')
# figsize = mpl.rcParams['figure.figsize']

ifig = 0

# -------------------------------------------------------------------------
# Fig 1 - Northern Polar Stereographic
#

if True:
    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig)
    fig = plt.figure(ifig)
    iplot += 1

    sub = fig.add_axes(jams.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))
    m   = Basemap(projection='npstere', boundinglat=30, lon_0=0, resolution='c')

    clab = r'$\mathrm{Maximum \; Snow \; Depth \; [m]}$'

    xx, yy = m(foster_lonh,foster_lath)
    zz     = foster_snow_max

    ncolor  = 8  # of colors in plot
    b1 = (1,1,1) #
    b2 = [ i/255. for i in (4, 90, 141) ]
    tmp  = jams.color.rgb_range(b1, b2, ncolor, cmap='MyBlue')
    cmap = mpl.cm.get_cmap('MyBlue')
    #cmap = jams.color.get_brewer('RdBu'+str(ncolor),reverse=True)

    mini  = 0.
    maxi  = 1.5
    cvals = mini + np.arange(ncolor+1)/np.float(ncolor)*(maxi-mini)

    norm = mpl.colors.BoundaryNorm(cvals, cmap.N)
    c    = m.pcolor(xx, yy, zz, norm=norm, cmap=cmap)
    cbar = plt.colorbar(c, orientation='vertical', ax=sub, shrink=0.5, pad=0.1)

    cnames = [r"${0:.1f}$".format(i) for i in cvals]
    cbar.set_ticks(cvals)
    cbar.set_ticklabels(cnames)
    for t in cbar.ax.get_yticklabels(): t.set_fontsize(cbsize)

    cbar.set_label(clab, fontsize=cbsize)
    # set continents
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries()
    #m.fillcontinents(color='coral',lake_color='aqua')
    #m.drawmapboundary(fill_color='0.9')
    # set grid
    parallels = np.arange(-80, 81, delat)
    m.drawparallels(parallels, labels=[1,1,0,0], fontsize=xsize, linewidth=alwidth)
    meridians = np.arange(-180, 181, delon)
    m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=xsize, linewidth=alwidth)

    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = pngbase+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)


# -------------------------------------------------------------------------
# Fig 2 - Germany
#

if True:
    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig)
    fig = plt.figure(ifig)

    # basemap lower left and upper right corners
    lllon = 4.0
    lllat = 47.0
    urlon = 16.0
    urlat = 56.0
    delon = 2
    delat = 1

    iplot += 1

    sub = fig.add_axes(jams.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))
    m    = (Basemap(projection='merc', resolution='i',
                    llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat))

    clab = r'$\mathrm{Maximum \; Snow \; Depth \; [m]}$'

    xx, yy = m(foster_lonh,foster_lath)
    zz     = foster_snow_max

    ncolor  = 8  # of colors in plot
    b1 = (1,1,1) #
    b2 = [ i/255. for i in (4, 90, 141) ]
    tmp  = jams.color.rgb_range(b1, b2, ncolor, cmap='MyBlue')
    cmap = mpl.cm.get_cmap('MyBlue')
    #cmap = jams.color.get_brewer('RdBu'+str(ncolor),reverse=True)

    mini  = 0.
    maxi  = 1.5
    cvals = mini + np.arange(ncolor+1)/np.float(ncolor)*(maxi-mini)

    norm = mpl.colors.BoundaryNorm(cvals, cmap.N)
    c    = m.pcolor(xx, yy, zz, norm=norm, cmap=cmap)
    cbar = plt.colorbar(c, orientation='vertical', ax=sub, shrink=0.5, pad=0.1)

    cnames = [r"${0:.1f}$".format(i) for i in cvals]
    cbar.set_ticks(cvals)
    cbar.set_ticklabels(cnames)
    for t in cbar.ax.get_yticklabels(): t.set_fontsize(cbsize)

    cbar.set_label(clab, fontsize=cbsize)
    # set continents
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries()
    #m.fillcontinents(color='coral',lake_color='aqua')
    #m.drawmapboundary(fill_color='0.9')
    # set grid
    parallels = np.arange(-80, 81, delat)
    m.drawparallels(parallels, labels=[1,1,0,0], fontsize=xsize, linewidth=0)#, linewidth=alwidth)
    meridians = np.arange(-180, 181, delon)
    m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=xsize, linewidth=0)#, linewidth=alwidth)

    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = pngbase+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)


# -------------------------------------------------------------------------
# Fig 3 - Globe
#

if True:
    ifig += 1
    iplot = 0
    print('Plot - Fig ', ifig)
    fig = plt.figure(ifig)

    # basemap lower left and upper right corners
    delon = 60
    delat = 30

    iplot += 1

    sub = fig.add_axes(jams.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))
    m    = Basemap(projection='robin', lon_0=0, resolution='c')

    clab = r'$\mathrm{Maximum \; Snow \; Depth \; [m]}$'

    xx, yy = m(foster_lonh,foster_lath)
    zz     = foster_snow_max

    ncolor  = 8  # of colors in plot
    b1 = (1,1,1) #
    b2 = [ i/255. for i in (4, 90, 141) ]
    tmp  = jams.color.rgb_range(b1, b2, ncolor, cmap='MyBlue')
    cmap = mpl.cm.get_cmap('MyBlue')
    #cmap = jams.color.get_brewer('RdBu'+str(ncolor),reverse=True)

    mini  = 0.
    maxi  = 1.5
    cvals = mini + np.arange(ncolor+1)/np.float(ncolor)*(maxi-mini)

    norm = mpl.colors.BoundaryNorm(cvals, cmap.N)
    c    = m.pcolor(xx, yy, zz, norm=norm, cmap=cmap)
    cbar = plt.colorbar(c, orientation='horizontal', ax=sub, shrink=0.5, pad=0.05)

    cnames = [r"${0:.1f}$".format(i) for i in cvals]
    cbar.set_ticks(cvals)
    cbar.set_ticklabels(cnames)
    for t in cbar.ax.get_yticklabels(): t.set_fontsize(cbsize)

    cbar.set_label(clab, fontsize=cbsize)
    # set continents
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries()
    #m.fillcontinents(color='coral',lake_color='aqua')
    #m.drawmapboundary(fill_color='0.9')
    # set grid
    parallels = np.arange(-80, 81, delat)
    m.drawparallels(parallels, labels=[1,1,0,0], fontsize=xsize, linewidth=alwidth)
    meridians = np.arange(-180, 181, delon)
    m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=xsize, linewidth=alwidth)

    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = pngbase+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)


# -------------------------------------------------------------------------
# Finished
#

if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()

t2    = time.time()
strin = '[m]: '+jams.astr((t2-t1)/60.,1) if (t2-t1)>60. else '[s]: '+jams.astr(t2-t1,0)
print('Time ', strin)
