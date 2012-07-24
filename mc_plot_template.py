#!/usr/bin/env python2.7
"""
The following information is the output of: python mc_plot_template.py -h


Usage: mc_plot_template.py [options]

This is the python plot template of Matthias Cuntz.

Options:
  -h, --help            show this help message and exit
  -p File, --pdffile=File
                        Name of pdf output file (default: open X-window).
  -t, --usetex          Use LaTeX to render text in pdf.



License
-------
This file is part of the UFZ Python library.

The UFZ Python library is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The UFZ Python library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with The UFZ Python library.  If not, see <http://www.gnu.org/licenses/>.

Copyright 2012 Matthias Cuntz


History
-------
Written, MC, Jul 2012
"""
import numpy as np
import ufz

# -------------------------------------------------------------------------
# Command line arguments
#

pdffile    = ''
usetex     = False
import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="This is the python plot template of Matthias Cuntz.")
parser.add_option('-p', '--pdffile', action='store', dest='pdffile', type='string',
                  default=pdffile, metavar='File',
                  help='Name of pdf output file (default: open X-window).')
parser.add_option('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                  help="Use LaTeX to render text in pdf.")
(opts, args) = parser.parse_args()

pdffile  = opts.pdffile
usetex   = opts.usetex
del parser, opts, args


# -------------------------------------------------------------------------
# Customize plots
#

if (pdffile == ''):
    outtype = 'x'
else:
    outtype = 'pdf'

# Main plot
nrow       = 3           # # of rows of subplots per figure
ncol       = 2           # # of columns of subplots per figure
hspace     = 0.10        # x-space between subplots
wspace     = 0.05        # y-space between subplots
textsize   = 13          # standard text size
dxabc      = 0.90        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc      = 0.05        # % of (max-min) shift up from lower x-axis for a,b,c,... labels

lwidth     = 1.5         # linewidth
elwidth    = 1.0         # errorbar line width
alwidth    = 1.0         # axis line width
msize      = 1.0         # marker size
mwidth     = 1.0         # marker edge width
# color: 'b'|'g'|'r'|'c'|'m'|'y'|'k'|'w'
#        'blue'|'green'|'red'|'cyan'|'magenta'|'yellow'|'black'|'white'
#        hex string '#eeefff' | RGB tuple (1,0.5,1) | html names 'burlywod', 'chartreuse', ...
#        grayscale intensity, e.g. '0.7', 'k'='0.0'
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
    mpl.rc('text.latex', unicode=True)
    #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    mpl.rc('font',**{'family':'serif','serif':['times']})
  mpl.rc('font', size=textsize)
else:
  import matplotlib.pyplot as plt
  mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
  mpl.rc('font', size=textsize)
mpl.rc('lines', linewidth=lwidth, color='black')
mpl.rc('axes', linewidth=alwidth, labelcolor='black')
mpl.rc('path', simplify=False) # do not remove


# -------------------------------------------------------------------------
# Prepare some psudo data
#

print 'Prepare data data'

ndata = 1000
mu, sigma = 0, 0.1 # mean and standard deviation
yy = np.random.normal(mu, sigma, ndata)

xx = np.random.normal(0, 1, ndata)

# -------------------------------------------------------------------------
# Plot
#

if (outtype == 'pdf'):
    print 'Plot PDF ', pdffile
    pdf_pages = PdfPages(pdffile)
else:
    print 'Plot X'
figsize = mpl.rcParams['figure.figsize']

ifig = 0


# -------------------------------------------------------------------------
# Fig 1 - xx vs. yy
#

ifig += 1
iplot = 0
print 'Plot - Fig ', ifig
fig = plt.figure(ifig)

# Set to None for free scaling at first
xlim = None
ylim = None
# xlim = [0.,1.]
# ylim = [-0.4,0.4]
xtick = 0.2

iplot += 1
xxplot = xx
xlab   = ur'$(0,1)$'
yyplot = yy
ylab   = ur'$N('+ufz.astr(mu)+','+ufz.astr(sigma,1)+')$'
sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,wspace=wspace))
mark1  = sub.plot(xxplot, yyplot)
plt.setp(mark1, linestyle='None', marker='o', markeredgecolor=mcol1, markerfacecolor='None',
         markersize=msize, markeredgewidth=mwidth, label=ur'UniNorm')
sub.grid(False)
if xlim != None:
    plt.setp(sub, xlabel=xlab, xlim=lim)
else:
    plt.setp(sub, xlabel=xlab)
if ylim != None:
    plt.setp(sub, ylabel=ylab, ylim=lim)
else:
    plt.setp(sub, ylabel=ylab)
sub.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xtick))
ll = sub.legend(mark1, [ur'Uniform vs. Gauss'], frameon=frameon, ncol=1,
                labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                loc='upper right', bbox_to_anchor=(1+llxbbox,1+llybbox), scatterpoints=1, numpoints=1)
ufz.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True)

iplot += 1
xlab   = ur'$N('+ufz.astr(mu)+','+ufz.astr(sigma,1)+')$'
yyplot = yy
if usetex:
    ylab   = ur'Normed \#/bin'
else:
    ylab   = ur'Normed #/bin'
sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,wspace=wspace))
nbins  = ndata/30
#   inbins, ibins, ipatches = sub.hist(msiboot[i,:], bins=nbins)
count, bins, ipatches = sub.hist(yyplot, bins=nbins, normed=True)
for k in ipatches: plt.setp(k, facecolor=mcol1, edgecolor='none', linewidth=0.)
sub.grid(False)
gauss = 1./(sigma * np.sqrt(2.*np.pi)) * np.exp( -(bins - mu)**2 / (2 * sigma**2))
line1 = sub.plot(bins, gauss)
plt.setp(line1, linestyle='-', linewidth=lwidth, marker=None, color=lcol1, label='Gauss')
sub.grid(False)
if xlim != None:
    plt.setp(sub, xlabel=xlab, xlim=lim)
else:
    plt.setp(sub, xlabel=xlab)
if ylim != None:
    plt.setp(sub, ylabel=ylab, ylim=lim)
else:
    plt.setp(sub, ylabel=ylab)
ufz.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True)



if (outtype == 'pdf'):
  pdf_pages.savefig(fig)
  plt.close()


  
if (outtype == 'pdf'):
  pdf_pages.close()
else:
  plt.show()
