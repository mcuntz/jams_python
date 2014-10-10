#!/usr/bin/env python
from __future__ import print_function
"""
usage: mc_plot_template.py [-h] [-f flag] [-g pngbase] [-p pdffile] [-t]
                           [additional_args [additional_args ...]]

description:
  This is the python plot template of Matthias Cuntz.
  It provides quick view on screen or PDF or PNG output.
  It is designed for publication-ready PDF output so screen output is not optimal.

  The template also serves as a lookup of Python tips and tricks as well as a gallery
  of plots, graphs, maps.

positional arguments:
  additional_args       Any additional arguments such as filenames, etc.

optional arguments:
  -h, --help            show this help message and exit
  -f flag, --flag flag  Bit flag (default: 3): 1: do this; 2: do that.
  -g pngbase, --pngbase pngbase
                        Name basis for png output files (default: open screen
                        window).
  -p pdffile, --pdffile pdffile
                        Name of pdf output file (default: open screen window).
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
along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
If not, see <http://www.gnu.org/licenses/>.

Copyright 2012-2013 Matthias Cuntz


History
-------
Written,  MC, Jul 2012
Modified, MC, Jul 2013 - optparse->argparse
          MC, Jul 2013 - extended to be lookup and gallery
          MC, Dec 2013 - add png support
"""

# Hardcoded switches for plot gallery
dorandom  = True  # Scatter and histogram
dobasemap = True  # map with basemap
docartopy = True  # map with cartopy

# # python debugger
# import pdb
# pdb.set_trace() # set breakpoint

# -------------------------------------------------------------------------
# Command line arguments
#

import argparse
import textwrap

addargs = []
doflag  = 3
pngbase = ''
pdffile = ''
usetex  = False
parser  = argparse.ArgumentParser(
                                  formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=textwrap.dedent('''\
          description:
            This is the python plot template of Matthias Cuntz.
            It provides quick view on screen or PDF or PNG output.
            It is designed for publication-ready PDF output so screen output is not optimal.

            The template also serves as a lookup of Python tips and tricks as well as a gallery
            of plots, graphs, maps.
          '''))
parser.add_argument('-f', '--flag', action='store',
                    default=doflag, dest='doflag', metavar='flag', type=int,
                    help="Bit flag (default: 3):"
                         " 1: do this;"
                         " 2: do that.")
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
addargs = args.further
doflag  = args.doflag
pdffile = args.pdffile
pngbase = args.pngbase
usetex  = args.usetex

if (pdffile != '') & (pngbase != ''):
    print('\nError: PDF and PNG are mutually exclusive. Only either -p or -g possible.\n')
    parser.print_usage()
    import sys
    sys.exit()

del parser, args

# import packages after help so that help with command line -h is fast
import numpy as np
import ufz
import scipy.optimize as opt
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
mcol1       = ufz.colours('red')        # primary marker colour
mcol2       = '0.0'                     # secondary
mcol3       = (202/255.,0/255.,32/255.) # third
mcols       = ufz.colours(['blue','red','darkgray','orange','darkblue','black'])
lcol1       = ufz.colours('blue')   # primary line colour
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

if dobasemap:
    try:
        from mpl_toolkits.basemap import Basemap, shiftgrid
    except:
        dobasemap = False
        
if docartopy:
    try:
        import cartopy.crs as ccrs
    except:
        docartopy = False


# -------------------------------------------------------------------------
# Analyse flag
#

# -2 is -number of possible flags, there should be the same amount of zeros in the '00' string
dothis = False
dothat = False
bflag = (('00' + np.binary_repr(doflag))[-2:])[::-1] # revert for easy extension
if bflag[0]=='1': dothis = True
if bflag[1]=='1': dothat = True


# -------------------------------------------------------------------------
# Prepare some pseudo data
#

import os
print('Prepare data')

if dorandom:
  # Random data
  ndata     = 1000
  mu, sigma = 0, 0.1 # mean, std
  yy        = np.random.normal(mu, sigma, ndata)
  xx        = np.random.rand(ndata)

  # Pseudo temperature
  lon        =   0.5 + np.arange(360, dtype=np.float)
  lat        = -89.5 + np.arange(180, dtype=np.float)
  lon2, lat2 = np.meshgrid(lon,lat)
  temp       = -263. + np.cos(np.deg2rad(lat2)) + np.sin(np.deg2rad(lon2))

  # Write ascii data to file with print
  ofile='mc_plot_template.csv'
  print('Write ascii data ', ofile)
  f = open(ofile,'w')
  print('uniform gauss', file=f)
  for i in range(xx.shape[0]):
      print(ufz.astr(xx[i],9), ufz.astr(yy[i],9), file=f)
  f.close()
  del xx, yy

  # Read in ascii data and header. Define variables with names of header
  ifile = ofile
  print('Read ascii data ', ifile)
  dat   = ufz.fread(ifile, skip=1, transpose=True)
  head  = ufz.fread(ifile, skip=1, header=True)
  ihead = dict(zip(head, range(len(head))))
  for ii in range(len(head)): exec(head[ii] + ' = ' + 'dat[ihead["'+head[ii]+'"],:]')
  del dat
  os.remove(ifile)

  # Write to numpy zipped file
  ofile   ='mc_plot_template.npz'
  allvars = ['uniform', 'gauss']
  print('Write npz-file', ofile)
  p = open(ofile, 'wb')
  savearg = 'p'
  for j in allvars: savearg = savearg + ', '+j+'='+j
  exec("np.savez_compressed("+savearg+")")
  p.close()

  # Read numpy zipped file
  ifile = ofile
  print('Read npz-file', ifile)
  p = open(ifile, 'rb')
  pp = np.load(p)
  for i in pp.files: exec(i+" = pp['"+i+"']")
  p.close()
  os.remove(ifile)

  # Fit polynomial to Random Gauss
  print('Fit Gauss')
  nbins  = ndata/30
  hist, bin_edges = np.histogram(gauss, bins=nbins, density=True)
  bin_sizes   = np.diff(bin_edges)
  bin_centres = bin_edges[0]-0.5*bin_sizes[0] + np.cumsum(bin_sizes)
  gauss_p = opt.fmin(ufz.functions.cost_gauss, np.array([1,1]), args=(bin_centres,hist), disp=False)

if dobasemap | docartopy:
  # Pseudo temperature
  nlon       = 360//5
  nlat       = 180//5
  lon        = 360./(2.*nlon)        + np.arange(nlon)/np.float(nlon) * 360.
  lat        = -90. + 180./(2.*nlat) + np.arange(nlat)/np.float(nlat) * 180.
  lon2, lat2 = np.meshgrid(lon,lat)
  temp       = (np.cos(np.deg2rad(lat2)) + np.sin(np.deg2rad(lon2*10)))*20.

  # grid cell edges
  lonh = np.empty((nlat+1,nlon+1), dtype=np.float)
  lath = np.empty((nlat+1,nlon+1), dtype=np.float)
  dlon = np.diff(lon2,axis=1)
  dlat = np.diff(lat2,axis=0)
  lonh[0:nlat,0]      = lon2[:,0]        - 0.5*dlon[:,0]          # 0
  lonh[0:nlat,1:nlon] = lonh[0:nlat,0:1] + np.cumsum(dlon,axis=1) # 1:nlon+1
  lonh[0:nlat,-1]     = lon2[:,-1]       + 0.5*dlon[:,-1]         # nlon
  lonh[nlat,0:nlon+1] = lonh[nlat-1,0:nlon+1]                    # lower corner of box: assign last upper boundary
  lath[0,0:nlon]      = lat2[0,:]        - 0.5*dlat[0,:]
  lath[1:nlat,0:nlon] = lath[0:1,0:nlon] + np.cumsum(dlat,axis=0)
  lath[-1,0:nlon]     = lat2[-1,:]       + 0.5*dlat[-1,:]
  lath[0:nlat+1,nlon] = lath[0:nlat+1,nlon-1]

  slon         = np.roll(lon,nlon//2) # np.where(lon>180, lon-360., lon)
  slat         = lat
  stemp        = np.roll(temp,nlon//2,1)
  slon2, slat2 = np.meshgrid(slon,slat)
  slonh = np.empty((nlat+1,nlon+1), dtype=np.float)
  slath = np.empty((nlat+1,nlon+1), dtype=np.float)
  dlon = np.diff(slon2,axis=1)
  dlat = np.diff(slat2,axis=0)
  slonh[0:nlat,0]      = slon2[:,0]        - 0.5*dlon[:,0]          # 0
  slonh[0:nlat,1:nlon] = slonh[0:nlat,0:1] + np.cumsum(dlon,axis=1) # 1:nlon+1
  slonh[0:nlat,-1]     = slon2[:,-1]       + 0.5*dlon[:,-1]         # nlon
  slonh[nlat,0:nlon+1] = slonh[nlat-1,0:nlon+1]                    # lower corner of box: assign last upper boundary
  slath[0,0:nlon]      = slat2[0,:]        - 0.5*dlat[0,:]
  slath[1:nlat,0:nlon] = slath[0:1,0:nlon] + np.cumsum(dlat,axis=0)
  slath[-1,0:nlon]     = slat2[-1,:]       + 0.5*dlat[-1,:]
  slath[0:nlat+1,nlon] = slath[0:nlat+1,nlon-1]

  mlon         = np.roll(np.where(lon>180, lon-360., lon),nlon//2)
  mlat         = lat
  mtemp        = np.roll(temp,nlon//2,1)
  mlon2, mlat2 = np.meshgrid(mlon,mlat)
  mlonh = np.empty((nlat+1,nlon+1), dtype=np.float)
  mlath = np.empty((nlat+1,nlon+1), dtype=np.float)
  dlon = np.diff(mlon2,axis=1)
  dlat = np.diff(mlat2,axis=0)
  mlonh[0:nlat,0]      = mlon2[:,0]        - 0.5*dlon[:,0]          # 0
  mlonh[0:nlat,1:nlon] = mlonh[0:nlat,0:1] + np.cumsum(dlon,axis=1) # 1:nlon+1
  mlonh[0:nlat,-1]     = mlon2[:,-1]       + 0.5*dlon[:,-1]         # nlon
  mlonh[nlat,0:nlon+1] = mlonh[nlat-1,0:nlon+1]                    # lower corner of box: assign last upper boundary
  mlath[0,0:nlon]      = mlat2[0,:]        - 0.5*dlat[0,:]
  mlath[1:nlat,0:nlon] = mlath[0:1,0:nlon] + np.cumsum(dlat,axis=0)
  mlath[-1,0:nlon]     = mlat2[-1,:]       + 0.5*dlat[-1,:]
  mlath[0:nlat+1,nlon] = mlath[0:nlat+1,nlon-1]

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

if dorandom:
  # -------------------------------------------------------------------------
  # Fig 1 - gauss vs. uniform
  #

  ifig += 1
  iplot = 0
  print('Plot - Fig ', ifig)
  fig = plt.figure(ifig)

  # Set to None for free scaling at first, then set limits
  xlim = None
  ylim = None
  # xlim = [0.,1.]
  # ylim = [-0.4,0.4]
  xtick = 0.2

  # Scatter
  iplot += 1
  xlab   = r'$(0,1)$'
  ylab   = r'$\aleph(\mu = '+ufz.astr(mu)+', \sigma = '+ufz.astr(sigma,1)+')$'
  sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

  # possible parameters
  # linestyle: ['-'|'--'|'-.'|':'|'None'|' '|'' ] and any drawstyle in combination with a linestyle, e.g. 'steps--'.
  # marker         description
  # '.'            point
  # ','            pixel
  # 'o'            circle
  # 'v'            triangle_down
  # '^'            triangle_up
  # '<'            triangle_left
  # '>'            triangle_right
  # '1'            tri_down
  # '2'            tri_up
  # '3'            tri_left
  # '4'            tri_right
  # 's'            square
  # 'p'            pentagon
  # '*'            star
  # 'h'            hexagon1
  # 'H'            hexagon2
  # '+'            plus
  # 'x'            x
  # 'D'            diamond
  # 'd'            thin_diamond
  # '|'            vline
  # '_'            hline
  # TICKLEFT       tickleft
  # TICKRIGHT      tickright
  # TICKUP         tickup
  # TICKDOWN       tickdown
  # CARETLEFT      caretleft
  # CARETRIGHT     caretright
  # CARETUP        caretup
  # CARETDOWN      caretdown
  # 'None'         nothing
  # ' '            nothing
  # ''             nothing
  # color: 'b'|'g'|'r'|'c'|'m'|'y'|'k'|'w'
  #        'blue'|'green'|'red'|'cyan'|'magenta'|'yellow'|'black'|'white'
  #        hex string '#eeefff' | RGB tuple (1,0.5,1) | html names 'burlywod', 'chartreuse', ...
  #        grayscale intensity '0.7'
  mark1  = sub.plot(uniform, gauss) # marker
  plt.setp(mark1, linestyle='None', marker='o', markeredgecolor=mcol1, markerfacecolor='None',
           markersize=msize, markeredgewidth=mwidth, label=r'UniNorm')

  plt.setp(sub, xlabel=xlab) # axis labels
  plt.setp(sub, ylabel=ylab)
  sub.grid(False)

  sub.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xtick)) # tick locations

  if xlim != None: plt.setp(sub, xlim=xlim) # set axis limit if wanted
  if ylim != None: plt.setp(sub, ylim=ylim)

  ll = sub.legend(mark1, [r'$\mathrm{\ddot{U}niform \; vs. \; G\ddot{a}uss}$'], frameon=frameon, ncol=1,
                  labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                  loc='upper right', bbox_to_anchor=(1+llxbbox,1+llybbox), scatterpoints=1, numpoints=1)
  plt.setp(ll.get_texts(), fontsize='small')

  ufz.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=usetex, mathrm=True, parenthesis='close')

  # redundant, only for demonstration
  xmin, xmax = sub.get_xlim()
  ymin, ymax = sub.get_ylim()
  plt.setp(sub, xlim=[xmin,xmax], ylim=[ymin,ymax])


  # Gauss curve
  iplot += 1
  xlab   = r'$N('+ufz.astr(mu)+','+ufz.astr(sigma,1)+')$'
  ylab   = r'$\mathrm{Rel. \; frequency}$'
  sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

  nbins  = ndata/30
  count, bins, ipatches = sub.hist(gauss, bins=nbins, normed=True)
  for k in ipatches: plt.setp(k, facecolor=mcol1, edgecolor='none', linewidth=0.)

  line1  = sub.plot(bins, ufz.functions.gauss(bins,mu,sigma))
  plt.setp(line1, linestyle='-', linewidth=lwidth, marker=None, color=lcol1, label='Gauss')

  line2  = sub.plot(bin_centres, ufz.functions.gauss_p(bin_centres,gauss_p))
  plt.setp(line2, linestyle='--', linewidth=lwidth, marker=None, color=lcols[4], label='Fit')

  plt.setp(sub, xlabel=xlab) # axis labels
  plt.setp(sub, ylabel=ylab)
  sub.grid(False)

  if xlim != None: plt.setp(sub, xlim=xlim) # set axis limit if wanted
  if ylim != None: plt.setp(sub, ylim=ylim)

  ll = sub.legend(line1+line2, [r'$\mathrm{Theoretical}$',r'$\mathrm{Fitted}$'], frameon=frameon, ncol=1,
                  labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                  loc='upper right', bbox_to_anchor=(1+llxbbox,1+llybbox), scatterpoints=1, numpoints=1)
  plt.setp(ll.get_texts(), fontsize='small')

  ufz.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=usetex, mathrm=True, parenthesis='close')


  # Same as scatter but set xlabel names
  iplot += 1
  xlab   = r'$(0,1)$'
  ylab   = r'$\aleph(\mu = '+ufz.astr(mu)+', \sigma = '+ufz.astr(sigma,1)+')$'
  sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

  mark1  = sub.plot(uniform, gauss) # marker
  plt.setp(mark1, linestyle='None', marker='o', markeredgecolor=mcol1, markerfacecolor='None',
           markersize=msize, markeredgewidth=mwidth, label=r'UniNorm')

  plt.setp(sub, xlabel=xlab) # axis labels
  plt.setp(sub, ylabel=ylab)
  sub.grid(False)

  sub.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xtick)) # tick locations

  if xlim != None: plt.setp(sub, xlim=xlim) # set axis limit if wanted
  if ylim != None: plt.setp(sub, ylim=ylim)

  ll = sub.legend(mark1, [r'$\mathrm{Uniform \; vs. \; Gauss}$'], frameon=frameon, ncol=1,
                  labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                  loc='upper right', bbox_to_anchor=(1+llxbbox,1+llybbox), scatterpoints=1, numpoints=1)
  plt.setp(ll.get_texts(), fontsize='small')

  ufz.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=usetex, mathrm=True, parenthesis='close')

  ixticks = sub.get_xticks()
  xnames  = [ r'$\mathrm{U'+str(i)+'}$' for i in ixticks]
  #plt.setp(sub, xticks=ixticks)
  xticknames = plt.setp(sub, xticklabels=xnames)
  plt.setp(xticknames, rotation=45, fontsize='small')


  # Same as scatter but with second axis to the right
  iplot += 1
  xlab   = r'$(0,1)$'
  ylab   = r'$\aleph(\mu = '+ufz.astr(mu)+', \sigma = '+ufz.astr(sigma,1)+')$'
  sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

  mark1  = sub.plot(uniform, gauss) # marker
  plt.setp(mark1, linestyle='None', marker='o', markeredgecolor=mcol1, markerfacecolor='None',
           markersize=msize, markeredgewidth=mwidth, label=r'UniNorm')

  plt.setp(sub, xlabel=xlab) # axis labels
  plt.setp(sub, ylabel=ylab)
  sub.grid(False)

  sub.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xtick)) # tick locations

  if xlim != None: plt.setp(sub, xlim=xlim) # set axis limit if wanted
  if ylim != None: plt.setp(sub, ylim=ylim)

  ll = sub.legend(mark1, [r'$\mathrm{Uniform \; vs. \; Gauss}$'], frameon=frameon, ncol=1,
                  labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                  loc='upper right', bbox_to_anchor=(1+llxbbox,1+llybbox), scatterpoints=1, numpoints=1)
  plt.setp(ll.get_texts(), fontsize='small')

  ufz.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=usetex, mathrm=True, parenthesis='close')

  ylab2 = r'$\aleph(\mu = '+ufz.astr(mu)+', \sigma = '+ufz.astr(sigma/2.,2)+')$'
  sub2  = sub.twinx()
  mark2 = sub2.plot(uniform, gauss*0.5)
  plt.setp(mark2, linestyle='None', marker='^', markeredgecolor=mcol2, markerfacecolor='None',
           markersize=msize, markeredgewidth=mwidth, label=r'UniNorm2')
  plt.setp(sub2, ylabel=ylab2)


  # Another scatter plot with the same y-axis scale as right axis before -> re-sets the axis limits above
  # Give points different colors
  iplot += 1
  xlab   = r'$(0,1)$'
  ylab   = r'$\aleph(\mu = '+ufz.astr(mu)+', \sigma = '+ufz.astr(sigma*2,1)+')$'
  sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace),sharey=sub2)

  cm = ufz.get_brewer('RdYlBu'+str(ncolor))
  ii = np.argsort(uniform)
  mark1  = sub.scatter(uniform[ii], gauss[ii]*2., marker='o', c=np.arange(uniform.size)%(uniform.size//3), cmap=cm,
                       linewidth=mwidth, s=msize) # colored markers in 3 repeating colour sequence
  plt.setp(mark1, edgecolor=mcol1, facecolor='None')

  plt.setp(sub, xlabel=xlab) # axis labels
  plt.setp(sub, ylabel=ylab)
  sub.grid(False)

  sub.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xtick)) # tick locations

  if xlim != None: plt.setp(sub, xlim=xlim) # set axis limit if wanted
  if ylim != None: plt.setp(sub, ylim=ylim)

  ufz.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=usetex, mathrm=True, parenthesis='close')


  if (outtype == 'pdf'):
      pdf_pages.savefig(fig)
      plt.close(fig)
  elif (outtype == 'png'):
      pngfile = pngbase+"{0:04d}".format(ifig)+".png"
      fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
      plt.close(fig)


if dobasemap:
  # -------------------------------------------------------------------------
  # Fig 2 - Maps with basemap
  #

  ifig += 1
  iplot = 0
  print('Plot - Fig ', ifig)
  fig = plt.figure(ifig)

  nrow = 3
  ncol = 1
  cm = ufz.get_brewer('RdYlBu'+str(ncolor),reverse=True)
  delon = 60
  delat = 30

  # Mollweide
  iplot += 1

  sub  = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))
  m    = Basemap(projection='moll',lon_0=0,resolution='c')

  mx, my = m(mlon2,mlat2)
  sx, sy = m(slon2,slat2)

  cvals  = np.amin(temp) + np.arange(ncolor+1)/np.float(ncolor)*(np.amax(temp)-np.amin(temp))
  norm   = mpl.colors.BoundaryNorm(cvals, cm.N)
  if nlon > 100:
      c = m.contourf(mx, my, temp, norm=norm, cmap=cm)
  else:
      c = m.pcolormesh(sx, sy, stemp, norm=norm, cmap=cm)
  cbar   = plt.colorbar(c, orientation='vertical')
  cnames = [r"${0:.0f}$".format(i) for i in cvals]
  cbar.set_ticks(cvals)
  cbar.set_ticklabels(cnames)
  for t in cbar.ax.get_yticklabels(): t.set_fontsize(cbsize)
  clab = r'$\mathrm{T_{surf} \; [^oC]}$'
  cbar.set_label(clab, fontsize=cbsize)
  # set continents
  m.drawcoastlines(linewidth=0.5)
  m.drawcountries()
  # set grid
  parallels = np.arange(-90., 90.+delat, delat)
  m.drawparallels(parallels, labels=[1,1,0,0], fontsize=xsize, linewidth=0)
  # meridians = np.arange(0., 360.+delon, delon)
  # m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=xsize, linewidth=0)

  # LatLon
  iplot += 1

  sub  = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))
  m    = Basemap(projection='cyl', lon_0=0, resolution='c')

  mx, my = m(mlon2,mlat2)

  cvals  = np.amin(temp) + np.arange(ncolor+1)/np.float(ncolor)*(np.amax(temp)-np.amin(temp))
  norm   = mpl.colors.BoundaryNorm(cvals, cm.N)
  if nlon > 100:
      c = m.contourf(mx, my, mtemp, norm=norm, cmap=cm)
  else:
      c = m.pcolormesh(mx, my, mtemp, norm=norm, cmap=cm)
  cbar   = plt.colorbar(c, orientation='vertical')
  cnames = [r"${0:.0f}$".format(i) for i in cvals]
  cbar.set_ticks(cvals)
  cbar.set_ticklabels(cnames)
  for t in cbar.ax.get_yticklabels(): t.set_fontsize(cbsize)
  clab = r'$\mathrm{T_{surf} \; [^oC]}$'
  cbar.set_label(clab, fontsize=cbsize)
  # set continents
  m.drawcoastlines(linewidth=0.5)
  m.drawcountries()
  # set grid
  parallels = np.arange(-90., 90.+delat, delat)
  m.drawparallels(parallels, labels=[1,0,0,0], fontsize=xsize, linewidth=0)
  meridians = np.arange(-180., 180.+delon, delon)
  m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=xsize, linewidth=0)

  # Germany
  inrow  = nrow
  incol  = ncol*2
  iplot = iplot*2 + 1

  xlim = [5.5,12.5]
  ylim = [47,55.5]
  delon = 2
  delat = 2

  sub  = fig.add_axes(ufz.position(inrow,incol,iplot,hspace=hspace,vspace=vspace))
  m    = Basemap(projection='cyl', lon_0=0, resolution='c',
                 llcrnrlon=xlim[0], llcrnrlat=ylim[0], urcrnrlon=xlim[1], urcrnrlat=ylim[1])

  mx, my = m(mlon2,mlat2)

  cvals  = np.amin(temp) + np.arange(ncolor+1)/np.float(ncolor)*(np.amax(temp)-np.amin(temp))
  norm   = mpl.colors.BoundaryNorm(cvals, cm.N)
  if nlon > 100:
      c = m.contourf(mx, my, mtemp, norm=norm, cmap=cm)
  else:
      c = m.pcolormesh(mx, my, mtemp, norm=norm, cmap=cm)
  cbar   = plt.colorbar(c, orientation='vertical')
  cnames = [r"${0:.0f}$".format(i) for i in cvals]
  cbar.set_ticks(cvals)
  cbar.set_ticklabels(cnames)
  for t in cbar.ax.get_yticklabels(): t.set_fontsize(cbsize)
  clab = r'$\mathrm{T_{surf} \; [^oC]}$'
  cbar.set_label(clab, fontsize=cbsize)
  # set continents
  m.drawcoastlines(linewidth=0.5)
  m.drawcountries()
  # set grid
  parallels = np.arange(-90., 90.+delat, delat)
  m.drawparallels(parallels, labels=[1,0,0,0], fontsize=xsize, linewidth=0)
  meridians = np.arange(-180., 180.+delon, delon)
  m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=xsize, linewidth=0)

  if (outtype == 'pdf'):
      pdf_pages.savefig(fig)
      plt.close(fig)
  elif (outtype == 'png'):
      pngfile = pngbase+"{0:04d}".format(ifig)+".png"
      fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
      plt.close(fig)


if docartopy:
  # -------------------------------------------------------------------------
  # Fig 2 - Maps with basemap
  #

  ifig += 1
  iplot = 0
  print('Plot - Fig ', ifig)
  fig = plt.figure(ifig)

  nrow = 3
  ncol = 1
  cm = ufz.get_brewer('RdYlBu'+str(ncolor),reverse=True)
  delon = 60
  delat = 60

  # Mollweide projection - no axis labels possible yet
  print('Plot Mollweide')
  iplot += 1

  sub  = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace),
                      projection=ccrs.Mollweide(central_longitude=0))

  cvals  = np.amin(temp) + np.arange(ncolor+1)/np.float(ncolor)*(np.amax(temp)-np.amin(temp))
  norm   = mpl.colors.BoundaryNorm(cvals, cm.N)
  if nlon > 100:
      c = sub.contourf(lon2, lat2, temp,
                       norm=norm, cmap=cm, transform=ccrs.PlateCarree())
  else:
      c = sub.pcolormesh(np.roll(lon2,nlon//2,1), np.roll(lat2,nlon//2,1), np.roll(temp,nlon//2,1),
                         norm=norm, cmap=cm, transform=ccrs.PlateCarree())
  cbar   = plt.colorbar(c, orientation='vertical')
  cnames = [r"${0:.0f}$".format(i) for i in cvals]
  cbar.set_ticks(cvals)
  cbar.set_ticklabels(cnames)
  for t in cbar.ax.get_yticklabels(): t.set_fontsize(cbsize)
  clab = r'$\mathrm{T_{surf} \; [^oC]}$'
  cbar.set_label(clab, fontsize=cbsize)

  sub.coastlines(resolution='110m')
  sub.gridlines(draw_labels=False)
  sub.set_global()


  # Lat/Lon projection with labels etc.
  print('Plot LatLon')
  iplot += 1

  sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace),
                        projection=ccrs.PlateCarree(central_longitude=0))

  cvals  = np.amin(temp) + np.arange(ncolor+1)/np.float(ncolor)*(np.amax(temp)-np.amin(temp))
  norm   = mpl.colors.BoundaryNorm(cvals, cm.N)
  if nlon > 100:
      c = sub.contourf(lon2, lat2, temp,
                       norm=norm, cmap=cm, transform=ccrs.PlateCarree())
  else:
      c = sub.pcolormesh(np.roll(lon2,nlon//2,1), np.roll(lat2,nlon//2,1), np.roll(temp,nlon//2,1),
                         norm=norm, cmap=cm, transform=ccrs.PlateCarree())
  cbar   = plt.colorbar(c, orientation='vertical')
  cnames = [r"${0:.0f}$".format(i) for i in cvals]
  cbar.set_ticks(cvals)
  cbar.set_ticklabels(cnames)
  for t in cbar.ax.get_yticklabels(): t.set_fontsize(cbsize)
  clab = r'$\mathrm{T_{surf} \; [^oC]}$'
  cbar.set_label(clab, fontsize=cbsize)

  sub.coastlines(resolution='110m')
  sub.gridlines(draw_labels=False) # DIY below
  sub.set_global()

  sub.set_xticks(np.arange(-180., 180.+delon, delon))
  sub.set_yticks(np.arange(-90., 90.+delat, delat))
  ixticks = sub.get_xticks()
  xnames  = [ r'$\mathrm{'+str(int(i))+'\,^oN}$' for i in ixticks]
  xticknames = plt.setp(sub, xticklabels=xnames)
  plt.setp(xticknames, fontsize='medium')
  iyticks = sub.get_yticks()
  ynames  = [ r'$\mathrm{'+str(int(i))+'\,^oE}$' for i in iyticks]
  yticknames = plt.setp(sub, yticklabels=ynames)
  plt.setp(yticknames, fontsize='medium')



  # Zoom on Germany
  print('Plot Germany')
  inrow  = nrow
  incol  = ncol*2
  iplot = iplot*2 + 1

  xlim = [5.5,12.5]
  ylim = [47,55.5]
  delon = 2
  delat = 2

  sub    = fig.add_axes(ufz.position(inrow,incol,iplot,hspace=hspace,vspace=vspace),
                        projection=ccrs.PlateCarree())

  cvals  = np.amin(temp) + np.arange(ncolor+1)/np.float(ncolor)*(np.amax(temp)-np.amin(temp))
  norm   = mpl.colors.BoundaryNorm(cvals, cm.N)
  if nlon > 100:
      c = sub.contourf(lon2, lat2, temp,
                       norm=norm, cmap=cm, transform=ccrs.PlateCarree())
  else:
      c = sub.pcolormesh(np.roll(lon2,nlon//2,1), np.roll(lat2,nlon//2,1), np.roll(temp,nlon//2,1),
                         norm=norm, cmap=cm, transform=ccrs.PlateCarree())
  cbar   = plt.colorbar(c, orientation='vertical')
  cnames = [r"${0:.0f}$".format(i) for i in cvals]
  cbar.set_ticks(cvals)
  cbar.set_ticklabels(cnames)
  for t in cbar.ax.get_yticklabels(): t.set_fontsize(cbsize)
  clab = r'$\mathrm{T_{surf} \; [^oC]}$'
  cbar.set_label(clab, fontsize=cbsize)

  sub.coastlines(resolution='110m')
  sub.gridlines(draw_labels=False) # DIY below

  sub.set_xticks(np.arange(-180., 180.+delon, delon))
  sub.set_yticks(np.arange(-90., 90.+delat, delat))
  ixticks = sub.get_xticks()
  xnames  = [ r'$\mathrm{'+str(int(i))+'\,^oN}$' for i in ixticks]
  xticknames = plt.setp(sub, xticklabels=xnames)
  plt.setp(xticknames, fontsize='medium')
  iyticks = sub.get_yticks()
  ynames  = [ r'$\mathrm{'+str(int(i))+'\,^oE}$' for i in iyticks]
  yticknames = plt.setp(sub, yticklabels=ynames)
  plt.setp(yticknames, fontsize='medium')

  sub.set_xlim(xlim)
  sub.set_ylim(ylim)


  if (outtype == 'pdf'):
      pdf_pages.savefig(fig)
      plt.close(fig)
  elif (outtype == 'png'):
      pngfile = pngbase+"{0:04d}".format(ifig)+".png"
      fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
      plt.close(fig)


if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()

t2    = time.time()
strin = '[m]: '+ufz.astr((t2-t1)/60.,1) if (t2-t1)>60. else '[s]: '+ufz.astr(t2-t1,0)
print('Time for demonstration', strin)
