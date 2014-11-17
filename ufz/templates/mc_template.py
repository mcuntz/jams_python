#!/usr/bin/env python
from __future__ import print_function
"""
usage: mc_template.py [-h] [-u] [-p plotname] [-t outtype] [infile]

This is the Python template for any new program of Matthias Cuntz.

positional arguments:
  infile                Mandatory input file.

optional arguments:
  -h, --help            show this help message and exit
  -p plotname, --plotname plotname
                        Name of plot output file for types pdf, html or d3,
                        and name basis for type png (default: mc_template).
  -t outtype, --type outtype
                        Output type is pdf, png, html, or d3 (default: open
                        screen windows).
  -u, --usetex          Use LaTeX to render text in pdf, png and html.


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

Copyright 2012-2014 Matthias Cuntz


History
-------
Written,  MC, Jul 2012
Modified, MC, Jul 2013 - optparse->argparse
          MC, Jul 2013 - extended to be lookup and gallery
          MC, Dec 2013 - add png support
          MC, Mar 2014 - split into individual templates
          MC, Nov 2014 - script -> function
                       - pdf, png, html, or d3
"""
if __name__ == '__main__':
    # -------------------------------------------------------------------------
    # Command line arguments
    #

    import argparse

    plotname = ''
    outtype  = ''
    usetex   = False
    parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description='''This is the Python template for any new program of Matthias Cuntz.''')
    parser.add_argument('-p', '--plotname', action='store',
                        default=plotname, dest='plotname', metavar='plotname',
                        help='Name of plot output file for types pdf, html or d3, '
                        'and name basis for type png (default: '+__file__[0:__file__.rfind(".")]+').')
    parser.add_argument('-t', '--type', action='store',
                        default=outtype, dest='outtype', metavar='outtype',
                        help='Output type is pdf, png, html, or d3 (default: open screen windows).')
    parser.add_argument('-u', '--usetex', action='store_true', default=usetex, dest="usetex",
                        help="Use LaTeX to render text in pdf, png and html.")
    parser.add_argument('file', nargs='?', default=None, metavar='infile',
                       help='Mandatory input file.')

    args     = parser.parse_args()
    infile   = args.file
    plotname = args.plotname
    outtype  = args.outtype
    usetex   = args.usetex

    del parser, args

# If function: uncomment section # FUNCTION and
#              move the if __name__ == '__main__': part above to end of file
# # FUNCTION
#     # Call routine
#     mc_template(infile, plotname=plotname, outtype=outtype, usetex=usetex)


# def mc_template(infile=None, plotname='', outtype='', usetex=False):
#     """
#     This is the Python template for any new program of Matthias Cuntz.

#     positional arguments:
#       infile                Mandatory input file.

#     optional arguments:
#       plotname              Name of plot output file if type is pdf, html or d3, and
#                             name basis if type is png (default: mc_template).
#       outtype               Output type is pdf, png, html, or d3 (default: open screen windows).
#       usetex                True: Use LaTeX to render text in pdf, png, and html (default: False)
#     """
# # FUNCTION
    # Check input
    if (infile is None):
        print('\nError: Input file must be given.\n')
        import sys
        sys.exit()

    outtype = outtype.lower()
    outtypes = ['', 'pdf', 'png', 'html', 'd3']
    if outtype not in outtypes:
        print('\nError: output type must be in ', outtypes)
        import sys
        sys.exit()

    import numpy as np
    import ufz
    import time
    t1 = time.time()

    if (outtype == 'd3'):
        try:
            import mpld3
        except:
            print("No mpld3 found. Use output type html instead of d3.")
            outtype = 'html'

    # -------------------------------------------------------------------------
    # Customize plots
    #

    # Main plot
    nrow        = 3           # # of rows of subplots per figure
    ncol        = 2           # # of columns of subplots per figure
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
            mpl.rc('text.latex', unicode=True)
            mpl.rcParams['text.latex.preamble']=r'\usepackage{wasysym}'
        else:
            #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
            mpl.rc('font',**{'family':'serif','serif':['times']})
        mpl.rc('text.latex', unicode=True)
    elif (outtype == 'png') or (outtype == 'html') or (outtype == 'd3'):
        mpl.use('Agg') # set directly after import matplotlib
        import matplotlib.pyplot as plt
        mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
        if usetex:
            mpl.rc('text', usetex=True)
            mpl.rc('text.latex', unicode=True)
            mpl.rcParams['text.latex.preamble']=r'\usepackage{wasysym}'
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


    # -------------------------------------------------------------------------
    # Plot
    #

    if (outtype == 'pdf'):
        if plotname == '':
            pdffile = __file__[0:__file__.rfind(".")] + '.pdf'
        else:
            pdffile = plotname
        print('Plot PDF ', pdffile)
        pdf_pages = PdfPages(pdffile)
    elif (outtype == 'png'):
        if plotname == '':
            pngbase = __file__[0:__file__.rfind(".")] + '_'
        else:
            pngbase = plotname
        print('Plot PNG ', pngbase)
    elif (outtype == 'd3'):
        if plotname == '':
            htmlfile = __file__[0:__file__.rfind(".")] + '.html'
        else:
            htmlfile = plotname
        print('Plot D3 ', htmlfile)
    elif (outtype == 'html'):
        if plotname == '':
            htmlfile = __file__[0:__file__.rfind(".")] + '.html'
        else:
            htmlfile = plotname
        print('Plot HTML ', htmlfile)
    else:
        print('Plot X')
    # figsize = mpl.rcParams['figure.figsize']

    if (outtype == 'html') or (outtype == 'd3'):
        print('Write html file ', htmlfile)
        fhtml = open(htmlfile,'w')
        print('<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">', file=fhtml)
        print("<html>", file=fhtml)
        print("<head>", file=fhtml)
        print("<title>"+__file__[0:__file__.rfind(".")]+"</title>", file=fhtml)
        print("</head>", file=fhtml)
        print("<body>", file=fhtml)

    ifig = 0

    # -------------------------------------------------------------------------
    # Fig 1
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

    iplot += 1
    xlab   = r'$(0,100)$'
    ylab   = r'$\sin(x)$'
    # if (iplot == 0) | (outtype != 'pdf') | (outtype == 'png') | (outtype == 'html'):
    #     sub  = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))
    #     sub1 = sub
    # else:
    #     # special if windows or d3: zoom one panel zooms all panels
    #     sub = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), sharex=sub1)
    sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

    mark1  = sub.plot(np.sin(np.arange(100)))
    plt.setp(mark1, linestyle='None', marker='o', markeredgecolor=mcol1, markerfacecolor='None',
             markersize=msize, markeredgewidth=mwidth, label=r'UniNorm')

    plt.setp(sub, xlabel=xlab)
    plt.setp(sub, ylabel=ylab)
    sub.grid(False)

    if xlim != None: plt.setp(sub, xlim=xlim)
    if ylim != None: plt.setp(sub, ylim=ylim)

    ll = sub.legend(mark1, [r'$\mathrm{sin\ Nothing}$'], frameon=frameon, ncol=1,
                    labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                    loc='upper left', bbox_to_anchor=(llxbbox,llybbox), scatterpoints=1, numpoints=1)
    plt.setp(ll.get_texts(), fontsize='small')

    ufz.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=usetex, mathrm=True, parenthesis='close')

    # save pages
    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = pngbase+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)
    elif (outtype == 'html'):
        pngfile = htmlfile[0:htmlfile.rfind(".")] + "_" + "{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)
        print('<p><img src='+pngfile+'/></p>', file=fhtml)
    elif (outtype == 'd3'):
        #Does not work: mpld3.plugins.connect(fig, mpld3.plugins.LinkedBrush(line1))
        d3str = mpld3.fig_to_html(fig)
        plt.close(fig)
        print(d3str, file=fhtml)

    # -------------------------------------------------------------------------
    # Fig 2
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

    iplot += 1
    xlab   = r'$(0,1)$'
    ylab   = r'$2\sin(x)$'
    sub    = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))


    line1 = sub.plot(2.*np.sin(np.arange(100)))
    plt.setp(line1, linestyle='-', linewidth=lwidth, marker=None, color=lcol1, label='Gauss')

    plt.setp(sub, xlabel=xlab)
    plt.setp(sub, ylabel=ylab)
    sub.grid(False)

    if xlim != None: plt.setp(sub, xlim=xlim) # set axis limit if wanted
    if ylim != None: plt.setp(sub, ylim=ylim)

    ll = sub.legend(line1, [r'$\mathrm{sin\ Nothing}$'], frameon=frameon, ncol=1,
                    labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                    loc='upper left', bbox_to_anchor=(llxbbox,llybbox), scatterpoints=1, numpoints=1)
    plt.setp(ll.get_texts(), fontsize='small')

    ufz.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=usetex, mathrm=True, parenthesis='close')

    # save pages
    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = pngbase+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)
    elif (outtype == 'html'):
        pngfile = htmlfile[0:htmlfile.rfind(".")] + "_" + "{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)
        print('<p><img src='+pngfile+'/></p>', file=fhtml)
    elif (outtype == 'd3'):
        #Does not work mpld3.plugins.connect(fig, mpld3.plugins.LinkedBrush(line1))
        #mpld3.save_html(fig, d3file)
        d3str = mpld3.fig_to_html(fig)
        plt.close(fig)
        print(d3str, file=fhtml)

    # -------------------------------------------------------------------------
    # Finished
    #

    # close files or show windows
    if (outtype == 'pdf'):
        pdf_pages.close()
    elif (outtype == 'png'):
        pass
    elif (outtype == 'html') or (outtype == 'd3'):
        print("</body>\n</html>", file=fhtml)
        fhtml.close()
    else:
        plt.show()

    t2    = time.time()
    strin = '[m]: '+ufz.astr((t2-t1)/60.,1) if (t2-t1)>60. else '[s]: '+ufz.astr(t2-t1,0)
    print('Time ', strin)
