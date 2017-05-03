#!/usr/bin/env python
from __future__ import print_function
"""
usage: mc_template.py [-h] [-u] [-p plotname] [-s] [-t outtype] [-w] [infile]

This is the Python template for any new program of Matthias Cuntz.

positional arguments:
  infile                Mandatory input file.

optional arguments:
  -h, --help            show this help message and exit
  -p plotname, --plotname plotname
                        Name of plot output file for types pdf, html, d3, bokeh, or plotly,
                        and name basis for type png (default: mc_template).
  -s, --serif           Use serif font; default sans serif.
  -t outtype, --type outtype
                        Output type is pdf, png, html, d3, bokeh, or plotly (default: open
                        screen windows).
  -u, --usetex          Use LaTeX to render text in pdf, png and html.
  -w, --white           White lines on transparent or black background; default: black lines on transparent or white background.


License
-------
This file is part of the JAMS Python library.

The JAMS Python library is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The JAMS Python library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the JAMS Python library.  If not, see <http://www.gnu.org/licenses/>.

Copyright 2012-2016 Matthias Cuntz


History
-------
Written,  MC, Jul 2012
Modified, MC, Jul 2013 - optparse->argparse
          MC, Jul 2013 - extended to be lookup and gallery
          MC, Dec 2013 - add png support
          MC, Mar 2014 - split into individual templates
          MC, Nov 2014 - script -> function
                       - pdf, png, html, or d3
          MC, Sep 2015 - Serif and sans serif fonts
          MC, Dec 2015 - white lines on black or transparent background
          MC, Mar 2017 - bokeh, plotly
"""

# -------------------------------------------------------------------------
# Command line arguments - if script
#

def filebase(f):
    f1 = f
    if f.startswith('..'):
        f1 = f[2:]
    elif f.startswith('.'):
        f1 = f[1:]
    else:
        f1 = f
    if '.' in f1:
        return f[0:f.rfind(".")]
    else:
        return f

# Comment|Uncomment - Begin
if __name__ == '__main__':

    import argparse

    plotname = ''
    outtype  = ''
    usetex   = False
    serif    = False
    dowhite  = False
    parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                       description='''This is the Python template for any new program of Matthias Cuntz.''')
    parser.add_argument('-p', '--plotname', action='store',
                        default=plotname, dest='plotname', metavar='plotname',
                        help='Name of plot output file for types pdf, html, d3, bokeh, or plotly, '
                        'and name basis for type png (default: '+filebase(__file__)+').')
    parser.add_argument('-s', '--serif', action='store_true', default=serif, dest='serif',
                        help='Use serif font; default sans serif.')
    parser.add_argument('-t', '--type', action='store',
                        default=outtype, dest='outtype', metavar='outtype',
                        help='Output type is pdf, png, html, d3, bokeh, or plotly (default: open screen windows).')
    parser.add_argument('-u', '--usetex', action='store_true', default=usetex, dest='usetex',
                        help='Use LaTeX to render text in pdf, png and html.')
    parser.add_argument('-w', '--white', action='store_true', default=dowhite, dest='dowhite',
                        help='White lines on transparent or black background; '
                        'default: black lines on transparent or white background.')
    parser.add_argument('file', nargs='?', default=None, metavar='infile',
                        help='Mandatory input file.')

    args     = parser.parse_args()
    infile   = args.file
    plotname = args.plotname
    outtype  = args.outtype
    serif    = args.serif
    usetex   = args.usetex
    dowhite  = args.dowhite

    del parser, args
# Comment|Uncomment - End

# -------------------------------------------------------------------------
# Function definition - if function
#

# Comment|Uncomment - Begin
# def mc_template(infile=None, plotname='', outtype='', serif=False, usetex=False, dowhite=False):
#     """
#     This is the Python template for any new program of Matthias Cuntz.

#     positional arguments:
#       infile                Mandatory input file.

#     optional arguments:
#       -h, --help            show this help message and exit
#       -p plotname, --plotname plotname
#                             Name of plot output file for types pdf, html, d3, bokeh, or plotly,
#                             and name basis for type png (default: mc_template).
#       -s, --serif           Use serif font; default sans serif.
#       -t outtype, --type outtype
#                             Output type is pdf, png, html, d3, bokeh, or plotly (default: open
#                             screen windows).
#       -u, --usetex          Use LaTeX to render text in pdf, png and html.
#       -w, --white           White lines on transparent or black background; default: black lines on transparent or white background.
#     """
# Comment|Uncomment - End

    # Check input
    if (infile is None):
        raise IOError('\nInput file must be given.\n')

    outtype = outtype.lower()
    outtypes = ['', 'pdf', 'png', 'html', 'd3', 'bokeh', 'plotly']
    if outtype not in outtypes:
        raise IOError('\nOutput '+outype+' type must be in: {:s}'.format(outtypes))

    import numpy as np
    import jams
    import time
    t1 = time.time()

    if dowhite:
        fgcolor = 'white'
        bgcolor = 'black'
    else:
        fgcolor = 'black'
        bgcolor = 'white'

    if (outtype == 'd3'):
        try:
            import mpld3
        except:
            print("No mpld3 found. Use output type html instead.")
            outtype = 'html'

    if (outtype == 'bokeh'):
        try:
            import bokeh.io
            import bokeh.mpl
        except:
            print("No bokeh found. Use output type html instead.")
            outtype = 'html'

    if (outtype == 'plotly'):
        try:
            import plotly.tools
            import plotly.offline
        except:
            print("No plotly found. Use output type html instead.")
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
    mcol1       = jams.color.colours('red')  # primary marker colour
    mcol2       = fgcolor                   # secondary
    mcol3       = jams.color.rgb2rgb01(202,0,32) # third
    if dowhite:
        mcols   = jams.color.colours(['blue','red','lightgray','orange','lightblue',fgcolor])
    else:
        mcols   = jams.color.colours(['blue','red','darkgray','orange','darkblue','black'])
    lcol1       = jams.color.colours('blue')   # primary line colour
    lcol2       = fgcolor
    lcol3       = fgcolor
    lcols       = mcols

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
            if not serif:
                #   r'\usepackage{helvet}',                             # use Helvetica
                mpl.rcParams['text.latex.preamble'] = [
                    r'\usepackage[math,lf,mathtabular,footnotefigures]{MyriadPro}', # use MyriadPro font
                    r'\renewcommand{\familydefault}{\sfdefault}',       # normal text font is sans serif
                    r'\figureversion{lining,tabular}',
                    r'\usepackage{wasysym}',                            # for permil symbol (load after MyriadPro)
                    ]
            else:
                mpl.rcParams['text.latex.preamble'] = [
                    r'\usepackage{wasysym}'                     # for permil symbol
                    ]
        else:
            if serif:
                mpl.rcParams['font.family']     = 'serif'
                mpl.rcParams['font.sans-serif'] = 'Times'
            else:
                mpl.rcParams['font.family']     = 'sans-serif'
                mpl.rcParams['font.sans-serif'] = 'Arial'       # Arial, Verdana
    elif (outtype == 'png') or (outtype == 'html') or (outtype == 'd3') or (outtype == 'bokeh') or (outtype == 'plotly'):
        mpl.use('Agg') # set directly after import matplotlib
        import matplotlib.pyplot as plt
        mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
        if usetex:
            mpl.rc('text', usetex=True)
            if not serif:
                #   r'\usepackage{helvet}',                             # use Helvetica
                mpl.rcParams['text.latex.preamble'] = [
                    r'\usepackage[math,lf,mathtabular,footnotefigures]{MyriadPro}', # use MyriadPro font
                    r'\renewcommand{\familydefault}{\sfdefault}',       # normal text font is sans serif
                    r'\figureversion{lining,tabular}',
                    r'\usepackage{wasysym}',                            # for permil symbol (load after MyriadPro)
                    ]
            else:
                mpl.rcParams['text.latex.preamble'] = [
                    r'\usepackage{wasysym}'                     # for permil symbol
                    ]
        else:
            if serif:
                mpl.rcParams['font.family']     = 'serif'
                mpl.rcParams['font.sans-serif'] = 'Times'
            else:
                mpl.rcParams['font.family']     = 'sans-serif'
                mpl.rcParams['font.sans-serif'] = 'Arial'       # Arial, Verdana
        mpl.rc('savefig', dpi=dpi, format='png')
    else:
        import matplotlib.pyplot as plt
        mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
    mpl.rc('text.latex', unicode=True)
    mpl.rc('font', size=textsize)
    mpl.rc('path', simplify=False) # do not remove
    # print(mpl.rcParams)
    mpl.rc('axes', linewidth=alwidth, edgecolor=fgcolor, facecolor=bgcolor, labelcolor=fgcolor)
    mpl.rc('figure', edgecolor=bgcolor, facecolor='grey')
    mpl.rc('grid', color=fgcolor)
    mpl.rc('lines', linewidth=lwidth, color=fgcolor)
    mpl.rc('patch', edgecolor=fgcolor)
    mpl.rc('savefig', edgecolor=bgcolor, facecolor=bgcolor)
    mpl.rc('patch', edgecolor=fgcolor)
    mpl.rc('text', color=fgcolor)
    mpl.rc('xtick', color=fgcolor)
    mpl.rc('ytick', color=fgcolor)


    # -------------------------------------------------------------------------
    # Data
    #

    pass


    # -------------------------------------------------------------------------
    # Plot
    #

    outtype_ends = ['', '.pdf', '_', '.html', '.html']
    if plotname == '':
        plotfile = filebase(__file__) + outtype_ends[outtypes.index(outtype)]
    else:
        plotfile = plotname
    if outtype == '':
        print('    Plot X')
    else:
        print('    Plot ', plotfile)

    if (outtype == 'pdf'):
        pdf_pages = PdfPages(plotfile)
    # figsize = mpl.rcParams['figure.figsize']

    if outtype in ['html','d3']:
        print('    Write html file ', plotfile)
        fhtml = open(plotfile,'w')
        print('<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">', file=fhtml)
        print("<html>", file=fhtml)
        print("<head>", file=fhtml)
        print("<title>"+plotfile+"</title>", file=fhtml)
        print("</head>", file=fhtml)
        print("<body>", file=fhtml)

    if (outtype == 'bokeh'):
        pass

    if (outtype == 'plotly'):
        htmlfiles = []

    ifig = 0

    # Uncomment for xkcd-style
    # plt.xkcd()

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
    if outtype in ['d3', 'bokeh', 'plotly']:
        xlab   = 'f(x) (0,100)'
        ylab   = 'delta Delta sin(x)'
    else:
        xlab   = jams.str2tex('f(x) (0,100)', usetex=usetex)
        ylab   = jams.str2tex(r'$\delta \Delta \sin(x)$', usetex=usetex)
    # if (iplot == 0) | (outtype == 'pdf') | (outtype == 'png') | (outtype == 'html'):
    #     sub  = fig.add_axes(jams.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))
    #     sub1 = sub
    # else:
    #     # special if windows or d3: zoom one panel zooms all panels
    #     sub = fig.add_axes(jams.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), sharex=sub1)
    sub    = fig.add_axes(jams.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

    mark1  = sub.plot(np.sin(np.arange(100)))
    plt.setp(mark1, linestyle='None', marker='o', markeredgecolor=mcol1, markerfacecolor=None,
             markersize=msize, markeredgewidth=mwidth)

    plt.setp(sub, xlabel=xlab)
    plt.setp(sub, ylabel=ylab)
    sub.grid(False)

    if xlim != None: plt.setp(sub, xlim=xlim)
    if ylim != None: plt.setp(sub, ylim=ylim)

    larr = mark1
    if outtype in ['d3', 'bokeh', 'plotly']:
        tarr = ['sin(x) sin Nothing 100']
    else:
        tarr = [jams.str2tex(r'$\sin(x)$ sin Nothing 100', usetex=usetex)]
    ll = sub.legend(larr, tarr, frameon=frameon, ncol=1,
                    labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                    loc='upper left', bbox_to_anchor=(llxbbox,llybbox), scatterpoints=1, numpoints=1)
    plt.setp(ll.get_texts(), fontsize='small')

    if outtype in ['d3', 'bokeh', 'plotly']:
        jams.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=False, mathrm=True, parenthesis='close')
    else:
        jams.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=usetex, mathrm=True, parenthesis='close')

    # save pages
    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
        plt.close(fig)
    elif (outtype == 'png'):
        pngfile = plotfile+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        plt.close(fig)
    elif (outtype == 'html'):
        pngfile = filebase(plotfile) + "_" + "{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        print('<p><img src='+pngfile+'/></p>', file=fhtml)
        plt.close(fig)
    elif (outtype == 'd3'):
        #Does not work: mpld3.plugins.connect(fig, mpld3.plugins.LinkedBrush(line1))
        d3str = mpld3.fig_to_html(fig)
        print(d3str, file=fhtml)
        plt.close(fig)
    elif (outtype == 'bokeh'):
        htmlfile = plotfile+"{0:04d}".format(ifig)+".html"
        bk = bokeh.mpl.to_bokeh(fig)
        bokeh.io.save(bk, filename=htmlfile, title=htmlfile)
        plt.close(fig)
    elif (outtype == 'plotly'):
        htmlfile = plotfile+"{0:04d}".format(ifig)+".html"
        plotly_fig = plotly.tools.mpl_to_plotly(fig)
        ff = plotly.offline.plot(plotly_fig, filename=htmlfile, auto_open=False)
        htmlfiles.append(htmlfile)
        plt.close(fig)

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
    if outtype in ['d3', 'bokeh', 'plotly']:
        xlab   = '(0,10)'
        ylab   = '2sin(x)'
    else:
        xlab   = jams.str2tex('(0,10)', usetex=usetex)
        ylab   = jams.str2tex(r'2$\sin(x)$', usetex=usetex)
    sub    = fig.add_axes(jams.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

    line1 = sub.plot(2.*np.sin(np.arange(100)))
    plt.setp(line1, linestyle='-', linewidth=lwidth, marker=None, color=lcol1)

    plt.setp(sub, xlabel=xlab)
    plt.setp(sub, ylabel=ylab)
    sub.grid(False)

    if xlim != None: plt.setp(sub, xlim=xlim) # set axis limit if wanted
    if ylim != None: plt.setp(sub, ylim=ylim)

    larr = line1
    if outtype in ['d3', 'bokeh', 'plotly']:
        tarr = ['sin Nothing']
    else:
        tarr = [jams.str2tex(r'$\sin$ Nothing', usetex=usetex)]
    ll = sub.legend(larr, tarr, frameon=frameon, ncol=1,
                    labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                    loc='upper left', bbox_to_anchor=(llxbbox,llybbox), scatterpoints=1, numpoints=1)
    plt.setp(ll.get_texts(), fontsize='small')

    if outtype in ['d3', 'bokeh', 'plotly']:
        jams.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=False, mathrm=True, parenthesis='close')
    else:
        jams.abc2plot(sub, dxabc, dyabc, iplot, lower=True, bold=True, usetex=usetex, mathrm=True, parenthesis='close')

    # save pages
    if (outtype == 'pdf'):
        pdf_pages.savefig(fig)
    elif (outtype == 'png'):
        pngfile = plotfile+"{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    elif (outtype == 'html'):
        pngfile = plotfile[0:plotfile.rfind(".")] + "_" + "{0:04d}".format(ifig)+".png"
        fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
        print('<p><img src='+pngfile+'/></p>', file=fhtml)
    elif (outtype == 'd3'):
        #Does not work mpld3.plugins.connect(fig, mpld3.plugins.LinkedBrush(line1))
        #mpld3.save_html(fig, d3file)
        d3str = mpld3.fig_to_html(fig)
        print(d3str, file=fhtml)
    elif (outtype == 'bokeh'):
        htmlfile = plotfile+"{0:04d}".format(ifig)+".html"
        bk = bokeh.mpl.to_bokeh(fig)
        bokeh.io.save(bk, filename=htmlfile, title=htmlfile)
    elif (outtype == 'plotly'):
        htmlfile = plotfile+"{0:04d}".format(ifig)+".html"
        plotly_fig = plotly.tools.mpl_to_plotly(fig)
        ff = plotly.offline.plot(plotly_fig, filename=htmlfile, auto_open=False)
        htmlfiles.append(htmlfile)
    plt.close(fig)

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
    elif (outtype == 'bokeh'):
        pass
    elif (outtype == 'plotly'):
        if ifig > 1:
            import fileinput
            htmlfile = plotfile+"0001-{0:04d}".format(ifig)+".html"
            with open(htmlfile, 'w') as fout:
                fin = fileinput.input(htmlfiles)
                for line in fin:
                    fout.write(line)
                fin.close()
            import os
            for ff in htmlfiles: os.remove(ff)
    else:
        plt.show()

    t2    = time.time()
    strin = '[m]: '+jams.astr((t2-t1)/60.,1) if (t2-t1)>60. else '[s]: '+jams.astr(t2-t1,0)
    print('Time ', strin)


# -------------------------------------------------------------------------
# Command line usage if function
#

# Comment|Uncomment - Begin
# if __name__ == '__main__':

#     import argparse

#     plotname = ''
#     outtype  = ''
#     usetex   = False
#     serif    = False
#     dowhite  = False
#     parser   = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
#                                        description='''This is the Python template for any new program of Matthias Cuntz.''')
#     parser.add_argument('-p', '--plotname', action='store',
#                         default=plotname, dest='plotname', metavar='plotname',
#                         help='Name of plot output file for types pdf, html, d3, bokeh, or plotly, '
#                         'and name basis for type png (default: '+filebase(__file__)+').')
#     parser.add_argument('-s', '--serif', action='store_true', default=serif, dest="serif",
#                         help="Use serif font; default sans serif.")
#     parser.add_argument('-t', '--type', action='store',
#                         default=outtype, dest='outtype', metavar='outtype',
#                         help='Output type is pdf, png, html, d3, bokeh, or plotly (default: open screen windows).')
#     parser.add_argument('-u', '--usetex', action='store_true', default=usetex, dest="usetex",
#                         help="Use LaTeX to render text in pdf, png and html.")
#     parser.add_argument('-w', '--white', action='store_true', default=dowhite, dest="dowhite",
#                         help="White lines on transparent or black background; default: black lines on transparent or white background.")
#     parser.add_argument('file', nargs='?', default=None, metavar='infile',
#                         help='Mandatory input file.')

#     args     = parser.parse_args()
#     infile   = args.file
#     plotname = args.plotname
#     outtype  = args.outtype
#     serif    = args.serif
#     usetex   = args.usetex
#     dowhite  = args.dowhite

#     del parser, args
    
#     # Call function
#     mc_template(infile, plotname=plotname, outtype=outtype, serif=serif, usetex=usetex, dowhite=dowhite)

# Comment|Uncomment - End
