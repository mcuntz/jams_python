#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

__all__ = ['mc_set_outtype', 'mc_set_matplotlib',
           'mc_plot_begin', 'mc_plot_start', 'mc_plot_save', 'mc_plot_end', 'mc_plot_stop']

mc_set_outtype = '''
outtype = outtype.lower()
outtypes = ['', 'pdf', 'png', 'html', 'd3', 'bokeh', 'plotly']
if outtype not in outtypes:
    raise IOError('\\nOutput '+outype+' type must be in: {:s}'.format(outtypes))

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
    if (outtype == 'plotly') and (plotname != ''):
        assert plotname.endswith('html'), 'Plotly plotnames must end with .html'
'''

mc_set_matplotlib = '''
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
            #   r'\\usepackage{helvet}',                             # use Helvetica
            mpl.rcParams['text.latex.preamble'] = [
                r'\\usepackage[math,lf,mathtabular,footnotefigures]{MyriadPro}', # use MyriadPro font
                r'\\renewcommand{\\familydefault}{\\sfdefault}',       # normal text font is sans serif
                r'\\figureversion{lining,tabular}',
                r'\\usepackage{wasysym}',                            # for permil symbol (load after MyriadPro)
                ]
        else:
            mpl.rcParams['text.latex.preamble'] = [
                r'\\usepackage{wasysym}'                     # for permil symbol
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
            #   r'\\usepackage{helvet}',                             # use Helvetica
            mpl.rcParams['text.latex.preamble'] = [
                r'\\usepackage[math,lf,mathtabular,footnotefigures]{MyriadPro}', # use MyriadPro font
                r'\\renewcommand{\\familydefault}{\\sfdefault}',       # normal text font is sans serif
                r'\\figureversion{lining,tabular}',
                r'\\usepackage{wasysym}',                            # for permil symbol (load after MyriadPro)
                ]
        else:
            mpl.rcParams['text.latex.preamble'] = [
                r'\\usepackage{wasysym}'                     # for permil symbol
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
'''

mc_plot_begin = '''
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
'''

mc_plot_start = mc_plot_begin

mc_plot_save = '''
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
'''
        
mc_plot_end = '''
# close files or show windows
if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
elif (outtype == 'html') or (outtype == 'd3'):
    print("</body>\\n</html>", file=fhtml)
    fhtml.close()
elif (outtype == 'bokeh'):
    pass
elif (outtype == 'plotly'):
    import os
    if ifig > 1:
        import fileinput
        htmlfile = plotfile
        with open(htmlfile, 'w') as fout:
            fin = fileinput.input(htmlfiles)
            for line in fin:
                fout.write(line)
            fin.close()
        for ff in htmlfiles: os.remove(ff)
    else:
        os.rename(htmlfiles[0], plotfile)
else:
    plt.show()
'''

mc_plot_stop = mc_plot_end
