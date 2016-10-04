#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import time

def signature2plot(handle, dx=None, dy=None, name=None,
                   small=False, medium=None, large=False,
                   bold=False, italic=False,
                   usetex=False, mathrm=False,
                   **kwargs):
    """
        Write a copyright notice on a plot:
            (C) name\nCHS/JAMS Leipzig, YYYY
        or if no name given
            (C) CHS/JAMS Leipzig, YYYY
        where YYYY is the 4-digit year.


        Definition
        ----------
        def signature2plot(handle, dx=1, dy=1.05, name=None,
                           small=False, medium=None, large=False,
                           bold=False, italic=False,
                           usetex=False, mathrm=False,
                           **kwargs):

                     
        Input
        -----
        None


        Optional Input
        --------------
        dx          % of xlim from min(xlim) (default: 1)
        dy          % of ylim from min(ylim) (default: 1.05)
        name        Name after (C) (default: None)
        small       True:    fontsize='small'
        medium      True:    fontsize='medium' (default)
        large       True:    fontsize='large'
        bold        True:    fontweight='bold'
                    False:   fontsize='normal' (default)
        italic      True:    fontstyle='italic'
                    False:   fontstyle='normal' (default)
        usetex      True:    Embed into LaTeX math environment
                    False:   No LaTeX math mode
        mathrm      True:    If usetex=True, surround by \mathrm{}, \mathit{} it italic=true or \mathbf{} if bold=True
                    False:   If usetex=True, standard math font.
        **kwargs             All additional parameters passed to axes.text()


        Output
        ------
        (C) name\nCHS/JAMS Leipzig, YYYY
          or
        (C) CHS/JAMS Leipzig, YYYY
          on plot, where YYYY is the 4-digit year.


        Restrictions
        ------------
        None


        Examples
        --------
        >>> None


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, May 2014 - from abc2plot
    """
    
    # Check input
    if medium is None:
        if (small+large) > 1:
            raise ValueError('only one of small, medium or large font size can be chosen (1).')
        elif (small+large) == 1:
            medium = False
        else:
            medium = True
    else:
        assert (small+medium+large) <= 1, 'only one of small, medium or large font size can be chosen (2).'
        assert (small+medium+large) > 0, 'If medium=false then another size has to be chosen explicitly: small or large.'
    if usetex & mathrm:
        assert (bold+italic) <= 1, 'if usetex and mathrm then bold and italic are mutually exclusive.'

    # Construct signature
    # dx/dy given as name
    if type(dx) is str or type(dx) is unicode:
        iname = dx
        dx = None
    elif type(dy) is str or type(dy) is unicode:
        iname = dy
        dy = None
    else:
        iname = name
    
    # default dx/dy
    xmin, xmax = handle.get_xlim()
    ymin, ymax = handle.get_ylim()
    if dx is None:
        if 'transform' in kwargs:
            if kwargs['transform'] is handle.transAxes:
                idx = 0.95
        else:
            idx = xmin + 0.95*(xmax-xmin)
    else:
        idx = dx
    if dy is None:
        if 'transform' in kwargs:
            if kwargs['transform'] is handle.transAxes:
                idy = 0.05
        else:
            idy = ymin + 0.05*(ymax-ymin)
    else:
        idy = dy

    # name
    year = str(time.localtime().tm_year)
    if iname is not None:
        if usetex:
            s1 = r'$\mathrm{\copyright\,'+iname.strip()+'}$'
            if ' ' in s1:
                s1 = s1.replace(' ', '\,')
            if '&' in s1:
                s1 = s1.replace('&', '\&')
            s2 = r'$\mathrm{CHS/JAMS\, Leipzig,\, '+year+'}$'
        else:
            s1 = r'$\copyright$ '+iname.strip()
            s2 = r'CHS/JAMS Leipzig, '+year
    else:
        if usetex:
            s1 = r'$\mathrm{\copyright\, CHS/JAMS\, Leipzig,\, '+year+'}$'
            s2 = None
        else:
            s1 = r'$\copyright$ CHS/JAMS Leipzig, '+year
            s2 = None

    if 'horizontalalignment' not in kwargs:
        kwargs['horizontalalignment'] = 'right'

    if usetex:
        if mathrm:
            if bold:
                s1 = s1.replace('mathrm','mathbf')
                if s2 is not None:
                    s2 = s2.replace('mathrm','mathbf')
            elif italic:
                s1 = s1.replace('mathrm','mathit')
                if s2 is not None:
                    s2 = s2.replace('mathrm','mathit')

    # Size
    if small:  fs='small'
    if medium: fs='medium'
    if large:  fs='large'

    # Weight
    if bold:
        fw='bold'
    else:
        fw='normal'

    # Style
    if italic:
        fst='italic'
    else:
        fst='normal'

    # shift y
    if s2 is None:
        label = handle.text(idx, idy, s1,
                            fontsize=fs, fontweight=fw, fontstyle=fst,
                            **kwargs)
    else:
        # the frame in pixel units
        # transform does not work for polar axis -> take transform_affine
        # x1, y1 = handle.transData.transform(np.array([xmin,ymin]))
        # x2, y2 = handle.transData.transform(np.array([xmax,ymax]))
        x1, y1 = handle.transData.transform_affine(np.array([xmin,ymin]))
        x2, y2 = handle.transData.transform_affine(np.array([xmax,ymax]))
        dpi = handle.figure.dpi         # pixels per inch
        import matplotlib as mpl
        ss  = mpl.rcParams['font.size'] # text size in points: 1 pt = 1/72 inch
        scale_ylim = ymax-ymin
        if 'transform' in kwargs:
            if kwargs['transform'] is handle.transAxes: scale_ylim=1
        if type(handle) is mpl.projections.polar.PolarAxes:
            shifty = 0.25*scale_ylim * ss/72.*float(dpi)/float(y2-y1) # shift by half a textsize
        else:
            shifty = 0.5*scale_ylim * ss/72.*float(dpi)/float(y2-y1)  # shift by half a textsize
        idy1 = idy + shifty
        label = handle.text(idx, idy1, s1,
                            fontsize=fs, fontweight=fw, fontstyle=fst,
                            **kwargs)
        idy2 = idy - shifty
        label = handle.text(idx, idy2, s2,
                            fontsize=fs, fontweight=fw, fontstyle=fst,
                            **kwargs)

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # #outtype = ''
    # outtype = 'pdf'
    # pdffile = 'signature2plot.pdf'
    # usetex  = True
    # textsize = 18

    # import numpy as np
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
    #     mpl.rc('text.latex', unicode=True)
    # else:
    #   import matplotlib.pyplot as plt
    #   mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
    # mpl.rc('path', simplify=False) # do not remove

    # if (outtype == 'pdf'):
    #     print('Plot PDF ', pdffile)
    #     pdf_pages = PdfPages(pdffile)
    # else:
    #     print('Plot X')
    # figsize = mpl.rcParams['figure.figsize']

    # fig = plt.figure()
    # sub = fig.add_axes([0.05,0.05,0.4,0.4])
    # mulx = 10.
    # muly = 20.
    # m = plt.plot(mulx*np.arange(100)/99.,muly*np.arange(100)/99.,'k:')
    # signature2plot(sub,mulx*0.5,muly*0.5, usetex=usetex)
    # signature2plot(sub,mulx*0.6,muly*0.6,'MC',small=True, usetex=usetex) #-
    # signature2plot(sub,mulx*0.7,muly*0.7,large=True,usetex=usetex,italic=True)
    # signature2plot(sub,mulx*0.8,muly*0.8,'MC',bold=True)
    # signature2plot(sub,mulx*0.9,muly*0.9,large=True,usetex=usetex,italic=True, mathrm=True)
    # signature2plot(sub,mulx*1.0,muly*1.0,'MC',usetex=usetex, bold=True, horizontalalignment='right') #-
    # signature2plot(sub,0.9,0.1,'MM',transform=sub.transAxes, usetex=usetex, bold=True, horizontalalignment='right') #-
    # signature2plot(sub,'std',transform=sub.transAxes, usetex=usetex, bold=True, horizontalalignment='right') #-
    # signature2plot(sub,'std', usetex=usetex) #-

    # if (outtype == 'pdf'):
    #   pdf_pages.savefig(fig)
    #   plt.close()

    # if (outtype == 'pdf'):
    #   pdf_pages.close()
    # else:
    #   plt.show()
