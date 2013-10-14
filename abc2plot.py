#!/usr/bin/env python
from __future__ import print_function
from romanliterals import int2roman

def abc2plot(handle, dx, dy, iplot, integer=False, roman=False, lower=False,
             small=False, medium=None, large=False, bold=False,
             parenthesis=None, brackets=None, braces=None, 
             usetex=False, mathrm=False):
    """
        Write a, b, c, ...
              A, B, C, ...
              i, ii, iii, ...
              I, II, III, ...
              a), b), c), ... on plots.

        Definition
        ----------
        def abc2plot(handle, dx, dy, iplot, integer=False, roman=False, lower=False,
                     small=False, medium=True, large=False, bold=False,
                     parenthesis=None, brackets=None, braces=None, 
                     usetex=False, mathrm=False):


        Input
        -----
        dx       % of xlim from min(xlim)
        dy       % of ylim from min(ylim)
        iplot    1=a, 2=b, ...


        Optional Input
        --------------
        integer     True:    use integers
                    False:   use letters (default)
        roman       True:    use roman literals
                    False:   use a, b, c (default)
        lower       True:    use lowercase letters
                    False:   use uppercase letters (default)
        small       True:    fontsize='small'
        medium      True:    fontsize='medium' (default)
        large       True:    fontsize='large'
        bold        True:    fontweight='bold'
                    False:   fontsize='normal' (default)
        parenthesis 'open':  opening parenthesis in front of number
                    'close': closing  parenthesis after number
                    'None':  no parenthesis
                    None:    no parenthesis (default)
        brackets    'open':  opening brackets in front of number
                    'close': closing  brackets after number
                    'None':  no brackets
                    None:    no brackets (default)
        braces      'open':  opening braces in front of number
                    'close': closing  braces after number
                    'None':  no braces
                    None:    no braces (default)
        usetex      True:    Embed into LaTeX math environment
                    False:   No LaTeX math mode
        mathrm      True:    If usetex=True, surround by \mathrm{} or \mathbf{} if bold=True
                    False:   If usetex=True, standard math font.


        Output
        ------
        Letter/number on plot.


        Restrictions
        ------------
        If output is letter then iplot>26 gives unexpected results.

        
        Examples
        --------
        >>> None


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

        Copyright 2012-2013 Matthias Cuntz

        
        History
        -------
        Written,  MC, May 2012
        Modified, AP, Feb 2013 - added parenthesis option
                  MC, Feb 2013 - ported to Python 3
                  MC, Feb 2013 - opening and closing parenthesis, brackets and braces
                  MC, Feb 2013 - usetex, mathrm
                  MC, Oct 2013 - \mathbf, large, default medium
    """
    # Check input
    if (roman & integer):
        raise ValueError('either roman literals or integers can be chosen.')
    if medium is None:
        medium = True
        if (small & large):
            raise ValueError('only one of small, medium or large font size can be chosen.')
        if small: medium = False
        if large: medium = False
    else:
        if (small | large):
            raise ValueError('only one of small, medium or large font size can be chosen.')

    # Number or letter
    if roman:
        t = int2roman(iplot, lower=lower)
    elif integer:
        t = str(iplot)
    else:
        if lower:
            t = chr(96+iplot)
        else:
            t = chr(64+iplot)

    # Size
    if small:  fs='small'
    if medium: fs='medium'
    if large:  fs='large'

    # Weight
    if bold:
        fw='bold'
    else:
        fw='regular'

    # parenthesis, brackets, braces
    if (parenthesis is not None) | (brackets is not None) | (braces is not None):
        if (((parenthesis is not None) & (brackets is not None)) |
            ((parenthesis is not None) & (braces is not None)) |
            ((brackets is not None) & (braces is not None))):
            raise ValueError("parenthesis, brackets and braces mutually exclusive.")

    if parenthesis:
        if parenthesis.lower() == 'open':
            t = '(' + t
        elif parenthesis.lower() == 'close':
            t = t + ')'
        elif parenthesis.lower() == 'both':
            t = '(' + t + ')'
        elif parenthesis.lower() == 'none':
            pass
        else:
            raise ValueError("parenthesis must be either 'open', 'close', 'both', or 'none'.")

    if brackets:
        if brackets.lower() == 'open':
            t = '[' + t
        elif brackets.lower() == 'close':
            t = t + ']'
        elif brackets.lower() == 'both':
            t = '[' + t + ']'
        elif brackets.lower() == 'none':
            pass
        else:
            raise ValueError("brackets must be either 'open', 'close', 'both', or 'none'.")

    if braces:
        if braces.lower() == 'open':
            if usetex:
                t = '\{' + t
            else:
                t = '{' + t
        elif braces.lower() == 'close':
            if usetex:
                t = t + '\}'
            else:
                t = t + '}'
        elif braces.lower() == 'both':
            if usetex:
                t = '\{' + t + '\}'
            else:
                t = '{' + t + '}'
        elif braces.lower() == 'none':
            pass
        else:
            raise ValueError("braces must be either 'open', 'close', 'both', or 'none'.")

    if usetex:
        if mathrm:
            if bold:
                t = '\mathbf{' + t + '}'
            else:
                t = '\mathrm{' + t + '}'
        if small:  t = '\small ' + t
        if large:  t = '\large ' + t
        t = '$' + t + '$'

    xmin, xmax = handle.get_xlim()
    ymin, ymax = handle.get_ylim()
    handle.text(xmin+dx*(xmax-xmin), ymin+dy*(ymax-ymin), t, fontsize=fs, fontweight=fw,
                horizontalalignment='left', verticalalignment='bottom')


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # #outtype = ''
    # outtype = 'pdf'
    # pdffile = 'abc2plot.pdf'
    # usetex  = True
    # textsize = 12
    
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
    # m = plt.plot(range(100),'k:')
    # abc2plot(sub,0,0,2)
    # abc2plot(sub,0.1,0.1,2,parenthesis='close')
    # abc2plot(sub,0.2,0.2,2,lower=True,parenthesis='open')
    # abc2plot(sub,0.3,0.3,2,roman=True,parenthesis='both')
    # abc2plot(sub,0.4,0.4,2,roman=True,lower=True,parenthesis='none')
    # abc2plot(sub,0.5,0.5,2,integer=True,parenthesis='both',usetex=usetex)
    # abc2plot(sub,0.6,0.6,2,small=True,parenthesis='both',usetex=usetex,mathrm=False)
    # abc2plot(sub,0.7,0.7,2,medium=True,brackets='both',usetex=usetex,mathrm=True)
    # abc2plot(sub,0.8,0.8,2,medium=True,bold=True,braces='both',usetex=usetex,mathrm=True)

    # sub = fig.add_axes([0.5,0.5,0.4,0.4])
    # m = plt.plot(range(100),'k:')
    # abc2plot(sub,0.1,0.1,2,brackets='close')
    # abc2plot(sub,0.2,0.2,2,lower=True,brackets='open')
    # abc2plot(sub,0.3,0.3,2,roman=True,brackets='both',usetex=usetex,mathrm=False)
    # abc2plot(sub,0.4,0.4,2,roman=True,lower=True,brackets='none')
    # abc2plot(sub,0.5,0.5,2,integer=True,braces='close')
    # abc2plot(sub,0.6,0.6,2,small=True,braces='open',usetex=usetex,mathrm=True)
    # abc2plot(sub,0.7,0.7,2,medium=True,braces='both')
    # abc2plot(sub,0.8,0.8,2,medium=True,bold=True,braces='none')

    # if (outtype == 'pdf'):
    #   pdf_pages.savefig(fig)
    #   plt.close()
      
    # if (outtype == 'pdf'):
    #   pdf_pages.close()
    # else:
    #   plt.show()
