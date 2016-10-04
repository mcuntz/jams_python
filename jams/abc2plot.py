#!/usr/bin/env python
from __future__ import print_function
from jams.romanliterals import int2roman

def abc2plot(handle, dx, dy, iplot,
             integer=False, roman=False, lower=False,
             parenthesis=None, brackets=None, braces=None,
             small=False, medium=None, large=False,
             xsmall=False, xxsmall=False, xlarge=False, xxlarge=False,
             bold=False, italic=False,
             usetex=False, mathrm=False, string=False, **kwargs):
    """
        Write a, b, c, ...
              A, B, C, ...
              i, ii, iii, ...
              I, II, III, ...
              a), b), c), ... on plots.


        Definition
        ----------
        def abc2plot(handle, dx, dy, iplot,
                     integer=False, roman=False, lower=False,
                     parenthesis=None, brackets=None, braces=None,
                     small=False, medium=None, large=False,
                     xsmall=False, xxsmall=False, xlarge=False, xxlarge=False,
                     bold=False, italic=False,
                     usetex=False, mathrm=False, string=False, **kwargs):


        Input
        -----
        handle   axes handle
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
        small       True:    fontsize='small'
        medium      True:    fontsize='medium' (default)
        large       True:    fontsize='large'
        xsmall      True:    fontsize='x-small'
        xxsmall     True:    fontsize='xx-small'
        xlarge      True:    fontsize='x-large'
        xxlarge     True:    fontsize='xx-large'
        bold        True:    fontweight='bold'
                    False:   fontsize='normal' (default)
        italic      True:    fontstyle='italic'
                    False:   fontstyle='normal' (default)
        usetex      True:    Embed into LaTeX math environment
                    False:   No LaTeX math mode
        mathrm      True:    If usetex=True, surround by \mathrm{}, \mathit{} it italic=true or \mathbf{} if bold=True
                    False:   If usetex=True, standard math font.
        string      True:    Treat iplot as literal string and not as number. integer, roman and lower are disabled.
                    False:   iplot is integer (default)
        **kwargs             All additional parameters passed to axes.text()


        Output
        ------
        Letter/number on plot.


        Restrictions
        ------------
        If output is letter then iplot>26 gives unexpected results.


        Examples
        --------
        abc2plot(sub,0.7,0.6,2,large=True,parenthesis='both',usetex=usetex,mathrm=False)


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

        Copyright 2012-2014 Matthias Cuntz


        History
        -------
        Written,  MC, May 2012
        Modified, AP, Feb 2013 - added parenthesis option
                  MC, Feb 2013 - ported to Python 3
                  MC, Feb 2013 - opening and closing parenthesis, brackets and braces
                  MC, Feb 2013 - usetex, mathrm
                  MC, Oct 2013 - \mathbf, large, default medium
                  MC, Nov 2013 - string, corrected medium default
                               - usetex works with fontsize keyword of axis.text() now (matplotlib v1.1.0)
                  MC, Feb 2014 - if medium is not None then small or large must be set
                  MC, Apr 2014 - assert
                  MC, May 2014 - **kwargs replaces default horizontalalignment='left', verticalalignment='bottom'
                               - italic
                  MC, Oct 2015 - xlarge, xxlarge, xsmall, xxsmall
    """
    # Check input
    assert (roman+integer) < 2, 'either roman literals or integers can be chosen.'
    if medium is None:
        if (small+large+xsmall+xxsmall+xlarge+xxlarge) > 1:
            raise ValueError('only one of xxsmall, xsmall, small, medium, large, xlarge, or xxlarge font size can be chosen (1).')
        elif (small+large+xsmall+xxsmall+xlarge+xxlarge) == 1:
            medium = False
        else:
            medium = True
    else:
        assert (small+medium+large+xsmall+xxsmall+xlarge+xxlarge) <= 1, 'only one of xxsmall, xsmall, small, medium, large, xlarge, or xxlarge font size can be chosen (2).'
        assert (small+medium+large+xsmall+xxsmall+xlarge+xxlarge) > 0, 'If medium=false then another size has to be chosen explicitly: small or large.'
    if usetex & mathrm:
        assert (bold+italic) <= 1, 'if usetex and mathrm then bold and italic are mutually exclusive.'

    # Number or letter
    if string:
        t = str(iplot)
    else:
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
    if small:   fs='small'
    if medium:  fs='medium'
    if large:   fs='large'
    if xsmall:  fs='x-small'
    if xxsmall: fs='xx-small'
    if xlarge:  fs='x-large'
    if xxlarge: fs='xx-large'

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
            elif italic:
                t = '\mathit{' + t + '}'
            else:
                t = '\mathrm{' + t + '}'
        # Works with fontsize keyword now (matplotlib v1.1.0)
        # if small:  t = '\small ' + t
        # if large:  t = '\large ' + t
        t = r'$' + t + r'$'

    xmin, xmax = handle.get_xlim()
    ymin, ymax = handle.get_ylim()
    handle.text(xmin+dx*(xmax-xmin), ymin+dy*(ymax-ymin), t,
                fontsize=fs, fontweight=fw, fontstyle=fst,
                **kwargs) #,
                # horizontalalignment='left', verticalalignment='bottom')


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # # outtype = ''
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
    # abc2plot(sub,0.5,0.6,2,small=True,medium=False,large=False,parenthesis='both',usetex=usetex,mathrm=False)
    # abc2plot(sub,0.6,0.6,2,small=False,medium=True,large=False,parenthesis='both',usetex=usetex,mathrm=False)
    # abc2plot(sub,0.7,0.6,2,small=False,medium=False,large=True,parenthesis='both',usetex=usetex,mathrm=False)
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
    # abc2plot(sub,0.7,0.7,2,medium=True,braces='both', horizontalalignment='left', verticalalignment='bottom')
    # abc2plot(sub,0.8,0.8,2,medium=True,bold=True,braces='none', horizontalalignment='right', verticalalignment='top')

    # if (outtype == 'pdf'):
    #   pdf_pages.savefig(fig)
    #   plt.close()

    # if (outtype == 'pdf'):
    #   pdf_pages.close()
    # else:
    #   plt.show()
