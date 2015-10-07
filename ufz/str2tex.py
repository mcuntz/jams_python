#!/usr/bin/env python
from __future__ import print_function
import numpy as np

__all__ = ['str2tex']

def str2tex(strin, space2linebreak=False, bold=False, italic=False, usetex=True):
    """
        Convert strings to LaTeX strings in math environement used by matplotlib's usetex.

        Strings are embedded into $\mathrm{strin}$ be default but can be embedded into
        \mathbf and \mathit.
        Spaces are escaped but can be replaced by linebreaks.
        
        
        Definition
        ----------
        def str2tex(strin, space2linebreak=False, bold=False, italic=False):


        Input
        -----
        list/ND-array of strings


        Optional Input
        --------------
        space2linebreak   True:  replace space (' ') by linebreak ('\n')
        bold              True:  use \mathbf
        italic            True:  use \mathit
        usetex            False: do only space2linebreak otherwise nothing


        Output
        ------
        list/ND-array of strings that can be used in plots independent of usetex.


        Examples
        --------
        # replace all \ by \\ in docstring in- and outputs
        >>> strin = ['One', 'One-', 'One-Two', 'One Two', 'One\\nTwo', 'A $S_{Ti}$ is great\\nbut use-less']
        >>> print(str2tex(strin))
        ['$\\\\mathrm{One}$', '$\\\\mathrm{One}$$\\\\textrm{-}$', '$\\\\mathrm{One}$$\\\\textrm{-}$$\\\\mathrm{Two}$', '$\\\\mathrm{One\\\\ Two}$', '$\\\\mathrm{One}$ \\n $\\\\mathrm{Two}$', '$\\\\mathrm{A\\\\ }$$S_{Ti}$$\\\\mathrm{\\\\ is\\\\ great}$ \\n $\\\\mathrm{but\\\\ use}$$\\\\textrm{-}$$\\\\mathrm{less}$']
        >>> print(str2tex(strin, bold=True))
        ['$\\\\mathbf{One}$', '$\\\\mathbf{One}$$\\\\textbf{-}$', '$\\\\mathbf{One}$$\\\\textbf{-}$$\\\\mathbf{Two}$', '$\\\\mathbf{One\\\\ Two}$', '$\\\\mathbf{One}$ \\n $\\\\mathbf{Two}$', '$\\\\mathbf{A\\\\ }$$S_{Ti}$$\\\\mathbf{\\\\ is\\\\ great}$ \\n $\\\\mathbf{but\\\\ use}$$\\\\textbf{-}$$\\\\mathbf{less}$']
        >>> print(str2tex(strin, italic=True))
        ['$\\\\mathit{One}$', '$\\\\mathit{One}$$\\\\textit{-}$', '$\\\\mathit{One}$$\\\\textit{-}$$\\\\mathit{Two}$', '$\\\\mathit{One\\\\ Two}$', '$\\\\mathit{One}$ \\n $\\\\mathit{Two}$', '$\\\\mathit{A\\\\ }$$S_{Ti}$$\\\\mathit{\\\\ is\\\\ great}$ \\n $\\\\mathit{but\\\\ use}$$\\\\textit{-}$$\\\\mathit{less}$']
        >>> print(str2tex(strin, space2linebreak=True))
        ['$\\\\mathrm{One}$', '$\\\\mathrm{One}$$\\\\textrm{-}$', '$\\\\mathrm{One}$$\\\\textrm{-}$$\\\\mathrm{Two}$', '$\\\\mathrm{One}$ \\n $\\\\mathrm{Two}$', '$\\\\mathrm{One}$ \\n $\\\\mathrm{Two}$', '$\\\\mathrm{A}$ \\n $\\\\mathrm{}$$S_{Ti}$$\\\\mathrm{ \\n $\\\\mathrm{is \\n $\\\\mathrm{great}$ \\n $\\\\mathrm{but \\n $\\\\mathrm{use}$$\\\\textrm{-}$$\\\\mathrm{less}$']
        >>> print(str2tex(strin, space2linebreak=True, bold=True))
        ['$\\\\mathbf{One}$', '$\\\\mathbf{One}$$\\\\textbf{-}$', '$\\\\mathbf{One}$$\\\\textbf{-}$$\\\\mathbf{Two}$', '$\\\\mathbf{One}$ \\n $\\\\mathbf{Two}$', '$\\\\mathbf{One}$ \\n $\\\\mathbf{Two}$', '$\\\\mathbf{A}$ \\n $\\\\mathbf{}$$S_{Ti}$$\\\\mathbf{ \\n $\\\\mathbf{is \\n $\\\\mathbf{great}$ \\n $\\\\mathbf{but \\n $\\\\mathbf{use}$$\\\\textbf{-}$$\\\\mathbf{less}$']
        >>> print(str2tex(strin, usetex=False))
        ['One', 'One-', 'One-Two', 'One Two', 'One\\nTwo', 'A $S_{Ti}$ is great\\nbut use-less']
        >>> print(str2tex(strin, space2linebreak=True, usetex=False))
        ['One', 'One-', 'One-Two', 'One\\nTwo', 'One\\nTwo', 'A\\n$S_{Ti}$\\nis\\ngreat\\nbut\\nuse-less']


        License
        -------
        This file is part of the UFZ Python package.

        The UFZ Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, Feb 2015
        Modified, MC, Sep 2015 - removed date2dec and dec2date
    """

    # Input type and shape
    if isinstance(strin, list):
        from copy import copy
        istrin  = copy(strin)
    elif isinstance(strin, tuple):
        istrin  = list(strin)
    elif isinstance(strin, np.ndarray):
        istrin   = list(strin.flatten())
    else:
        istrin  = [strin]
    # nstrin = len(istrin)

    # font style
    if (bold+italic) > 1:
        raise ValueError('bold and italic are mutually exclusive.')
    else:
        if bold:
            mtex = r'$\mathbf{'
            ttex = r'$\textbf{'
        elif italic:
            mtex = r'$\mathit{'
            ttex = r'$\textit{'
        else:            
            mtex = r'$\mathrm{'
            ttex = r'$\textrm{'
    
    if usetex:
        a0 = chr(0)
        for j, s in enumerate(istrin):
            if '$' in s:
                ss = s.split('$')
                for ii in range(0,len(ss),2):
                    ss[ii] = mtex+ss[ii]+'}$'
                istrin[j] = '$'.join(ss)
                if s[0] == '$': istrin[j] = istrin[j][11:] # remove leading $\mathrm{}$
            else:
                istrin[j] = mtex+s+'}$'
            # - not minus sign
            if '-' in istrin[j]:
                istrin[j] = istrin[j].replace('-', '}$'+ttex+'-}$'+mtex)
                if istrin[j].endswith('{}$'): istrin[j] = istrin[j][:-11] # remove trailing $\mathrm{}$
            # \n not in tex mode but normal matplotlib
            if '\n' in istrin[j]:
                istrin[j] = istrin[j].replace('\n', '}$'+a0+'\n'+a0+mtex)
            # escape space or linebreak at space
            if ' ' in istrin[j]:
                if space2linebreak:
                    # line break
                    ic = istrin[j].split(' ')
                    for ii, iic in enumerate(ic):
                        if ii==0:
                            istrin[j] = iic + '}$'
                        else:
                            istrin[j] = istrin[j] + a0 + '\n' + a0 + mtex+ iic
                else:                
                    # escaped space 
                    istrin[j] = istrin[j].replace(' ', '\ ')
            # rm ascii character 0 around linebreaks introduced above
            if a0 in istrin[j]:
                istrin[j] = istrin[j].replace(a0, ' ')
    else:
        if space2linebreak:
            istrin = [ i.replace(' ','\n') for i in istrin ]

    # Return right type
    if isinstance(strin, list):
        return istrin
    elif isinstance(strin, tuple):
        return tuple(istrin)
    elif isinstance(strin, np.ndarray):
        return np.array(istrin).reshape(strin.shape)
    else:
        return istrin[0]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # strin = ['One', 'One-', 'One-Two', 'One Two', 'One\nTwo', 'A $S_{Ti}$ is great\nbut use-less']
    # print(str2tex(strin))
    # # ['$\\mathrm{One}$', '$\\mathrm{One}$$\\textrm{-}$', '$\\mathrm{One}$$\\textrm{-}$$\\mathrm{Two}$', '$\\mathrm{One\\ Two}$', '$\\mathrm{One}$ \n $\\mathrm{Two}$', '$\\mathrm{A\\ }$$S_{Ti}$$\\mathrm{\\ is\\ great}$ \n $\\mathrm{but\\ use}$$\\textrm{-}$$\\mathrm{less}$']
    # print(str2tex(strin, bold=True))
    # # ['$\\mathbf{One}$', '$\\mathbf{One}$$\\textbf{-}$', '$\\mathbf{One}$$\\textbf{-}$$\\mathbf{Two}$', '$\\mathbf{One\\ Two}$', '$\\mathbf{One}$ \n $\\mathbf{Two}$', '$\\mathbf{A\\ }$$S_{Ti}$$\\mathbf{\\ is\\ great}$ \n $\\mathbf{but\\ use}$$\\textbf{-}$$\\mathbf{less}$']
    # print(str2tex(strin, italic=True))
    # # ['$\\mathit{One}$', '$\\mathit{One}$$\\textit{-}$', '$\\mathit{One}$$\\textit{-}$$\\mathit{Two}$', '$\\mathit{One\\ Two}$', '$\\mathit{One}$ \n $\\mathit{Two}$', '$\\mathit{A\\ }$$S_{Ti}$$\\mathit{\\ is\\ great}$ \n $\\mathit{but\\ use}$$\\textit{-}$$\\mathit{less}$']
    # print(str2tex(strin, space2linebreak=True))
    # # ['$\\mathrm{One}$', '$\\mathrm{One}$$\\textrm{-}$', '$\\mathrm{One}$$\\textrm{-}$$\\mathrm{Two}$', '$\\mathrm{One}$ \n $\\mathrm{Two}$', '$\\mathrm{One}$ \n $\\mathrm{Two}$', '$\\mathrm{A}$ \n $\\mathrm{}$$S_{Ti}$$\\mathrm{ \n $\\mathrm{is \n $\\mathrm{great}$ \n $\\mathrm{but \n $\\mathrm{use}$$\\textrm{-}$$\\mathrm{less}$']
    # print(str2tex(strin, space2linebreak=True, bold=True))
    # # ['$\\mathbf{One}$', '$\\mathbf{One}$$\\textbf{-}$', '$\\mathbf{One}$$\\textbf{-}$$\\mathbf{Two}$', '$\\mathbf{One}$ \n $\\mathbf{Two}$', '$\\mathbf{One}$ \n $\\mathbf{Two}$', '$\\mathbf{A}$ \n $\\mathbf{}$$S_{Ti}$$\\mathbf{ \n $\\mathbf{is \n $\\mathbf{great}$ \n $\\mathbf{but \n $\\mathbf{use}$$\\textbf{-}$$\\mathbf{less}$']
    # print(str2tex(strin, usetex=False))
    # # ['One', 'One-', 'One-Two', 'One Two', 'One\nTwo', 'A $S_{Ti}$ is great\nbut use-less']
    # print(str2tex(strin, space2linebreak=True, usetex=False))
    # # ['One', 'One-', 'One-Two', 'One\nTwo', 'One\nTwo', 'A\n$S_{Ti}$\nis\ngreat\nbut\nuse-less']
