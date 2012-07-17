#!/usr/bin/env python
from romanliterals import int2roman

def abc2plot(handle, dx, dy, iplot, roman=False, lower=False, integer=False,
             medium=False, small=False, bold=False):
    """
        Write a, b, c, ...
              A, B, C, ...
              i, ii, iii, ...
              I, II, III, ... on plots.

        Definition
        ----------
        def abc2plot(handle, dx, dy, iplot, roman=False, lower=False, integer=False,
                     medium=False, small=False):


        Input
        -----
        dx       % of xlim from min(xlim)
        dy       % of ylim from min(ylim)
        iplot    1=a, 2=b, ...


        Optional Input
        --------------
        roman    True:  use roman literals
                 False: use a, b, c
        lower    True:  use lowercase letters
                 False: use uppercase letters
        integer  True:  use integers
                 False: use letters
        medium   True:  fontsize='medium'
                 False: fontsize='large'
        small    True:  fontsize='small'
                 False: fontsize='large'
        bold     True:  fontweight='bold'
                 False: fontsize='normal'


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

        Copyright 2012 Matthias Cuntz

        
        History
        -------
        Written, MC, May 2012
    """
    # Check input
    if (roman & integer):
        raise ValueError('either roman literals or integers can be chosen.')
    if (small & medium):
        raise ValueError('either small or medium font size can be chosen.')
    if roman:
        t = int2roman(iplot, lower=lower)
    elif integer:
        t = str(iplot)
    else:
        if lower:
            t = chr(96+iplot)
        else:
            t = chr(64+iplot)
    if small:
        fs='small'
    elif medium:
        fs='medium'
    else:
        fs='large'
    if bold:
        fw='bold'
    else:
        fw='regular'
    xmin, xmax = handle.get_xlim()
    ymin, ymax = handle.get_ylim()
    handle.text(xmin+dx*(xmax-xmin), ymin+dy*(ymax-ymin), t, fontsize=fs, fontweight=fw,
                horizontalalignment='left', verticalalignment='bottom')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # import matplotlib.pyplot as plt
    # fig = plt.figure()
    # sub = fig.add_axes([0.1,0.1,0.9,0.9])
    # m = plt.plot(range(100),'k:')
    # abc2plot(sub,0,0,2)
    # abc2plot(sub,0.1,0.1,2)
    # abc2plot(sub,0.2,0.2,2,lower=True)
    # abc2plot(sub,0.3,0.3,2,roman=True)
    # abc2plot(sub,0.4,0.4,2,roman=True,lower=True)
    # abc2plot(sub,0.5,0.5,2,integer=True)
    # abc2plot(sub,0.6,0.6,2,small=True)
    # abc2plot(sub,0.7,0.7,2,medium=True)
    # abc2plot(sub,0.7,0.7,2,medium=True,bold=True)
    # plt.show()
