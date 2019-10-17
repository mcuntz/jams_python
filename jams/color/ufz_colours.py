#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Provide UFZ colours.


    Provided colours
    ----------------
    Can be accessed also with colours(list of colours)
    ufzdarkblue, ufzblue, ufzlightblue,
    ufzred, ufzorange, ufzyellow,
    ufzdarkgreen, ufzgreen, ufzlightgreen,
    ufzgray1, ufzgray2, ufzgray3,
    ufzgrey1, ufzgrey2, ufzgrey3,
    ufzdarkgray, ufzgray, ufzlightgray,
    ufzdarkgrey, ufzgrey, ufzlightgrey,
    ufzblack, ufzwhite,
    darkblue=ufzdarkblue, blue=ufzblue, lightblue=ufzlightblue,
    red=ufzred, orange=ufzorange, yellow=ufzyellow,
    darkgreen=ufzdarkgreen, green=ufzgreen, lightgreen=ufzlightgreen,
    gray1=ufzgray1, gray2=ufzgray2, gray3=ufzgray3,
    grey1=ufzgrey1, grey2=ufzgrey2, grey3=ufzgrey3,
    darkgray=ufzdarkgray, gray=ufzgray, lightgray=ufzlightgray,
    darkgrey=ufzdarkgrey, grey=ufzgrey, lightgrey=ufzlightgrey,
    black=ufzblack, white=ufzwhite


    Example
    -------
    >>> import numpy as np
    >>> import jams
    >>> from jams.autostring import astr
    >>> from jams.color.colours import colours
    >>> print(astr(np.array(jams.color.ufzdarkblue), 4))
    ['0.0000' '0.2431' '0.4314']

    >>> print(colours('ufzdarkblue', rgb256=True))
    (0.0, 62.0, 110.0)

    >>> print(colours(names=True)[0:3])
    ['ufzdarkblue', 'ufzblue', 'ufzlightblue']

    >>> print(astr(np.array(jams.color.colours('UFZDARKBLUE')), 4))
    ['0.0000' '0.2431' '0.4314']

    >>> print(astr(np.array(jams.color.colours('DarkBlue')), 4))
    ['0.0000' '0.2431' '0.4314']

    >>> print(jams.color.colours(['orange','ufzdarkblue'], rgb256=True))
    [(207.0, 104.0, 0.0), (0.0, 62.0, 110.0)]


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2015 Matthias Cuntz - mc (at) macu (dot) de

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.


    History
    -------
    Written,  MC, Mar 2015
"""

__all__ = ['ufzdarkblue', 'ufzblue', 'ufzlightblue',
           'ufzred', 'ufzorange', 'ufzyellow',
           'ufzdarkgreen', 'ufzgreen', 'ufzlightgreen',
           'ufzgray1', 'ufzgray2', 'ufzgray3',
           'ufzgrey1', 'ufzgrey2', 'ufzgrey3',
           'ufzdarkgray', 'ufzgray', 'ufzlightgray',
           'ufzdarkgrey', 'ufzgrey', 'ufzlightgrey',
           'ufzblack', 'ufzwhite',
           'darkblue', 'blue', 'lightblue',
           'red', 'orange', 'yellow',
           'darkgreen', 'green', 'lightgreen',
           'gray1', 'gray2', 'gray3',
           'grey1', 'grey2', 'grey3',
           'darkgray', 'gray', 'lightgray',
           'darkgrey', 'grey', 'lightgrey',
           'black', 'white']

# ---------------------------------------------------------------------

ufzdarkblue   = (  0./255.,  62./255., 110./255.)
ufzblue       = (  0./255.,  88./255., 156./255.)
ufzlightblue  = (  0./255., 162./255., 224./255.)
ufzred        = (212./255.,  45./255.,  18./255.)
ufzorange     = (207./255., 104./255.,   0./255.)
ufzyellow     = (230./255., 175./255.,  17./255.)
ufzdarkgreen  = ( 20./255.,  77./255.,  40./255.)
ufzgreen      = (169./255., 181./255.,   9./255.)
ufzlightgreen = ufzgreen
ufzgray1      = ( 81./255.,  81./255.,  81./255.)
ufzgray2      = (156./255., 156./255., 156./255.)
ufzgray3      = (185./255., 185./255., 185./255.)
ufzgrey1      = ufzgray1
ufzgrey2      = ufzgray2
ufzgrey3      = ufzgray3
ufzdarkgray   = ufzgray1
ufzgray       = ufzgray2
ufzlightgray  = ufzgray3
ufzdarkgrey   = ufzgray1
ufzgrey       = ufzgray2
ufzlightgrey  = ufzgray3
ufzblack      = (  0./255.,   0./255.,   0./255.)
ufzwhite      = (255./255., 255./255., 255./255.)

# In Emacs copy ori and query-replace-regexp
# ^\([[:alnum:]]*\).*  ->  \1 = \1
# ^ufz  ->
darkblue   = ufzdarkblue
blue       = ufzblue
lightblue  = ufzlightblue
red        = ufzred
orange     = ufzorange
yellow     = ufzyellow
darkgreen  = ufzdarkgreen
green      = ufzgreen
lightgreen = ufzlightgreen
gray1      = ufzgray1
gray2      = ufzgray2
gray3      = ufzgray3
grey1      = ufzgrey1
grey2      = ufzgrey2
grey3      = ufzgrey3
darkgray   = ufzdarkgray
gray       = ufzgray
lightgray  = ufzlightgray
darkgrey   = ufzdarkgrey
grey       = ufzgrey
lightgrey  = ufzlightgrey
black      = ufzblack
white      = ufzwhite

# ---------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
