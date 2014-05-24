#!/usr/bin/env python
from __future__ import print_function

ufzdarkblue   = (0,62,110)
ufzblue       = (0,88,156)
ufzlightblue  = (0,162,224)
ufzred        = (212,45,18)
ufzorange     = (207,104,0)
ufzyellow     = (230,175,17)
ufzdarkgreen  = (20,77,40)
ufzgreen      = (169,181,9)
ufzlightgreen = ufzgreen
ufzgray1      = (81,81,81)
ufzgray2      = (156,156,156)
ufzgray3      = (185,185,185)
ufzgrey1      = ufzgray1
ufzgrey2      = ufzgray2
ufzgrey3      = ufzgray3
ufzdarkgray   = ufzgray1
ufzgray       = ufzgray2
ufzlightgray  = ufzgray3
ufzdarkgrey   = ufzgray1
ufzgrey       = ufzgray2
ufzlightgrey  = ufzgray3
ufzblack      = (0,0,0)
ufzwhite      = (255,255,255)

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

ufzall = ['ufzdarkblue', 'ufzblue', 'ufzlightblue',
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


def get_colour_tuple(name, rgb256=False):
    """
        Helper function for coulours.
        Returns the colour tuple from the global definition.


        Definition
        ----------
        def get_colour_tuple(name, rgb256=False):


        Input
        -----
        name    name of colour


        Optional Input
        --------------
        rgb256     if True: return RGB value tuple between 0 and 255


        Output
        ------
        rgb colour tuple


        Examples
        --------
        >>> print(get_colour_tuple('ufzdarkblue', rgb256=True))
        (0, 62, 110)

        >>> import numpy as np
        >>> from autostring import astr
        >>> cc = get_colour_tuple('UFZDARKBLUE')
        >>> print(astr(np.array(cc), 4))
        ['0.0000' '0.2431' '0.4314']


        History
        -------
        Written,  MC, Jun 2013 - extracted from colours
    """
    d = {}
    exec('out = '+name.lower(), globals(), d)
    out = d['out']
    if rgb256:
        return out
    else:
        return tuple([ i/255. for i in out ])


def colours(name=False, rgb=True, rgb256=False, names=False):
    """
        Defines the UFZ standard colours
            ufzdarkblue, ufzblue, ufzlightblue
            ufzred, ufzorange, ufzyellow
            ufzdarkgreen, ufzgreen, ufzlightgreen
            ufzgray1, ufzgray2, ufzgray3
            ufzgrey1, ufzgrey2, ufzgrey3
            ufzdarkgray, ufzgray, ufzlightgray
            ufzdarkgrey, ufzgrey, ufzlightgrey
        where grey is an alias of gray and 1, 2, 3 aliases dark, none and light.
        All colours exist also without the prefix ufz.


        Definition
        ----------
        def colours(name, rgb=True, rgb256=False, names=False):


        Input
        -----
        name    name(s) of colours


        Optional Input
        --------------
        rgb        if True: return RGB value tuple between 0 and 1 (default)
        rgb256     if True: return RGB value tuple between 0 and 255
        gray       if True: return gray equivalent
        grey       same as gray
        names      print all names


        Output
        ------
        rgb colour tuple(s)


        Examples
        --------
        >>> print(colours('ufzdarkblue', rgb256=True))
        (0, 62, 110)

        >>> print(colours(names=True)[0:3])
        ['ufzdarkblue', 'ufzblue', 'ufzlightblue']

        >>> import numpy as np
        >>> from autostring import astr
        >>> cc = colours('UFZDARKBLUE')
        >>> print(astr(np.array(cc), 4))
        ['0.0000' '0.2431' '0.4314']

        >>> cc = colours('DarkBlue')
        >>> print(astr(np.array(cc), 4))
        ['0.0000' '0.2431' '0.4314']

        >>> print(colours(['orange','ufzdarkblue'], rgb256=True))
        [(207, 104, 0), (0, 62, 110)]


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

        Copyright 2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jun 2013
        Modified, MC, Jul 2013 - name can be iterable
                               - added black and white
    """
    if names:
        return ufzall

    if name == False:
        raise ValueError('colour must be given.')

    if type(name) == type('string'):
        if name.lower() not in ufzall: raise ValueError('colour not known.')
        return get_colour_tuple(name, rgb256=rgb256)
    else:
        try:
            iterator = iter(name)
        except:
            raise TypeError('colours must be string or iterable.')
        out = list()
        for cname in name:
            if cname.lower() not in ufzall: raise ValueError('colour not known.')
            out += [get_colour_tuple(cname, rgb256=rgb256)]
        return out


def colors(*args, **kwargs):
    """
        Wrapper function for colours
        def colours(name=False, rgb=True, rgb256=False, names=False):


        Examples
        --------
        >>> print(colors('ufzdarkblue', rgb256=True))
        (0, 62, 110)

        >>> print(colors(names=True)[0:3])
        ['ufzdarkblue', 'ufzblue', 'ufzlightblue']

        >>> import numpy as np
        >>> from autostring import astr
        >>> cc = colors('UFZDARKBLUE')
        >>> print(astr(np.array(cc), 4))
        ['0.0000' '0.2431' '0.4314']

        >>> cc = colors('DarkBlue')
        >>> print(astr(np.array(cc), 4))
        ['0.0000' '0.2431' '0.4314']

        >>> print(colors(['orange','ufzdarkblue'], rgb256=True))
        [(207, 104, 0), (0, 62, 110)]
    """
    return colours(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # print(colours('ufzdarkblue', rgb256=True))
    # # (0, 62, 110)
    # print(colours(names=True)[0:3])
    # # ['ufzdarkblue', 'ufzblue', 'ufzlightblue']
    # import numpy as np
    # from autostring import astr
    # cc = colours('UFZDARKBLUE')
    # print(astr(np.array(cc), 4))
    # # ['0.0000' '0.2431' '0.4314']
    # cc = colours('DarkBlue')
    # print(astr(np.array(cc), 4))
    # # ['0.0000' '0.2431' '0.4314']
    # print(colours(['orange','ufzdarkblue'], rgb256=True))
    # #[(207, 104, 0), (0, 62, 110)]
