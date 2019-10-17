#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

from jams.color import ufzdarkblue, ufzblue, ufzlightblue
from jams.color import ufzred, ufzorange, ufzyellow, ufzdarkgreen, ufzgreen, ufzlightgreen
from jams.color import ufzgray1, ufzgray2, ufzgray3, ufzgrey1, ufzgrey2, ufzgrey3
from jams.color import ufzdarkgray, ufzgray, ufzlightgray, ufzdarkgrey, ufzgrey, ufzlightgrey
from jams.color import ufzblack, ufzwhite
from jams.color import darkblue, blue, lightblue, red, orange, yellow, darkgreen, green, lightgreen
from jams.color import gray1, gray2, gray3, grey1, grey2, grey3, darkgray, gray, lightgray, darkgrey
from jams.color import grey, lightgrey, black, white

__all__ = ['colours', 'colors']

# ---------------------------------------------------------------------

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

# ---------------------------------------------------------------------

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
        (0.0, 62.0, 110.0)

        >>> import numpy as np
        >>> from jams.autostring import astr
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
        return tuple([ i*255. for i in out ])
    else:
        return out

# ---------------------------------------------------------------------

def colours(name=False, rgb=True, rgb256=False, names=False):
    """
        Defines UFZ colours
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
        (0.0, 62.0, 110.0)

        >>> print(colours(names=True)[0:3])
        ['ufzdarkblue', 'ufzblue', 'ufzlightblue']

        >>> import numpy as np
        >>> from jams.autostring import astr
        >>> cc = colours('UFZDARKBLUE')
        >>> print(astr(np.array(cc), 4))
        ['0.0000' '0.2431' '0.4314']

        >>> cc = colours('DarkBlue')
        >>> print(astr(np.array(cc), 4))
        ['0.0000' '0.2431' '0.4314']

        >>> print(colours(['orange','ufzdarkblue'], rgb256=True))
        [(207.0, 104.0, 0.0), (0.0, 62.0, 110.0)]


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2013 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, Jun 2013
        Modified, MC, Jul 2013 - name can be iterable
                               - added black and white
    """
    if names:
        return ufzall

    if name == False:
        raise ValueError('colour must be given.')

    if type(name) == type('string'):
        if name.lower() not in ufzall: return name # raise ValueError('colour not known.')
        return get_colour_tuple(name, rgb256=rgb256)
    else:
        try:
            iterator = iter(name)
        except:
            raise TypeError('colours must be string or iterable.')
        out = list()
        for cname in name:
            if cname.lower() not in ufzall:
                out += [cname] # raise ValueError('colour not known.')
            else:
                out += [get_colour_tuple(cname, rgb256=rgb256)]
        return out

# ---------------------------------------------------------------------

def colors(*args, **kwargs):
    """
        Wrapper function for colours
        def colours(name=False, rgb=True, rgb256=False, names=False):


        Examples
        --------
        >>> print(colors('ufzdarkblue', rgb256=True))
        (0.0, 62.0, 110.0)

        >>> print(colors(names=True)[0:3])
        ['ufzdarkblue', 'ufzblue', 'ufzlightblue']

        >>> import numpy as np
        >>> from jams.autostring import astr
        >>> cc = colors('UFZDARKBLUE')
        >>> print(astr(np.array(cc), 4))
        ['0.0000' '0.2431' '0.4314']

        >>> cc = colors('DarkBlue')
        >>> print(astr(np.array(cc), 4))
        ['0.0000' '0.2431' '0.4314']

        >>> print(colors(['orange','ufzdarkblue'], rgb256=True))
        [(207.0, 104.0, 0.0), (0.0, 62.0, 110.0)]
    """
    return colours(*args, **kwargs)

# ---------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # print(colours('ufzdarkblue', rgb256=True))
    # # (0, 62, 110)
    # print(colours(names=True)[0:3])
    # # ['ufzdarkblue', 'ufzblue', 'ufzlightblue']
    # import numpy as np
    # from jams.autostring import astr
    # cc = colours('UFZDARKBLUE')
    # print(astr(np.array(cc), 4))
    # # ['0.0000' '0.2431' '0.4314']
    # cc = colours('DarkBlue')
    # print(astr(np.array(cc), 4))
    # # ['0.0000' '0.2431' '0.4314']
    # print(colours(['orange','ufzdarkblue'], rgb256=True))
    # #[(207, 104, 0), (0, 62, 110)]
