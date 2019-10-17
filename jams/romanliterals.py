#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
"""
        Convert integer to and from roman numerals.


        Routines
        ----------
        int2roman    integer to roman numeral
        roman2int    roman numral to integer


        References
        ------------
        http://code.activestate.com/recipes/81611-roman-numerals


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2014 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, May 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
"""

__all__ = ['int2roman', 'roman2int']

numeral_map = list(zip((1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1),
                  ('M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I')))

def int2roman(i, lower=False):
    """
        Convert an integer to a roman numeral.

        Definition
        ----------
        def int2roman(i, lower=False):


        Input
        -----
        i         integer


        Optional Input
        --------------
        lower    True:  output lowercase numerals
                 False: output uppercase numerals


        Output
        ------
        Roman numeral


        Restrictions
        ------------
        Input can be only single integer.


        References
        ------------
        http://code.activestate.com/recipes/81611-roman-numerals


        Examples
        --------
        >>> print(int2roman(1))
        I

        >>> print(int2roman(19))
        XIX

        >>> print(int2roman(159))
        CLIX

        >>> print(int2roman(159,lower=True))
        clix


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2014 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, May 2012 - Code from Tim Valenta
                                 http://code.activestate.com/recipes/81611-roman-numerals
        Modified, MC, May 2012 - added lower in int2roman
                  MC, Feb 2013 - ported to Python 3
                  MC, Apr 2014 - assert
    """
    assert i >= 1, 'integer must be > 0.'
    result = []
    for integer, numeral in numeral_map:
        count = int(i // integer)
        result.append(numeral * count)
        i -= integer * count
    if lower: result = [ i.lower() for i in result ]
    return ''.join(result)


def roman2int(n):
    """
        Convert a roman numeral to an integer.

        Definition
        ----------
        def roman2int(i):


        Input
        -----
        i         string with roman numeral


        Output
        ------
        Integer


        Restrictions
        ------------
        Input can be only single numeral.


        References
        ------------
        http://code.activestate.com/recipes/81611-roman-numerals


        Examples
        --------
        >>> print(roman2int('I'))
        1

        >>> print(roman2int('i'))
        1

        >>> print(roman2int('iv'))
        4

        >>> print(roman2int('MCCCLIV'))
        1354


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2012-2013 Matthias Cuntz - mc (at) macu (dot) de

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
        Written,  MC, May 2012 - Code from Tim Valenta
                                 http://code.activestate.com/recipes/81611-roman-numerals
        Modified, MC, Feb 2013 - ported to Python 3
    """
    n = str(n).upper()
    i = result = 0
    for integer, numeral in numeral_map:
        while n[i:i + len(numeral)] == numeral:
            result += integer
            i += len(numeral)
    return result


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    # print int2roman(1)
    # print int2roman(19)
    # print int2roman(159)
    # print int2roman(159,lower=True)
    # print roman2int('I')
    # print roman2int('i')
    # print roman2int('iv')
    # print roman2int('MCCCLIV')

