#!/usr/bin/env python
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


        History
        -------
        Written, MC, May 2012
"""
numeral_map = zip((1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1),
                  ('M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'))

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
        >>> print int2roman(1)
        I

        >>> print int2roman(19)
        XIX

        >>> print int2roman(159)
        CLIX

        >>> print int2roman(159,lower=True)
        clix


        History
        -------
        Written,  MC, May 2012 - Code from Tim Valenta
                                 http://code.activestate.com/recipes/81611-roman-numerals
        Modified, MC, May 2012 - added lower in int2roman
    """
    if i < 1:
        raise ValueError('integer must be > 0.')
    result = []
    for integer, numeral in numeral_map:
        count = int(i / integer)
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
        >>> print roman2int('I')
        1
        
        >>> print roman2int('i')
        1
        
        >>> print roman2int('iv')
        4
        
        >>> print roman2int('MCCCLIV')
        1354


        History
        -------
        Written,  MC, May 2012 - Code from Tim Valenta
                                 http://code.activestate.com/recipes/81611-roman-numerals
    """
    n = unicode(n).upper()
    i = result = 0
    for integer, numeral in numeral_map:
        while n[i:i + len(numeral)] == numeral:
            result += integer
            i += len(numeral)
    return result


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # print int2roman(1)
    # print int2roman(19)
    # print int2roman(159)
    # print int2roman(159,lower=True)
    # print roman2int('I')
    # print roman2int('i')
    # print roman2int('iv')
    # print roman2int('MCCCLIV')
