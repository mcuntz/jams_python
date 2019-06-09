#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
'''
    Extension of Python's ConfigParser: adding getlist method.

    Using getlist from NoahTheDuke:
        https://gist.github.com/NoahTheDuke/e6d282b421f6a126062e81696e4cfc2a
    but extending Python's ConfigParser standard module.


    Example
    -------
    From
        https://www.reddit.com/r/Python/comments/89c9b5/what_config_file_format_do_you_prefer/
    File myconfig.ini:
    [params]
    years: 08,09,10,11,12,13,14,15,16
    categories:
        5_000
        12_500
        17_500
        22_500
    names: [Steve John Brian]
    years and names: ${params:years} ${params:names}

    >>> from jconfigparser import jConfigParser
    >>> jconfig = jConfigParser()
    >>> jconfig.read('myconfig.ini')
    >>> jconfig.getiter('params', 'years')
    <filter at 0x19ef8641128>

    >>> jconfig.getlist('params', 'years')
    ['08', '09', '10', '11', '12', '13', '14', '15', '16']

    >>> list(map(int, jconfig.getiter('params', 'categories')))
    [5000, 12500, 17500, 22500]

    >>> jconfig.getlist('params', 'names')
    ['Steve', 'John', 'Brian']

    >>> jconfig.getlist('params', 'years and names')
    ['08', '09', '10', '11', '12', '13', '14', '15', '16', 'Steve', 'John', 'Brian']


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT License.

    Copyright (c) 2018 Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  Matthias Cuntz, Nov 2018
    Modified, Matthias Cuntz, Nov 2018 - getlist, getiter
                                       - sep
'''
try:    # Python 2
    from ConfigParser import ConfigParser
except: # Python 3
    from configparser import ConfigParser
import itertools as it
import re

__all__ = ['jConfigParser']

# https://stackoverflow.com/questions/1713038/super-fails-with-error-typeerror-argument-1-must-be-type-not-classobj-when
class jConfigParser(ConfigParser, object):
    def __init__(self, *args, **kwargs):
        super(jConfigParser, self).__init__(*args, **kwargs)

    # https://gist.github.com/NoahTheDuke/e6d282b421f6a126062e81696e4cfc2a
    def getiter(self, section, option, sep='[ ,\[\]]'):
        ''' Return an iterator over the splitted option.
            Separator sep can be given, otherwise pattern '[ ,\[\]]' will be used. '''
        lines = self.get(section, option).splitlines()
        if len(lines) == 1:
            return filter(None, re.split(sep, lines[0]))
        else:
            return filter(None, it.chain.from_iterable(map(lambda x: re.split(sep, x), lines)))

    def getlist(self, section, option, sep='[ ,\[\]]'):
        ''' Return a list of the splitted option.
            Separator sep can be given, otherwise pattern '[ ,\[\]]' will be used. '''
        lines = self.get(section, option).splitlines()
        if len(lines) == 1:
            return re.split(sep, lines[0])
        else:
            return list(map(lambda x: re.split(sep, x), lines))
