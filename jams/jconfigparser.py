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

    Copyright 2018 Matthias Cuntz


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
