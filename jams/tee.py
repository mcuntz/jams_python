#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
tee : Prints arguments on screen and in file, like *nix tee utility.

This module was by Matthias Cuntz while at Department of Computational
Hydrosystems, Helmholtz Centre for Environmental Research - UFZ,
Leipzig, Germany, and continued while at Institut National de Recherche
pour l'Agriculture, l'Alimentation et l'Environnement (INRAE), Nancy, France.

Copyright (c) 2014-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Dec 2014 by Matthias Cuntz (mc (at) macu (dot) de)
* Allow that file keyword not given, Nov 2016, Matthias Cuntz
* Make numpy doctsring format, Dec 2019, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided:

.. autosummary::
   tee
"""


__all__ = ['tee']


def tee(*args, **kwargs):
    """
    Prints arguments on screen and in file, like Unix/Linux tee utility.
        
    Calls print function twice, once with the keyword file and once without, i.e. prints on sys.stdout.

    Parameters
    ----------
    *args : iterable
        All arguments of the print function; will be passed to the `print` function.
    **kwargs : dict
        All keyword arguments of the `print` function; will be passed to the `print` function.
    file : object
        The file argument must be an object with a `write(string)` method.
        If it is not present or `None`, `*args` will be printed on `sys.stdout` only.
        Since printed arguments are converted to text strings, `print()` cannot be used with
        binary mode file objects.arguments.

    Returns
    -------
    None
        If `file` is given and not `None`, then print will be called with `*args` and `**kwargs`
        and a second time with the `file` keyword argument removed, so that `*args` will be
        written to `sys.stdout`.

        If `file` is not given or `None`, `*args` will only be written to `sys.stdout`,
        i.e. ``tee`` is a simple wrapper of the `print` function.

    Examples
    --------
    >>> st = 'Output on screen and in log file'
    >>> ff = 'tee_log.txt'

    >>> f = open(ff, 'w')
    >>> tee(st, file=f)
    Output on screen and in log file
    >>> f.close()

    >>> f = open(ff, 'r')
    >>> test = f.readline()
    >>> f.close()
    >>> test = test[:-1] # rm trailing newline character
    >>> if test == st:
    ...     print('Yes')
    ... else:
    ...     print('No')
    Yes

    >>> import os
    >>> os.remove(ff)

    >>> f=None
    >>> st = 'Output only on screen'
    >>> tee(st, file=f)
    Output only on screen


    History
    -------
    Written,  Matthias Cuntz, Oct 2014
    Modified, Matthias Cuntz, Nov 2016 - file=None
              Matthias Cuntz, Dec 2019 - numpy docstring format
    """
    if 'file' in kwargs:
        if kwargs['file'] is not None:
            print(*args, **kwargs)  # file
        del kwargs['file']
    print(*args, **kwargs)          # screen


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # st = 'Output on screen and in log file'
    # ff = 'tee_log.txt'
    # # write
    # f = open(ff, 'w')
    # tee(st, file=f)
    # f.close()

    # # test
    # f = open(ff, 'r')
    # test = f.readline()
    # f.close()
    # # rm trailing newline character
    # test = test[:-1]
    # if test == st:
    #     print('Yes')
    # else:
    #     print('No')

    # import os
    # os.remove(ff)
