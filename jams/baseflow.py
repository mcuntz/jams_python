#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

__all__ = ["hollickLyneFilter"]

def hollickLyneFilter(data,beta=0.925,invert=True):
    
    """
    Calculate baseflow from a given discharge timeseries using the 
    digital filter method proposed by Hollick and Lyne (1979).
    

    Definition
    ----------
    def hollickLyneFilter(data,beta=0.925,invert=True)

    Input
    -----
    data   : numpy.ndarray 1D
    
    Optional Input
    --------------
    beta   : float -> dimensionless filter parameter. The default value 
                      was proposed by Nathan and McMahon (1990)
    invert : True  -> add a backward filter pass
           : False -> no backward filterpass, will result in a phase shift

    Output
    ------
    numpy.ndarray 1D

    Literature
    ----------
    Lyne, V. & Hollick, M. 1979, “Stochastic time-variable rainfall-runoff modelling”, 
    Proceedings of the Hydrology and Water Resources Symposium, Perth, 10-12 September, 
    Institution of Engineers National Conference Publication, No. 79/10, pp. 89-92. 

    Nathan, R. J. & McMahon, T. A. 1990a, “Evaluation of automated techniques for base 
    flow and recession analysis”, Water Resources Research , Vol. 26, pp. 1465-1473. 

    Examples
    --------
    >>> import numpy as np
    >>> data = np.array([368.500, 408.100, 493.600, 502.500, 453.100, 439.600, 408.100,
    ...                  430.600, 417.100, 412.600, 390.100, 381.300, 360.300, 340.000,
    ...                  293.100, 276.600, 283.100, 279.800, 270.200, 263.900, 239.800])

    >>> bflow = hollickLyneFilter(data)
    >>> bflow
    array([336.1417221 , 333.15179923, 329.30263055, 324.48886637,
           318.83842976, 312.48099711, 305.46000297, 297.73275208,
           289.26593684, 280.58136362, 272.05823094, 264.05214155,
           257.07123411, 252.24863147, 249.60527727, 247.40663421,
           245.00583555, 242.67971875, 240.70375   ])

    >>> blow = hollickLyneFilter(data,invert=False)
    >>> bflow
    array([336.1417221 , 333.15179923, 329.30263055, 324.48886637,
           318.83842976, 312.48099711, 305.46000297, 297.73275208,
           289.26593684, 280.58136362, 272.05823094, 264.05214155,
           257.07123411, 252.24863147, 249.60527727, 247.40663421,
           245.00583555, 242.67971875, 240.70375   ])

    >>> blow = hollickLyneFilter(data,.8)
    >>> blow
    array([372.46      , 388.138     , 396.45444418, 391.34596522,
           382.70103452, 373.83941116, 364.99301394, 353.31626743,
           338.43908429, 322.71135536, 306.9641942 , 291.00524275,
           276.21905344, 266.1363168 , 261.457896  , 257.59112   ,
           252.5714    , 247.178     , 242.21      ])

    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2015 David Schaefer

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
    Written, David Schaefer, Jun 2015
 
    """
    out = data.copy()
    
    for i in range(1,len(data)):
        bflow = beta * out[i-1] + (1 - beta) * .5 * (data[i]+data[i-1])
        out[i] = min(bflow,data[i])

    if invert:
        out = hollickLyneFilter(out[::-1],beta,invert=False)[::-1]

    return out[1:]


  
if __name__== "__main__":

    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
