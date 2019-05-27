#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import jams.const as const

def tcherkez(Rstar, Phi=0.3, T=0.056,
             a2=1.0012, a3=1.0058, a4=1.0161,
             t1=0.9924, t2=1.0008, g=20e-3,
             RG=False, Rchl=False, Rcyt=False, fullmodel=True):
    """
        Calculates the Tcherkez model of 13C-discrimiantion in the Calvin cycle.


        Definition
        ----------
        def tcherkez(Rstar, Phi=0.3, T=0.056, a2=1.0012, a3=1.0058, a4=1.0161,
                     t1=0.9924, t2=1.0008, g=20e-3
                     RG=False, Rchl=False, Rcyt=False, fullmodel=True):


        Input
        -----
        Rstar      Isotope ratio of assimilated carbon, e.g. of Farquhar et al. (1982) model


        Optional Input
        --------------
        Phi       Vo/Vc: ratio of carboxylateion to oxygentation of Rubisco (default: 0.3)
        T         Relative flux of starch synthesis [mol(C6 of starch)/mol(CO2 assimilated)] (default: 0.056)
        a2        Inverse fractionation associated with aldolase
                  for the C-2 position of FBP (Fructose-1,6-bisphosphate) (default: 1.0012)
        a3        Same for C-3 of FBP (default: 1.0058)
        a4        Same for C-4 of FBP (default: 1.0161)
        t1        Inverse fractionation associated with trankelotase
                  for C-1 in E4P (erythrose-4-phosphate) and R5P (ribose-5-phosphate) (default: 0.9924)
        t2        Same for C-2 of X5P (xylulose-5-phosphate) (default: 1.0008)
        g         Isotope discrimination of photorespiratory decarboxylation of Gly (Glycine) (default: 20e-3)
        RG        If True, output isotope ratio of G3P (3-phosphoglyceraldehyde
                  or glyceraldehyde-3-phosphate) (default: False)
        Rchl      If True, output isotope ratio of chloroplastic hexoses and transitory starch (default: False)
        Rcyt      If True, output isotope ratio of cytoplasmic hexoses (default: False)
        fullmodel If True, output RG, Rchl and Rcyt (default: True)


        Output
        ------
        RG, Rchl, Rcyt    if fullmodel=True


        Restrictions
        ------------
        If at least one of RG, Rchl or Rcyt is given then fullmode=False.


        References
        ----------
        Tcherkez G, Farquhar GD, Badeck F & Ghashghaie J, Theoretical considerations about carbon isotope
            distribution in glucose of C3 plants, Functional Plant Biology 31, 857-877, 2004
        Gessler A, Tcherkez G, Peuke AD, Ghashghaie J & Farquhar GD, Experimental evidence for diel variations
            of the carbon isotope composition in leaf, stem and phloem sap organic matter in Ricinus communis,
            Plant, Cell and Environment 31, 941-953, 2004


        Examples
        --------
        >>> a       = -4.4e-3
        >>> b       = -27e-3
        >>> ca    = 353e-6
        >>> ci    = 0.7*ca
        >>> Delta = a+(b-a)*ci/ca
        >>> delta_a1 = -8e-3
        >>> Ra1      = (delta_a1+1.)*const.R13VPDB
        >>> Rstar1   = (1.-Delta)*Ra1
        >>> from autostring import astr
        >>> print(astr((np.array(tcherkez(Rstar1, Phi=0.3, T=0.056))/const.R13VPDB-1.)*1000.,3,pp=True))
        ['12.764' '17.125' '12.978']

        >>> delta_a2 = -7.8e-3
        >>> Ra2      = (delta_a2+1.)*const.R13VPDB
        >>> Rstar2   = (1.-Delta)*Ra2
        >>> R1 = (np.array(tcherkez([Rstar1, Rstar2], Rcyt=True))/const.R13VPDB-1.)*1000.
        >>> print(astr(R1,3,pp=True))
        [['12.978' '13.182']]

        >>> R1, R2 = tcherkez([Rstar1, Rstar2], Rchl=True, Rcyt=True)
        >>> print(astr((R1/const.R13VPDB-1.)*1000.,3,pp=True))
        ['17.125' '17.330']
        >>> print(astr((R2/const.R13VPDB-1.)*1000.,3,pp=True))
        ['12.978' '13.182']


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2012-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jan 2012
        Modified, MC, Feb 2013 - ported to Python 3
    """
    #
    if (RG | Rchl | Rcyt):
       fullmodel = False
    if fullmodel:
       RG   = True
       Rchl = True
       Rcyt = True
    #
    a2tilde = (1.+0.5*Phi-T) / ((2.*a2+1.)/3.+Phi*(2.*a2-0.5)/3.+T*(a2-2.))
    a3tilde = (1.+0.5*Phi-T) / ((2.*a3+1.)/3.+Phi*(2.*a3-0.5)/3.+T*(a3-2.))

    t1tilde = (1.+3.*T)/(t1+3.*T)*t1
    t2tilde = (1.+3.*T)/(t2+3.*T)*t2

    eps     = a3*a3tilde
    epsdash = (t1tilde+1.5*Phi)*a3*a3tilde/(3.*(1.+0.5*Phi-(1.+t2tilde)*a2*a2tilde/3.))

    iRG = np.array(Rstar) / (1.+Phi*(0.5-(1.+g)/(2.+g)*(eps+2.*a2*a2tilde*epsdash)/3.)+T*(a4-1.))

    iRchl = 1./6.*(epsdash*(1.+(a2*a2tilde*t2tilde)/t2)+eps*(2.+t1tilde/t1)+a4) * iRG
    iRcyt = 1./6.*(2.*eps+3.*(a2+1.)/(a2+2.)*epsdash*a2tilde+3.*a3tilde/(2.+a3)*(a3+2.*a4/(1.+a4))) * iRG

    out  = []
    if RG:
       out += [iRG]
    if Rchl:
       out += [iRchl]
    if Rcyt:
       out += [iRcyt]

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

