#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import jams.const as const

def cuntz_gleixner(idecdate, iGPP, iRd, iCa, iRa, igtot, sunrise, Vcyt=None,
                   date0=False,
                   V0starch=const.eps,
                   R0starch=const.R13VPDB,
                   R0cyt=const.R13VPDB,
                   daynight=None, daylength=57600,
                   Phi=0.3, s_resid=const.eps,
                   betas=None, betap=0.75,
                   epsa=4.4e-3, epsb=29.5e-3,
                   epsg=20.0e-3, epst=-4.4e-3,
                   epss=10.0e-3, epsp=1.0e-3,
                   steady=False,
                   Rass=False, Rm=False, Rchl=False,
                   Rcyt=False, Rstarch=False,
                   Rpyr=False, Rbio=False,
                   Rphloem=False,
                   Vstarch=False, ass13=False, disc=False,
                   Rnew_starch=False, Rnew_cyt=False,
                   fullmodel=True, julian=True, nocheck=False,
                   starch_mol2g=None, V0starchg=const.eps):
    """
        Calculates the Cuntz-Gleixner steady state and non-steady state models
        of 13C discrimiantion in the Calvin cycle.


        Definition
        ----------
        def cuntz_gleixner(idecdate, iGPP, iRd, iCa, iRa, igtot, sunrise, Vcyt=None,
                           date0=False,
                           V0starch=const.eps, R0starch=const.R13VPDB,
                           R0cyt=const.R13VPDB,
                           daynight=None, daylength=57600,
                           Phi=0.3, s_resid=const.eps,
                           betas=None, betap=0.75,
                           epsa=4.4e-3, epsb=29.5e-3,
                           epsg=20.0e-3, epst=-4.4e-3,
                           epss=10.0e-3, epsp=1.0e-3,
                           steady=False,
                           Rass=False, Rm=False, Rchl=False,
                           Rcyt=False, Rstarch=False,
                           Rpyr=False, Rbio=False,
                           Rphloem=False,
                           Vstarch=False, ass13=False, disc=False,
                           Rnew_starch=False, Rnew_cyt=False,
                           fullmodel=True, julian=True, nocheck=False,
                           starch_mol2g=None, V0starchg=const.eps):


        Input
        -----
        idecdate         decimal date
        iGPP             Gross Photosynthesis: GPP = A + Rd [umol(CO2)/m2s]
        iRd              Leaf respiration: GPP - A [umol(CO2)/m2s]
        iCa              Outside CO2 concentration [ppm=umol(CO2)/umol(air)]
        iRa              13C/12C ratio of outside CO2 concentration
        igtot            Total conductance for CO2 from outside air to chloroplast [mol(CO2)/m2s]
        sunrise         decial date of first sunrise in data set


        Input (only nss model)
        -----
        Vcyt            C-concentration of sucrose pool in cytoplasm [umol(C)/m2(leaf)]


        Optional Input
        --------------
        date0           Start date of 1st time step (default: False)
                        If False, take same time step as first time step in idecdate
        V0starch        Initial C-concentration in Starch [umol(C)/m2(leaf)] (default: 1e-6)
        V0starchg       Initial C-concentration in Starch [g(C)/gDW] (default: 1e-6)
        starch_mol2g    Conversion factor from [umol(C)/m2(leaf)] to [g(C)/gDW] used for starch (default: None)
                        If given and Vstarch==True then Vstarchg will be returned as well
        R0starch        Initial 13C/12C ratio of starch pool (default: PDB)
        R0cyt           Initial 13C/12C ratio of C in cytoplasm (default: PDB)
        daynight        1/0 array of day or night (default: False)
                        If False, day is when gpp>0
        daylength       length of daylight [s] (default: 57600 = 16h)
        Phi             Vc/Vo, ratio of carboxylation to oxygenation of Rudisco (default: 0.3)
        s_resid         Residual starch concentration at end of night [umol(C)/m2(leaf)] (default: 1e-6)
        betas           factor of leaf respiration transferred to biosynthesis (default: False)
                        If False, betas is 3*gpp/max(gpp)
                        Note: betas*(1-betap) <= 1: if betap=2/3 -> betas<3: if betap=5/6 -> betas < 6
        betap           fraction of respiration occuring during biosynthesis: min=2/3; max=4/5 (default: 3/4)
        epsa            effective fractionation along gtot (default: 4.4e-3)
        epsb            fractionation of Rubisco (default: 29.5e-3)
        epsg            fractionation of photorespiration (default: 20e-3)
        epst            equilibrium fractionation value for starch synthesis (default: -4.4e-3)
        epss            fractionation of biosynthesis production (default: 10e-3)
        epsp            fractionation of biosynthesis bifurcation (default: 1e-3)


        Parameter
        ---------
        steady        If True, steady-state instead of non-steady-state model (default: False)
        Rass          If True, output 13C/12C ratio of assimilated carbon (default: False)
        Rm            If True, output 13C/12C ratio of chloroplast CO2 (default: False)
        Rchl          If True, output 13C/12C ratio of sugars in chloroplast (default: False)
        Rcyt          If True, output 13C/12C ratio of sugars in cytoplasm (default: False)
        Rstarch       If True, output 13C/12C ratio of starch (default: False)
        Rpyr          If True, output 13C/12C ratio of sugars at pyruvate pathway (default: False)
        Rbio          If True, output 13C/12C ratio of biosynthesis products before bifurcation (default: False)
        Rphloem       If True, output 13C/12C ratio of new phloem products (sugars & biosynthesis products) (default: False)
        Vstarch       If True, output C-concentration in Starch [umol(C)/m2(leaf)] (default: False)
        ass13         If True, output 13C assimilation rate  [umol(13C)/m2s] (default: False)
        disc          If True, output Discrimination 1-Rass/Ra (default: False)
        Rnew_starch   If True, output 13C/12C ratio of newly produced starch (default: False)
        Rnew_cyt      If True, output 13C/12C ratio of newly produced sugars in cytoplasm (default: False)
        fullmodel     If True, output all in the above order (default: True)
        julian        If True, dates are given as Julian days, otherwise as decimal year (default: True)
        nocheck       If True, do not check betap and betas ranges (default: False)


        Output
        ------
        if fullmodel=True
          Rass, Rm, Rchl, Rcyt, Rstarch, Rpyr, Rbio, Rphloem, Vstarch, ass13, disc, Rnew_starch, Rnew_cyt
        and if Vstarch=True and starch_mol2g!=None
          Vstarchg


        Restrictions
        ------------
        If at least one individual output parameter is True then fullmode=False.


        References
        ----------
        Tcherkez G, Farquhar GD, Badeck F & Ghashghaie J, Theoretical considerations about carbon isotope
           distribution in glucose of C3 plants, Functional Plant Biology 31, 857-877, 2004
        Gessler A, Tcherkez G, Peuke AD, Ghashghaie J & Farquhar GD, Experimental evidence for diel variations
           of the carbon isotope composition in leaf, stem and phloem sap organic matter in Ricinus communis,
           Plant, Cell and Environment 31, 941-953, 2004


        Examples
        --------
        # steady state
        >>> adecdate = np.array([2008.918658925319050,2008.918772768671033,2008.918886612022106,
        ...                      2008.919000455374089,2008.919114298724935,2008.919228142076918,
        ...                      2008.919341985427991,2008.919455828779974,2008.919569672131956,
        ...                      2008.919683515483030,2008.919797358835012,2008.919911202186086,
        ...                      2008.920025045538068,2008.920138888888914,2008.920252732240897,
        ...                      2008.920366575591970,2008.920480418943953,2008.920556314511941,
        ...                      2008.920594262295026,2008.920708105647009,2008.920821948998992,
        ...                      2008.920935792350065,2008.921049635702047,2008.921163479052893,
        ...                      2008.921277322405103,2008.921391165755949,2008.921505009107932])
        >>> gpp = np.array([0.000000000000,23.700827991217,22.449718259243,21.253578109071,20.222525197027,
        ...                   19.503625355216,18.797132965271,18.102416224453,17.780887860470,17.491607940331,
        ...                   17.207072197663,17.089915139494,17.995854647885,18.901914959729,19.681631460738,
        ...                   19.681631460738,19.681631460738,0.000000000000,0.000000000000,0.000000000000,
        ...                   0.000000000000,0.000000000000,0.000000000000,0.000000000000,0.000000000000,
        ...                   0.000000000000] )
        >>> Rd = np.array([0.511900000000,2.361686144687,2.743373026721,3.180029474251,3.476842651940,
        ...                  3.259038512076,3.053641828083,2.860020793216,2.958750580931,3.083603451827,
        ...                  3.213200496886,3.331826587704,3.352936200975,3.374166608865,3.392531460738,
        ...                  3.392531460738,3.392531460738,1.025929405070,0.829676977143,0.633424549217,
        ...                  0.437172122858,0.303515021488,0.547877741613,0.792240464668,1.036603184794,
        ...                  1.280965907360] )
        >>> CO2air = np.array([620.902600000000,537.510500000000,608.806500000000,671.251000000000,
        ...                      652.204000000000,560.157800000000,427.130100000000,395.276000000000,
        ...                      427.000400000000,410.953300000000,386.943500000000,500.417500000000,
        ...                      552.776800000000,515.865800000000,542.450400000000,692.503500000000,
        ...                      656.423500000000,588.844100000000,675.156500000000,725.101900000000,
        ...                      664.837000000000,598.080600000000,610.713600000000,487.087000000000,
        ...                      531.921300000000,675.177700000000] )
        >>> Ra = np.array([0.011067265443,0.011083081802,0.011071245659,0.011060401761,0.011063313320,
        ...                  0.011080216316,0.011111970396,0.011122420992,0.011111174802,0.011116914764,
        ...                  0.011125605614,0.011097923896,0.011079382516,0.011087211473,0.011083896499,
        ...                  0.011057329511,0.011062335683,0.011072518834,0.011061590657,0.011053508863,
        ...                  0.011061281634,0.011071628848,0.011069690431,0.011093962783,0.011086022577,
        ...                  0.011059558971] )
        >>> gtot = np.array([0.064395001124,0.074054920058,0.078085762302,0.078484156864,0.078127160737,
        ...                    0.085209848990,0.088685679784,0.089611189047,0.088528095110,0.086087621579,
        ...                    0.081901616151,0.076984314568,0.080693530135,0.084173028182,0.087005780100,
        ...                    0.087005780100,0.087005780100,0.046798889383,0.042324852911,0.037583815518,
        ...                    0.032460459750,0.028193059760,0.031985237181,0.035564641600,0.038983725824,
        ...                    0.042282334176] )
        >>> ndecdate = 2008.918772768670806
        >>> V0starch = 60498.901260546168
        >>> R0starch = 0.010949362493
        >>> daynight = np.array([0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0])
        >>> daylength = 56400.
        >>> Phi = np.array([0.081107124487,0.228537911761,0.178664767959,0.161131279216,0.173373369802,
        ...                   0.195820137312,0.276226485398,0.291675965181,0.259480203313,0.285111049806,
        ...                   0.331368464761,0.228636707479,0.199019858196,0.224041328002,0.209848472684,
        ...                   0.147530335706,0.158874834596,0.090595209673,0.077545142771,0.070693599170,
        ...                   0.075277039328,0.082213420279,0.081562178107,0.102842605886,0.095940191411,
        ...                   0.077688517811] )
        >>> s_resid = 40532.561476983901
        >>> betas = np.array([0.000000000000,4.999999000000,4.736061020188,4.483719593421,4.266205481166,
        ...                     4.114544146364,3.965500325995,3.818940732998,3.751110114569,3.690082727626,
        ...                     3.630056187040,3.605340371581,3.796460173119,3.987605459727,4.152096865111,
        ...                     4.152096865111,4.152096865111,0.000000000000,0.000000000000,0.000000000000,
        ...                     0.000000000000,0.000000000000,0.000000000000,0.000000000000,0.000000000000,
        ...                     0.000000000000] )
        >>> betap = 0.8
        >>> epsa = np.array([0.002995512907,0.003039740417,0.003192366495,0.003375544906,0.003479516543,
        ...                    0.003318658489,0.003187204366,0.003078822327,0.003144673535,0.003235857882,
        ...                    0.003343598868,0.003447509132,0.003408830551,0.003373557382,0.003345605523,
        ...                    0.003345605523,0.003345605523,0.003500242606,0.003551923813,0.003615033761,
        ...                    0.003693231875,0.003766948014,0.003706681514,0.003655727827,0.003612282052,
        ...                    0.003574980671] )
        >>> epsb = 0.029
        >>> epsg = 0.0185
        >>> epst = -0.004
        >>> epss = 0.01
        >>> epsp = 0.003
        >>> [Vstarch] = cuntz_gleixner(adecdate[1:], gpp, Rd, CO2air, Ra, gtot, ndecdate,
        ...                            date0=adecdate[0], V0starch=V0starch, R0starch=R0starch, daynight=daynight,
        ...                            daylength=daylength, Phi=Phi, s_resid=s_resid,
        ...                            betas=betas, betap=betap, epsa=epsa, epsb=epsb,
        ...                            epsg=epsg, epst=epst, epss=epss, epsp=epsp,
        ...                            steady=True, Vstarch=True, julian=False)
        >>> from autostring import astr
        >>> print(astr(Vstarch[0:6],5,pp=True))
        [' 40532.56148' ' 70158.59647' ' 98220.74429' '124787.71693' '150065.87342' '174445.40512']
        >>> [Rass,Rm,Rchl,Rcyt,Rstarch,Rpyr,Rbio,Rphloem,Vstarch,ass13,disc,Rnew_starch,Rnew_cyt,Vstarchg] = cuntz_gleixner(
        ...                adecdate[1:], gpp, Rd, CO2air, Ra, gtot, ndecdate,
        ...                date0=adecdate[0], V0starch=V0starch, R0starch=R0starch, daynight=daynight,
        ...                daylength=daylength, Phi=Phi, s_resid=s_resid,
        ...                betas=betas, betap=betap, epsa=epsa, epsb=epsb,
        ...                epsg=epsg, epst=epst, epss=epss, epsp=epsp,
        ...                steady=True, fullmodel=True, julian=False,
        ...                starch_mol2g=1., V0starchg=V0starch)
        >>> print(astr(Rass[0:6],5,pp=True))
        ['0.01095' '0.01092' '0.01087' '0.01083' '0.01083' '0.01085']
        >>> print(astr(Rm[0:6],5,pp=True))
        ['0.01107' '0.01123' '0.01119' '0.01116' '0.01116' '0.01118']
        >>> print(astr(Rchl[0:6],5,pp=True))
        ['0.01095' '0.01092' '0.01087' '0.01084' '0.01084' '0.01086']
        >>> print(astr(Rcyt[0:6],5,pp=True))
        ['0.01095' '0.01092' '0.01087' '0.01084' '0.01084' '0.01086']
        >>> print(astr(Rstarch[0:6],5,pp=True))
        ['0.01095' '0.01095' '0.01094' '0.01093' '0.01092' '0.01092']
        >>> print(astr(Rpyr[0:6],5,pp=True))
        ['0.01095' '0.01105' '0.01100' '0.01097' '0.01097' '0.01099']
        >>> print(astr(Rbio[0:6],5,pp=True))
        ['0.01095' '0.01094' '0.01089' '0.01086' '0.01086' '0.01088']
        >>> print(astr(Rphloem[0:6],5,pp=True))
        ['0.01095' '0.01091' '0.01086' '0.01082' '0.01082' '0.01084']
        >>> print(astr(Vstarch[0:6],5,pp=True))
        [' 40532.56148' ' 70158.59647' ' 98220.74429' '124787.71693' '150065.87342' '174445.40512']
        >>> print(astr(ass13[0:6],5,pp=True))
        ['-0.00560' ' 0.23301' ' 0.21418' ' 0.19580' ' 0.18136' ' 0.17631']
        >>> print(astr(Rnew_starch[0:6],5,pp=True))
        ['0.01095' '0.01096' '0.01091' '0.01088' '0.01088' '0.01090']
        >>> print(astr(Rnew_cyt[0:6],5,pp=True))
        ['0.01095' '0.01092' '0.01087' '0.01084' '0.01084' '0.01086']
        >>> print(astr(Vstarchg[0:6],5,pp=True))
        [' 40532.56148' ' 70158.59647' ' 98220.74429' '124787.71693' '150065.87342' '174445.40512']

        # non-steady state
        >>> R0cyt = 0.010911449304
        >>> Vcyt = np.array([135000.,135000.,135000.,135000.,135000.,135000.,135000.,
        ...                   135000.,135000.,135000.,135000.,135000.,135000.,135000.,
        ...                   135000.,135000.,135000.,135000.,135000.,135000.,135000.,
        ...                   135000.,135000.,135000.,135000.,135000.])
        >>> [Rass,Rm,Rchl,Rcyt,Rstarch,Rpyr,Rbio,Rphloem,Vstarch,ass13,disc,Rnew_starch,Rnew_cyt] = cuntz_gleixner(
        ...                adecdate[1:], gpp, Rd, CO2air, Ra, gtot, ndecdate, Vcyt=Vcyt,
        ...                date0=adecdate[0], V0starch=V0starch, R0starch=R0starch, R0cyt=R0cyt,
        ...                daynight=daynight,
        ...                daylength=daylength, Phi=Phi, s_resid=s_resid,
        ...                betas=betas, betap=betap, epsa=epsa, epsb=epsb,
        ...                epsg=epsg, epst=epst, epss=epss, epsp=epsp,
        ...                steady=False, fullmodel=True, julian=False)
        >>> print(astr(Vstarch[0:6],5,pp=True))
        [' 40532.56148' ' 70158.59647' ' 98220.74429' '124787.71693' '150065.87342' '174445.40512']
        >>> print(astr(ass13[0:6],5,pp=True))
        ['-0.00559' ' 0.23302' ' 0.21413' ' 0.19571' ' 0.18129' ' 0.17629']
        >>> print(astr(Rass[0:6],5,pp=True))
        ['0.01092' '0.01092' '0.01087' '0.01083' '0.01083' '0.01085']
        >>> print(astr(Rm[0:6],5,pp=True))
        ['0.01107' '0.01123' '0.01119' '0.01116' '0.01116' '0.01118']
        >>> print(astr(Rchl[0:6],5,pp=True))
        ['0.01095' '0.01092' '0.01087' '0.01084' '0.01084' '0.01086']
        >>> print(astr(Rstarch[0:6],5,pp=True))
        ['0.01095' '0.01095' '0.01094' '0.01093' '0.01092' '0.01092']
        >>> print(astr(Rcyt[0:6],5,pp=True))
        ['0.01092' '0.01092' '0.01090' '0.01088' '0.01087' '0.01087']
        >>> print(astr(Rpyr[0:6],5,pp=True))
        ['0.01092' '0.01105' '0.01103' '0.01101' '0.01100' '0.01099']
        >>> print(astr(Rbio[0:6],5,pp=True))
        ['0.01092' '0.01094' '0.01092' '0.01090' '0.01089' '0.01088']
        >>> print(astr(Rnew_cyt[0:6],5,pp=True))
        ['0.01095' '0.01092' '0.01087' '0.01084' '0.01084' '0.01086']
        >>> print(astr(Rnew_starch[0:6],5,pp=True))
        ['0.01095' '0.01096' '0.01091' '0.01088' '0.01088' '0.01090']
        >>> print(astr(Rphloem[0:6],5,pp=True))
        ['0.01092' '0.01091' '0.01089' '0.01087' '0.01085' '0.01085']
        >>> from dec2date import dec2date
        >>> from date2dec import date2dec
        >>> aa = dec2date(adecdate, ascii=True, calendar='decimal')
        >>> jadecdate = date2dec(ascii=aa)
        >>> ndecdate = 2008.918772768670806
        >>> bb = dec2date(ndecdate, ascii=True, calendar='decimal')
        >>> jndecdate = date2dec(ascii=bb)
        >>> [Rass,Rm,Rchl,Rcyt,Rstarch,Rpyr,Rbio,Rphloem,Vstarch,ass13,disc,Rnew_starch,Rnew_cyt] = cuntz_gleixner(
        ...                jadecdate[1:], gpp, Rd, CO2air, Ra, gtot, jndecdate, Vcyt=Vcyt,
        ...                date0=jadecdate[0], V0starch=V0starch, R0starch=R0starch, R0cyt=R0cyt,
        ...                daynight=daynight,
        ...                daylength=daylength, Phi=Phi, s_resid=s_resid,
        ...                betas=betas, betap=betap, epsa=epsa, epsb=epsb,
        ...                epsg=epsg, epst=epst, epss=epss, epsp=epsp,
        ...                steady=False, fullmodel=True, julian=True)
        >>> # There are slight differences due to precision of dates
        >>> print(astr(Rphloem[0:6],5,pp=True))
        ['0.01092' '0.01091' '0.01089' '0.01087' '0.01085' '0.01085']
        >>> [Rass,Rm,Rchl,Rcyt,Rstarch,Rpyr,Rbio,Rphloem,Vstarch,ass13,disc,Rnew_starch,Rnew_cyt] = cuntz_gleixner(
        ...                jadecdate[1:], gpp, Rd, CO2air, Ra, gtot, jndecdate, Vcyt=Vcyt,
        ...                date0=jadecdate[0], V0starch=V0starch, R0starch=R0starch, R0cyt=R0cyt,
        ...                daynight=daynight,
        ...                daylength=daylength, Phi=Phi, s_resid=s_resid,
        ...                betas=betas, betap=2./3.-0.1, epsa=epsa, epsb=epsb,
        ...                epsg=epsg, epst=epst, epss=epss, epsp=epsp,
        ...                steady=False, fullmodel=True, julian=True, nocheck=True)
        >>> print(astr(Rphloem[0:6],5,pp=True))
        ['0.01092' '0.00956' '0.00932' '0.00906' '0.00889' '0.00909']


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2012-2016 Matthias Cuntz - mc (at) macu (dot) de

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
        Modified, MC, Mar 2012 - julian
                  MC, May 2012 - nocheck
                  MC, May 2012 - Vcyt, daynight and betas=None default
                  MC, Feb 2013 - starch_mol2g, V0starchg
                  MC, Feb 2013 - ported to Python 3
                  MC, Nov 2016 - const.tiny -> const.eps
    """
    #
    # Checks
    nss = False
    if not steady: nss = True
    nd = idecdate.size
    ng = iGPP.size
    nr = iRd.size
    nc = iCa.size
    na = iRa.size
    ns = igtot.size
    if ((nd != ng) | (nd != nr) | (nd != nc) | (nd != na) | (nd != ns)):
        raise ValueError('not all input sizes are equal')
    if (Rass | Rm | Rchl | Rcyt | Rstarch | Rpyr | Rbio | Rphloem |
        Vstarch | ass13 | disc | Rnew_starch | Rnew_cyt): fullmodel = False
    if fullmodel:
        Rass        = True
        Rm          = True
        Rchl        = True
        Rcyt        = True
        Rstarch     = True
        Rpyr        = True
        Rbio        = True
        Rphloem     = True
        Vstarch     = True
        ass13       = True
        disc        = True
        Rnew_starch = True
        Rnew_cyt    = True
    # Defaults
    # Day (1) or night (0)
    if not np.any(daynight is not None): daynight = np.where(iGPP > 0., 1, 0)
    isarr = np.ndim(daynight)
    if (isarr==0):
        idaynight = np.ones(nd, dtype=np.int) * daynight
    else:
        idaynight = daynight
        nn = idaynight.size
        if (nn != nd):
            raise ValueError('daynight must be size 1 or size(idecdate)')
    # length of day [s]
    isarr = np.ndim(daylength)
    if (isarr==0):
        idaylength = np.ones(nd) * daylength
    else:
        idaylength = daylength
        nn = idaylength.size
        if (nn != nd):
            raise ValueError('daylength must be size 1 or size(idecdate)')
    # Vo/Vc: ratio of carboxylateion to oxygentation of Rubisco
    isarr = np.ndim(Phi)
    if (isarr==0):
        iPhi = np.ones(nd) * Phi
    else:
        iPhi = Phi
        nn = iPhi.size
        if (nn != nd):
            raise ValueError('Phi=Vo/Vc must be size 1 or size(idecdate)')
    # fraction of respiration occuring during biosynthesis: min=2/3; max=4/5
    isarr = np.ndim(betap)
    if (isarr==0):
        ibetap = np.ones(nd) * betap
    else:
        ibetap = betap
        nn = ibetap.size
        if (nn != nd):
            raise ValueError('betap must be size 1 or size(idecdate)')
    mini = np.amin(ibetap)
    maxi = np.amax(ibetap)
    if not nocheck:
        if ((mini < 0.66) | (maxi > 0.8)):
            raise ValueError('betap must be: 2/3 < betap < 4/5')
    # betas the factor of leaf respiration that is transferred to biosynthesis (default: 3*gpp/max(gpp))
    # betas*(1-betap) <= 1: if betap=2/3 -> betas<3: if betap=5/6 -> betas < 6
    if not np.any(betas is not None): betas = np.maximum((1./(1.-ibetap) * iGPP/np.amax(iGPP)) - const.eps, 0.)
    isarr = np.ndim(betas)
    if (isarr==0):
        ibetas = np.ones(nd) * betas
    else:
        ibetas = betas
        nn = ibetas.size
        if (nn != nd):
            raise ValueError('betas must be size 1 or size(idecdate)')
    # effective fractionation through leaf boundary layer, stomata and mesophyll
    isarr = np.ndim(epsa)
    if (isarr==0):
        iepsa = np.ones(nd) * epsa
    else:
        iepsa = epsa
        nn = iepsa.size
        if (nn != nd):
            raise ValueError('epsa must be size 1 or size(idecdate)')
    # Vcyt
    if (nss & (~np.any(Vcyt is not None))):
        raise ValueError('Vcyt must be given if non-steady state')
    if nss:
        isarr = np.ndim(Vcyt)
        if (isarr==0):
            iVcyt = np.ones(nd, dtype=np.float) * Vcyt
        else:
            iVcyt = Vcyt
            nn = iVcyt.size
            if (nn != nd):
                raise ValueError('Vcyt must be size 1 or size(idecdate)')
    # starch_mol2g
    if np.any(starch_mol2g is not None):
        isarr = np.ndim(starch_mol2g)
        if (isarr==0):
            istarch_mol2g = np.ones(nd, dtype=np.float) * starch_mol2g
        else:
            istarch_mol2g = starch_mol2g
            nn = istarch_mol2g.size
            if (nn != nd):
                raise ValueError('starch_mol2g must be size 1 or size(idecdate)')
    else:
        istarch_mol2g  = np.ones(nd, dtype=np.float)
    #
    # derived variables
    iAss = np.abs(iGPP) - np.abs(iRd)  # Net assimilation
    iCc  = iCa - iAss/igtot             # CO2 @ chloroplast
    iVc  = np.abs(iGPP)/(1.-0.5*iPhi) # carboxylation
    # seconds
    if julian:
        decfac = np.ones(nd)
    else:
        yr     = np.floor(idecdate)
        leap   = ((np.fmod(yr,4) == 0) & (np.fmod(yr,100) != 0)) | (np.fmod(yr,400) == 0)
        decfac = 365. + np.array(leap, dtype=np.float)
    isecs = np.rint(np.maximum(idecdate-np.roll(idecdate,1), 0.) * 86400. * decfac)
    if np.size(date0) == 0:
        isecs[0] = np.rint((idecdate[1]-idecdate[0])*86400.*decfac[0])
    else:
        isecs[0] = np.rint((idecdate[0]-date0)*86400.*decfac[0])
    # seconds remainig until next sunrise
    nightlength     = 86400. - idaylength # length of night before
    nightlength     = np.roll(nightlength,-1)
    nightlength[-1] = nightlength[-2]
    dsecs = idaylength + nightlength
    nsecs = np.empty(nd)
    issun = sunrise
    sr    = np.rint((idecdate - issun) * 86400. * decfac)
    for i in range(nd):
        if sr[i] > 0.: sr[i:] -= dsecs[i]
        nsecs[i] = -sr[i]
    nsecs = np.roll(nsecs,1)
    if np.size(date0) == 0:
        nsecs[0] = nsecs[0] + (nsecs[0]-nsecs[1])
    else:
        sr0 = np.rint((date0-issun)*86400.*decfac[0])
        if sr0 > 0.: sr0 -= dsecs[0]
        nsecs[0] = -sr0
    #
    # fraction of starch synthesis, same as T of Tchekez but with different units
    tmp    = 1./86400.
    ibigT  = 1. - idaylength*tmp # =0.056*6 @ Tcherkez
    ibetar = 1. - ibetas*(1.-ibetap)
    if not nocheck:
        if np.any(ibetar <= 0.):
            raise ValueError('betar = 1-betas*(1-betap) <= 0.')
    iepsr  = np.where(ibetar > const.eps, epss * (1.-ibetar)/ibetar, 1.) # limit to eps=1000 permil
    iepsg1 = epsg*0.5*iPhi + ibigT*(1.-0.5*iPhi)*epst
    iepsp1 = (epss*ibetas + iepsr*ibetar) / (ibetas + ibetar)
    #
    # Define output
    iAss13       = np.empty(nd)      # 13CO2 flux
    iRass        = np.empty(nd)      # Input
    idisc        = np.empty(nd)      # Discrimination
    iRm          = np.empty(nd)      # at site of carboxylation
    iRchl        = np.empty(nd)      # sucrose in chloroplast
    iRstarch     = np.empty(nd)      # accumulated starch in chloroplast
    iVstarch     = np.empty(nd)      # Starch concentration
    iVstarchg    = np.empty(nd)      # Starch concentration in gC/gDW
    iRcyt        = np.empty(nd)      # sucrose in cytoplast
    iRpyr        = np.empty(nd)      # pyruvate in cytoplast
    iRbio        = np.empty(nd)      # biosynthesis in cytoplast
    iRnew_starch = np.empty(nd)      # new starch in chloroplast
    iRnew_cyt    = np.empty(nd)      # new sucrose in cytoplasm
    iRphloem     = np.empty(nd)      # new phloem
    #
    if nss:
        dVcytdt = (iVcyt-np.roll(iVcyt,1))/isecs
        dVcytdt[0] = dVcytdt[1]
    #
    # Calc model
    if nss: # non-steady-state
        for i in range(nd):
            if idaynight[i] == 1: # nss day
                tmp1     = (1.-ibigT[i]) * (1.-0.5*iPhi[i]) * iVc[i]
                tmp2     = (1.-epsb) * iVc[i] * (1.-(1.-epsg)/(1.-iepsg1[i])*0.5*iPhi[i])
                tmp3     = (1.-epsb) / (1.-iepsg1[i])
                zaehler  = (1.-iepsa[i]) * igtot[i] * iCa[i] * iRa[i]
                nenner   = tmp2 + (1.-iepsa[i]) * igtot[i] * iCc[i]
                k1       = tmp1 * tmp3 * zaehler/nenner / iVcyt[i]
                k2       = ( (1. - tmp3*iRd[i]/(1.-iepsp1[i])/nenner) * tmp1 + dVcytdt[i] ) / iVcyt[i]
                Rcytm1   = R0cyt if i == 0 else iRcyt[i-1]
                iRcyt[i]        = k1/k2 + (Rcytm1-k1/k2) * np.exp(-k2*isecs[i])
                iRm[i]          = zaehler/nenner + iRd[i]/(1.-iepsp1[i])/nenner * iRcyt[i]
                iRchl[i]        = (1.-epsb) / (1.-iepsg1[i])*iRm[i]
                iRpyr[i]        = 1./(1.-iepsp1[i]) * iRcyt[i]
                iRbio[i]        = (1.-epss) / (1.-iepsp1[i]) * iRcyt[i]
                iRnew_starch[i] = (1.-epst) * iRchl[i]
                iRnew_cyt[i]    = iRchl[i]
                fphloem1        = (1.-ibigT[i]) * np.abs(iGPP[i]) - (ibetar[i]+ibetas[i]) * iRd[i]
                iRphloem[i]     = fphloem1 * iRcyt[i]
                fphloem2        = ibetap[i] * ibetas[i] * iRd[i]
                iRphloem[i]    += fphloem2 * (1.-epsp) * iRbio[i]
                iRphloem[i]    /= fphloem1 + fphloem2
                # Integrate starch
                fstarch      = ibigT[i] * np.abs(iGPP[i])
                Vstarchm1    = V0starch  if i == 0 else iVstarch[i-1]
                Vstarchgm1   = V0starchg if i == 0 else iVstarchg[i-1]
                Rstarchm1    = R0starch  if i == 0 else iRstarch[i-1]
                iVstarch[i]  = Vstarchm1  + fstarch*isecs[i]
                iVstarchg[i] = Vstarchgm1 + fstarch*isecs[i]*istarch_mol2g[i]
                iRstarch[i]  = (Rstarchm1*Vstarchm1 + fstarch*isecs[i]*iRnew_starch[i]) / iVstarch[i]
            else: # nss night
                Vstarchm1       = V0starch  if i == 0 else iVstarch[i-1]
                Vstarchgm1      = V0starchg if i == 0 else iVstarchg[i-1]
                fstarch         = -(Vstarchm1-s_resid)/nsecs[i]
                iVstarch[i]     = Vstarchm1  + fstarch*isecs[i]
                iVstarchg[i]    = Vstarchgm1 + fstarch*isecs[i]*istarch_mol2g[i]
                Rstarchm1       = R0starch  if i == 0 else iRstarch[i-1]
                iRstarch[i]     = Rstarchm1
                iRnew_starch[i] = Rstarchm1
                iRchl[i]        = iRnew_starch[i]
                fs              = np.abs(fstarch)
                k1              = fs/iVcyt[i]*iRchl[i]
                k2              = (fs+np.abs(dVcytdt[i])) / iVcyt[i]
                Rcytm1          = R0cyt if i == 0 else iRcyt[i-1]
                iRcyt[i]        = k1/k2 + (Rcytm1-k1/k2) * np.exp(-k2*isecs[i])
                iRpyr[i]        = iRcyt[i]
                iRbio[i]        = iRpyr[i] if i == 0 else iRbio[i-1]
                iRm[i]          = iCa[i]/iCc[i]*iRa[i] + (1.-iCa[i]/iCc[i])*iRpyr[i]/(1.-iepsa[i])
                iRnew_cyt[i]    = iRchl[i]
                iRphloem[i]     = iRcyt[i]
    else: # steady-state
        for i in range(nd):
            if idaynight[i] == 1: # ss day
                tmp1 = (1.-iepsa[i]) * igtot[i] * iRa[i] * iCa[i]
                tmp2 = (1.-iepsa[i]) * igtot[i] * iCc[i]
                tmp3 = ( (1.-epsb) * iVc[i] * (1.-(1.-epsg)/(1.-iepsg1[i])*0.5*iPhi[i]) -
                         (1.-epsb) / ((1.-iepsp1[i])*(1.-iepsg1[i])) * iRd[i] )
                iRm[i]          = tmp1 / (tmp2+tmp3)
                iRchl[i]        = (1.-epsb)/(1.-iepsg1[i]) * iRm[i]
                iRcyt[i]        = (1.-epsb)/(1.-iepsg1[i]) * iRm[i]
                iRpyr[i]        = (1.-epsb)/(1.-iepsg1[i])/(1.-iepsp1[i]) * iRm[i]
                iRbio[i]        = (1.-epsb)/(1.-iepsg1[i])/(1.-iepsp1[i])*(1.-epss) * iRm[i]
                iRnew_starch[i] = (1.-epst) * iRchl[i]
                iRnew_cyt[i]    = iRchl[i]
                fphloem1        = (1.-ibigT[i]) * np.abs(iGPP[i]) - (ibetar[i]+ibetas[i]) * iRd[i]
                iRphloem[i]     = fphloem1 * iRcyt[i]
                fphloem2        = ibetap[i] * ibetas[i] * iRd[i]
                iRphloem[i]    += fphloem2 * (1.-epsp) * iRbio[i]
                iRphloem[i]    /= fphloem1 + fphloem2
                # Integrate starch
                fstarch      = ibigT[i] * np.abs(iGPP[i])
                Vstarchm1    = V0starch  if i == 0 else iVstarch[i-1]
                Vstarchgm1   = V0starchg if i == 0 else iVstarchg[i-1]
                Rstarchm1    = R0starch  if i == 0 else iRstarch[i-1]
                iVstarch[i]  = Vstarchm1  + fstarch*isecs[i]
                iVstarchg[i] = Vstarchgm1 + fstarch*isecs[i]*istarch_mol2g[i]
                iRstarch[i]  = (Rstarchm1*Vstarchm1 + fstarch*isecs[i]*iRnew_starch[i]) / iVstarch[i]
            else: # ss night
                Vstarchm1       = V0starch  if i == 0 else iVstarch[i-1]
                Vstarchgm1      = V0starchg if i == 0 else iVstarchg[i-1]
                fstarch         = -(Vstarchm1-s_resid)/nsecs[i]
                iVstarch[i]     = Vstarchm1  + fstarch*isecs[i]
                iVstarchg[i]    = Vstarchgm1 + fstarch*isecs[i]*istarch_mol2g[i]
                Rstarchm1       = R0starch  if i == 0 else iRstarch[i-1]
                iRstarch[i]     = Rstarchm1
                iRnew_starch[i] = Rstarchm1
                iRchl[i]        = Rstarchm1
                iRcyt[i]        = iRchl[i]
                iRpyr[i]        = iRcyt[i]
                iRbio[i]        = iRpyr[i] if i == 0 else iRbio[i-1]
                iRm[i]          = iCa[i]/iCc[i]*iRa[i] + (1.-iCa[i]/iCc[i])*iRpyr[i]/(1.-iepsa[i])
                iRnew_cyt[i]    = iRchl[i]
                iRphloem[i]     = iRcyt[i]
    #
    iAss13 = (1.-iepsa) * igtot * (iRa*iCa-iRm*iCc)
    iRass  = iAss13 / iAss
    idisc  = 1. - iRass / iRa
    #
    out  = []
    if Rass:        out += [iRass]
    if Rm:          out += [iRm]
    if Rchl:        out += [iRchl]
    if Rcyt:        out += [iRcyt]
    if Rstarch:     out += [iRstarch]
    if Rpyr:        out += [iRpyr]
    if Rbio:        out += [iRbio]
    if Rphloem:     out += [iRphloem]
    if Vstarch:     out += [iVstarch]
    if ass13:       out += [iAss13]
    if disc:        out += [idisc]
    if Rnew_starch: out += [iRnew_starch]
    if Rnew_cyt:    out += [iRnew_cyt]
    if Vstarch & np.any(starch_mol2g is not None): out += [iVstarchg]

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


