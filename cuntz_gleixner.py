#!/usr/bin/env python
import numpy as np
import const # from ufz

def cuntz_gleixner(idecdate, iGPP, iRd, iCa, iRa, igtot, sunrise, Vcyt=None,
                   date0=False,
                   V0starch=const.tiny,
                   R0starch=const.RPDB,
                   R0cyt=const.RPDB,
                   daynight=None, daylength=57600,
                   Phi=0.3, s_resid=const.tiny,
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
                   fullmodel=True, julian=True, nocheck=False):
    """
       Calculates the Cuntz-Gleixner steady state and non-steady state models
       of 13C discrimiantion in the Calvin cycle.

       Definition
       ----------
       def cuntz_gleixner(idecdate, iGPP, iRd, iCa, iRa, igtot, sunrise, Vcyt=False,
                          date0=False,
                          V0starch=const.tiny, R0starch=const.RPDB,
                          R0cyt=const.RPDB,
                          daynight=False, daylength=57600,
                          Phi=0.3, s_resid=const.tiny,
                          betas=False, betap=0.75,
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
                          fullmodel=True, julian=True, nocheck=False):


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
       Rass, Rm, Rchl, Rcyt, Rstarch, Rpyr, Rbio, Rphloem, Vstarch, ass13,
       disc, Rnew_starch, Rnew_cyt    if fullmodel=True


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
       >>> adecdate = np.array([2008.918658925319050,2008.918772768671033,2008.918886612022106, \
                                2008.919000455374089,2008.919114298724935,2008.919228142076918, \
                                2008.919341985427991,2008.919455828779974,2008.919569672131956, \
                                2008.919683515483030,2008.919797358835012,2008.919911202186086, \
                                2008.920025045538068,2008.920138888888914,2008.920252732240897, \
                                2008.920366575591970,2008.920480418943953,2008.920556314511941, \
                                2008.920594262295026,2008.920708105647009,2008.920821948998992, \
                                2008.920935792350065,2008.921049635702047,2008.921163479052893, \
                                2008.921277322405103,2008.921391165755949,2008.921505009107932])
       >>> gpp = np.array([0.000000000000,23.700827991217,22.449718259243,21.253578109071,20.222525197027, \
                             19.503625355216,18.797132965271,18.102416224453,17.780887860470,17.491607940331, \
                             17.207072197663,17.089915139494,17.995854647885,18.901914959729,19.681631460738, \
                             19.681631460738,19.681631460738,0.000000000000,0.000000000000,0.000000000000, \
                             0.000000000000,0.000000000000,0.000000000000,0.000000000000,0.000000000000, \
                             0.000000000000] )
       >>> Rd = np.array([0.511900000000,2.361686144687,2.743373026721,3.180029474251,3.476842651940, \
                            3.259038512076,3.053641828083,2.860020793216,2.958750580931,3.083603451827, \
                            3.213200496886,3.331826587704,3.352936200975,3.374166608865,3.392531460738, \
                            3.392531460738,3.392531460738,1.025929405070,0.829676977143,0.633424549217, \
                            0.437172122858,0.303515021488,0.547877741613,0.792240464668,1.036603184794, \
                            1.280965907360] )
       >>> CO2air = np.array([620.902600000000,537.510500000000,608.806500000000,671.251000000000, \
                                652.204000000000,560.157800000000,427.130100000000,395.276000000000, \
                                427.000400000000,410.953300000000,386.943500000000,500.417500000000, \
                                552.776800000000,515.865800000000,542.450400000000,692.503500000000, \
                                656.423500000000,588.844100000000,675.156500000000,725.101900000000, \
                                664.837000000000,598.080600000000,610.713600000000,487.087000000000, \
                                531.921300000000,675.177700000000] )
       >>> Ra = np.array([0.011067265443,0.011083081802,0.011071245659,0.011060401761,0.011063313320, \
                            0.011080216316,0.011111970396,0.011122420992,0.011111174802,0.011116914764, \
                            0.011125605614,0.011097923896,0.011079382516,0.011087211473,0.011083896499, \
                            0.011057329511,0.011062335683,0.011072518834,0.011061590657,0.011053508863, \
                            0.011061281634,0.011071628848,0.011069690431,0.011093962783,0.011086022577, \
                            0.011059558971] )
       >>> gtot = np.array([0.064395001124,0.074054920058,0.078085762302,0.078484156864,0.078127160737, \
                              0.085209848990,0.088685679784,0.089611189047,0.088528095110,0.086087621579, \
                              0.081901616151,0.076984314568,0.080693530135,0.084173028182,0.087005780100, \
                              0.087005780100,0.087005780100,0.046798889383,0.042324852911,0.037583815518, \
                              0.032460459750,0.028193059760,0.031985237181,0.035564641600,0.038983725824, \
                              0.042282334176] )
       >>> ndecdate = 2008.918772768670806
       >>> V0starch = 60498.901260546168
       >>> R0starch = 0.010949362493
       >>> daynight = np.array([0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0])
       >>> daylength = 56400.
       >>> Phi = np.array([0.081107124487,0.228537911761,0.178664767959,0.161131279216,0.173373369802, \
                             0.195820137312,0.276226485398,0.291675965181,0.259480203313,0.285111049806, \
                             0.331368464761,0.228636707479,0.199019858196,0.224041328002,0.209848472684, \
                             0.147530335706,0.158874834596,0.090595209673,0.077545142771,0.070693599170, \
                             0.075277039328,0.082213420279,0.081562178107,0.102842605886,0.095940191411, \
                             0.077688517811] )
       >>> s_resid = 40532.561476983901
       >>> betas = np.array([0.000000000000,4.999999000000,4.736061020188,4.483719593421,4.266205481166, \
                               4.114544146364,3.965500325995,3.818940732998,3.751110114569,3.690082727626, \
                               3.630056187040,3.605340371581,3.796460173119,3.987605459727,4.152096865111, \
                               4.152096865111,4.152096865111,0.000000000000,0.000000000000,0.000000000000, \
                               0.000000000000,0.000000000000,0.000000000000,0.000000000000,0.000000000000, \
                               0.000000000000] )
       >>> betap = 0.8
       >>> epsa = np.array([0.002995512907,0.003039740417,0.003192366495,0.003375544906,0.003479516543, \
                              0.003318658489,0.003187204366,0.003078822327,0.003144673535,0.003235857882, \
                              0.003343598868,0.003447509132,0.003408830551,0.003373557382,0.003345605523, \
                              0.003345605523,0.003345605523,0.003500242606,0.003551923813,0.003615033761, \
                              0.003693231875,0.003766948014,0.003706681514,0.003655727827,0.003612282052, \
                              0.003574980671] )
       >>> epsb = 0.029
       >>> epsg = 0.0185
       >>> epst = -0.004
       >>> epss = 0.01
       >>> epsp = 0.003
       >>> [Vstarch] = cuntz_gleixner(adecdate[1:], gpp, Rd, CO2air, Ra, gtot, ndecdate, \
                                      date0=adecdate[0], V0starch=V0starch, R0starch=R0starch, daynight=daynight, \
                                      daylength=daylength, Phi=Phi, s_resid=s_resid, \
                                      betas=betas, betap=betap, epsa=epsa, epsb=epsb, \
                                      epsg=epsg, epst=epst, epss=epss, epsp=epsp, \
                                      steady=True, Vstarch=True, julian=False)
       >>> print Vstarch
       [  40532.56147698   70158.59646601   98220.74429006  124787.7169264
         150065.87342268  174445.4051167   197941.82132329  220569.84160386
         242795.95142944  264660.46135486  286169.30160194  307531.6955263
         330026.51383616  353653.90753582  378255.94686174  402857.98618767
         419259.34573828  404110.27436783  358663.06025647  313215.84614512
         267768.63203376  222321.41792241  176874.20381105  131426.9896997
          85979.77558834   40532.56147698]

       >>> [Rass,Rm,Rchl,Rcyt,Rstarch,Rpyr,Rbio,Rphloem,Vstarch,ass13,disc,Rnew_starch,Rnew_cyt] = \
           cuntz_gleixner(adecdate[1:], gpp, Rd, CO2air, Ra, gtot, ndecdate, \
                          date0=adecdate[0], V0starch=V0starch, R0starch=R0starch, daynight=daynight, \
                          daylength=daylength, Phi=Phi, s_resid=s_resid, \
                          betas=betas, betap=betap, epsa=epsa, epsb=epsb, \
                          epsg=epsg, epst=epst, epss=epss, epsp=epsp, \
                          steady=True, fullmodel=True, julian=False)
       >>> print Rass
       [ 0.01094936  0.01092075  0.01086872  0.01083349  0.01083054  0.01085329
         0.01091213  0.01092848  0.01090373  0.01091399  0.0109335   0.01087484
         0.01084721  0.01086573  0.01085855  0.01080656  0.0108166   0.01092572
         0.01092572  0.01092572  0.01092572  0.01092572  0.01092572  0.01092572
         0.01092572  0.01092572]
       >>> print Rm
       [ 0.01106619  0.0112322   0.01118999  0.01115974  0.01115875  0.01117865
         0.01122928  0.01124342  0.01122279  0.01123199  0.01124851  0.01120095
         0.01117471  0.01119003  0.01118345  0.0111368   0.0111459   0.01106863
         0.01105886  0.01105151  0.01105939  0.01106978  0.01106687  0.01108836
         0.01108028  0.01105549]
       >>> print Rchl
       [ 0.01094936  0.01091611  0.0108697   0.01083842  0.01083878  0.01086053
         0.01091846  0.01093389  0.01091033  0.01092206  0.01094316  0.01088575
         0.01085705  0.01087464  0.01086671  0.01081468  0.01082473  0.01092572
         0.01092572  0.01092572  0.01092572  0.01092572  0.01092572  0.01092572
         0.01092572  0.01092572]
       >>> print Rcyt
       [ 0.01094936  0.01091611  0.0108697   0.01083842  0.01083878  0.01086053
         0.01091846  0.01093389  0.01091033  0.01092206  0.01094316  0.01088575
         0.01085705  0.01087464  0.01086671  0.01081468  0.01082473  0.01092572
         0.01092572  0.01092572  0.01092572  0.01092572  0.01092572  0.01092572
         0.01092572  0.01092572]
       >>> print Rstarch
       [ 0.01094936  0.01095376  0.01094216  0.01092931  0.01092136  0.01091893
         0.01092406  0.01092956  0.01093179  0.0109346   0.01093853  0.01093789
         0.01093534  0.01093419  0.01093263  0.01092807  0.01092572  0.01092572
         0.01092572  0.01092572  0.01092572  0.01092572  0.01092572  0.01092572
         0.01092572  0.01092572]
       >>> print Rpyr
       [ 0.01094936  0.01102637  0.01100025  0.01096706  0.01096599  0.01098694
         0.01104442  0.01105887  0.01103448  0.01104583  0.01106664  0.01100837
         0.01098097  0.01100027  0.01099346  0.01094082  0.01095099  0.01092572
         0.01092572  0.01092572  0.01092572  0.01092572  0.01092572  0.01092572
         0.01092572  0.01092572]
       >>> print Rbio
       [ 0.01094936  0.01091611  0.01089024  0.01085739  0.01085633  0.01087707
         0.01093397  0.01094828  0.01092414  0.01093537  0.01095598  0.01089829
         0.01087116  0.01089026  0.01088352  0.01083141  0.01084148  0.01084148
         0.01084148  0.01084148  0.01084148  0.01084148  0.01084148  0.01084148
         0.01084148  0.01084148]
       >>> print Rphloem
       [ 0.01094936  0.01089251  0.01085912  0.01082391  0.01082045  0.0108423
         0.01090029  0.01091589  0.01089086  0.01090076  0.01091983  0.01086099
         0.0108346   0.01085417  0.01084784  0.01079589  0.01080593  0.01092572
         0.01092572  0.01092572  0.01092572  0.01092572  0.01092572  0.01092572
         0.01092572  0.01092572]
       >>> print Vstarch
       [  40532.56147698   70158.59646601   98220.74429006  124787.7169264
         150065.87342268  174445.4051167   197941.82132329  220569.84160386
         242795.95142944  264660.46135486  286169.30160194  307531.6955263
         330026.51383616  353653.90753582  378255.94686174  402857.98618767
         419259.34573828  404110.27436783  358663.06025647  313215.84614512
         267768.63203376  222321.41792241  176874.20381105  131426.9896997
          85979.77558834   40532.56147698]
       >>> print ass13
       [-0.00560498  0.23303935  0.2141828   0.19579957  0.18136481  0.17630714
         0.17179506  0.1665762   0.16161657  0.15724885  0.15300197  0.14961704
         0.15883485  0.16872032  0.17687605  0.17602914  0.17619275 -0.01120902
        -0.00906482 -0.00692062 -0.00477642 -0.00331612 -0.00598596 -0.0086558
        -0.01132564 -0.01399547]
       >>> print Rnew_starch
       [ 0.01094936  0.01095977  0.01091318  0.01088177  0.01088213  0.01090398
         0.01096213  0.01097763  0.01095397  0.01096575  0.01098693  0.0109293
         0.01090048  0.01091814  0.01091018  0.01085794  0.01086803  0.01092572
         0.01092572  0.01092572  0.01092572  0.01092572  0.01092572  0.01092572
         0.01092572  0.01092572]
       >>> print Rnew_cyt
       [ 0.01094936  0.01091611  0.0108697   0.01083842  0.01083878  0.01086053
         0.01091846  0.01093389  0.01091033  0.01092206  0.01094316  0.01088575
         0.01085705  0.01087464  0.01086671  0.01081468  0.01082473  0.01092572
         0.01092572  0.01092572  0.01092572  0.01092572  0.01092572  0.01092572
         0.01092572  0.01092572]

       # non-steady state
       >>> R0cyt = 0.010911449304
       >>> Vcyt = np.array([135000.,135000.,135000.,135000.,135000.,135000.,135000., \
                             135000.,135000.,135000.,135000.,135000.,135000.,135000., \
                             135000.,135000.,135000.,135000.,135000.,135000.,135000., \
                             135000.,135000.,135000.,135000.,135000.])
       >>> [Rass,Rm,Rchl,Rcyt,Rstarch,Rpyr,Rbio,Rphloem,Vstarch,ass13,disc,Rnew_starch,Rnew_cyt] = \
           cuntz_gleixner(adecdate[1:], gpp, Rd, CO2air, Ra, gtot, ndecdate, Vcyt=Vcyt, \
                          date0=adecdate[0], V0starch=V0starch, R0starch=R0starch, R0cyt=R0cyt, \
                          daynight=daynight, \
                          daylength=daylength, Phi=Phi, s_resid=s_resid, \
                          betas=betas, betap=betap, epsa=epsa, epsb=epsb, \
                          epsg=epsg, epst=epst, epss=epss, epsp=epsp, \
                          steady=False, fullmodel=True, julian=False)
        >>> print Vstarch
        [  40532.56147698   70158.59646601   98220.74429006  124787.7169264
          150065.87342268  174445.4051167   197941.82132329  220569.84160386
          242795.95142944  264660.46135486  286169.30160194  307531.6955263
          330026.51383616  353653.90753582  378255.94686174  402857.98618767
          419259.34573828  404110.27436783  358663.06025647  313215.84614512
          267768.63203376  222321.41792241  176874.20381105  131426.9896997
           85979.77558834   40532.56147698]
        >>> print ass13
        [-0.00558824  0.23303895  0.21413257  0.19570929  0.1812934   0.17629181
          0.171858    0.16663691  0.1616368   0.1572795   0.15305138  0.14957457
          0.15875531  0.16869062  0.17684164  0.17591122  0.17611879 -0.01114552
         -0.00902824 -0.00690076 -0.00476668 -0.00331133 -0.00597985 -0.00864958
         -0.01131996 -0.01399062]
        >>> print Rass
        [ 0.01091666  0.01092073  0.01086617  0.01082849  0.01082628  0.01085234
          0.01091613  0.01093246  0.01090509  0.01091612  0.01093703  0.01087175
          0.01084178  0.01086382  0.01085644  0.01079932  0.01081206  0.01086382
          0.01088164  0.01089436  0.01090345  0.01090993  0.01091457  0.01091788
          0.01092024  0.01092193]
        >>> print Rm
        [ 0.01106578  0.01123222  0.0111918   0.01116236  0.01116084  0.01117914
          0.01122643  0.01124041  0.01122191  0.01123053  0.01124571  0.01120267
          0.01117738  0.0111911   0.01118457  0.0111395   0.01114772  0.0110664
          0.01105761  0.01105079  0.01105895  0.0110695   0.01106657  0.01108801
          0.01108001  0.01105533]
        >>> print Rchl
        [ 0.01094936  0.01091613  0.01087146  0.01084096  0.01084081  0.01086101
          0.01091568  0.01093096  0.01090947  0.01092064  0.01094044  0.01088743
          0.01085964  0.01087568  0.0108678   0.0108173   0.0108265   0.01092615
          0.01092615  0.01092615  0.01092615  0.01092615  0.01092615  0.01092615
          0.01092615  0.01092615]
        >>> print Rstarch
        [ 0.01094936  0.01095377  0.01094268  0.01093025  0.01092249  0.01091997
          0.01092464  0.01092978  0.01093191  0.01093459  0.01093832  0.01093781
          0.01093544  0.01093436  0.01093285  0.01092844  0.01092615  0.01092615
          0.01092615  0.01092615  0.01092615  0.01092615  0.01092615  0.01092615
          0.01092615  0.01092615]
        >>> print Rcyt
        [ 0.01091666  0.01091648  0.01090203  0.01088329  0.01087079  0.010868
          0.01088118  0.01089452  0.01089846  0.01090423  0.0109135   0.01090685
          0.01089426  0.0108891   0.01088297  0.01086405  0.01085642  0.01086382
          0.01088164  0.01089436  0.01090345  0.01090993  0.01091457  0.01091788
          0.01092024  0.01092193]
        >>> print Rpyr
        [ 0.01091666  0.01102675  0.01103297  0.01101246  0.01099838  0.01099449
          0.01100671  0.01101904  0.01102248  0.01102779  0.01103665  0.01102971
          0.0110186   0.01101489  0.01100991  0.01099076  0.01098305  0.01086382
          0.01088164  0.01089436  0.01090345  0.01090993  0.01091457  0.01091788
          0.01092024  0.01092193]
        >>> print Rbio
        [ 0.01091666  0.01091648  0.01092264  0.01090234  0.0108884   0.01088454
          0.01089665  0.01090885  0.01091225  0.01091751  0.01092628  0.01091941
          0.01090842  0.01090474  0.01089981  0.01088086  0.01087322  0.01087322
          0.01087322  0.01087322  0.01087322  0.01087322  0.01087322  0.01087322
          0.01087322  0.01087322]
        >>> print Rnew_cyt
        [ 0.01094936  0.01091613  0.01087146  0.01084096  0.01084081  0.01086101
          0.01091568  0.01093096  0.01090947  0.01092064  0.01094044  0.01088743
          0.01085964  0.01087568  0.0108678   0.0108173   0.0108265   0.01092615
          0.01092615  0.01092615  0.01092615  0.01092615  0.01092615  0.01092615
          0.01092615  0.01092615]
        >>> print Rnew_starch
        [ 0.01094936  0.0109598   0.01091494  0.01088433  0.01088417  0.01090445
          0.01095935  0.01097468  0.01095311  0.01096432  0.0109842   0.01093098
          0.01090308  0.01091918  0.01091127  0.01086057  0.01086981  0.01092615
          0.01092615  0.01092615  0.01092615  0.01092615  0.01092615  0.01092615
          0.01092615  0.01092615]
        >>> print Rphloem
        [ 0.01091666  0.01089289  0.01089142  0.01086872  0.01085241  0.01084975
          0.01086308  0.01087657  0.01087901  0.01088297  0.01089023  0.01088204
          0.01087173  0.0108686   0.01086407  0.01084518  0.01083756  0.01086382
          0.01088164  0.01089436  0.01090345  0.01090993  0.01091457  0.01091788
          0.01092024  0.01092193]
       >>> from dec2date import *
       >>> from date2dec import *
       >>> aa = dec2date(adecdate, ascii=True, calendar='decimal')
       >>> jadecdate = date2dec(ascii=aa)
       >>> ndecdate = 2008.918772768670806
       >>> bb = dec2date(ndecdate, ascii=True, calendar='decimal')
       >>> jndecdate = date2dec(ascii=bb)
       >>> [Rass,Rm,Rchl,Rcyt,Rstarch,Rpyr,Rbio,Rphloem,Vstarch,ass13,disc,Rnew_starch,Rnew_cyt] = \
           cuntz_gleixner(jadecdate[1:], gpp, Rd, CO2air, Ra, gtot, jndecdate, Vcyt=Vcyt, \
                          date0=jadecdate[0], V0starch=V0starch, R0starch=R0starch, R0cyt=R0cyt, \
                          daynight=daynight, \
                          daylength=daylength, Phi=Phi, s_resid=s_resid, \
                          betas=betas, betap=betap, epsa=epsa, epsb=epsb, \
                          epsg=epsg, epst=epst, epss=epss, epsp=epsp, \
                          steady=False, fullmodel=True, julian=True)
       >>> # There are slight differences due to precision of dates
       >>> print Rphloem
       [ 0.01091666  0.01089289  0.01089142  0.01086873  0.01085241  0.01084975
         0.01086308  0.01087658  0.01087901  0.01088297  0.01089023  0.01088204
         0.01087173  0.0108686   0.01086407  0.01084517  0.01083756  0.01086382
         0.01088164  0.01089436  0.01090345  0.01090993  0.01091457  0.01091788
         0.01092024  0.01092193]
       >>> [Rass,Rm,Rchl,Rcyt,Rstarch,Rpyr,Rbio,Rphloem,Vstarch,ass13,disc,Rnew_starch,Rnew_cyt] = \
           cuntz_gleixner(jadecdate[1:], gpp, Rd, CO2air, Ra, gtot, jndecdate, Vcyt=Vcyt, \
                          date0=jadecdate[0], V0starch=V0starch, R0starch=R0starch, R0cyt=R0cyt, \
                          daynight=daynight, \
                          daylength=daylength, Phi=Phi, s_resid=s_resid, \
                          betas=betas, betap=2./3.-0.1, epsa=epsa, epsb=epsb, \
                          epsg=epsg, epst=epst, epss=epss, epsp=epsp, \
                          steady=False, fullmodel=True, julian=True, nocheck=True)
       >>> print Rphloem
       [ 0.01091666  0.00956084  0.00932304  0.00905517  0.0088888   0.00908849
         0.00928323  0.00946748  0.00943847  0.00939182  0.00934281  0.00927611
         0.00915699  0.00905204  0.00896914  0.00896053  0.00895677  0.01072887
         0.01074818  0.01076197  0.01077181  0.01077884  0.01078386  0.01078745
         0.01079001  0.01079184]


       History
       -------
       Written,  MC, Jan 2012
       Modified, MC, Mar 2012 - julian
                 MC, May 2012 - nocheck
                 MC, May 2012 - Vcyt, daynight and betas=None default
    """
    #
    # Checks
    nss = False
    if not steady: nss = True
    nd = np.size(idecdate)
    ng = np.size(iGPP)
    nr = np.size(iRd)
    nc = np.size(iCa)
    na = np.size(iRa)
    ns = np.size(igtot)
    if ((nd != ng) | (nd != nr) | (nd != nc) | (nd != na) | (nd != ns)):
        print 'CUNTZ_GLEIXNER: not all input sizes are equal'
        return False
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
    if ~np.any(daynight != None): daynight = np.where(gpp > 0., 1, 0)
    isarr = np.size(np.shape(daynight))
    if (isarr==0):
        idaynight = np.ones(nd, dtype=np.int) * daynight
    else:
        idaynight = daynight
        nn = np.size(idaynight)
        if (nn != nd):
            print 'CUNTZ_GLEIXNER: daynight must be size 1 or size(idecdate)'
            return False
    # length of day [s]
    isarr = np.size(np.shape(daylength))
    if (isarr==0):
        idaylength = np.ones(nd) * daylength
    else:
        idaylength = daylength
        nn = np.size(idaylength)
        if (nn != nd):
            print 'CUNTZ_GLEIXNER: daylength must be size 1 or size(idecdate)'
            return False
    # Vo/Vc: ratio of carboxylateion to oxygentation of Rubisco
    isarr = np.size(np.shape(Phi))
    if (isarr==0):
        iPhi = np.ones(nd) * Phi
    else:
        iPhi = Phi
        nn = np.size(iPhi)
        if (nn != nd):
            print 'CUNTZ_GLEIXNER: Phi=Vo/Vc must be size 1 or size(idecdate)'
            return False
    # fraction of respiration occuring during biosynthesis: min=2/3; max=4/5
    isarr = np.size(np.shape(betap))
    if (isarr==0):
        ibetap = np.ones(nd) * betap
    else:
        ibetap = betap
        nn = np.size(ibetap)
        if (nn != nd):
            print 'CUNTZ_GLEIXNER: betap must be size 1 or size(idecdate)'
            return False
    mini = np.amin(ibetap)
    maxi = np.amax(ibetap)
    if not nocheck:
        if ((mini < 0.66) | (maxi > 0.8)):
            print 'CUNTZ_GLEIXNER: betap must be: 2/3 < betap < 4/5'
            return False
    # betas the factor of leaf respiration that is transferred to biosynthesis (default: 3*gpp/max(gpp))
    # betas*(1-betap) <= 1: if betap=2/3 -> betas<3: if betap=5/6 -> betas < 6
    if ~np.any(betas != None): betas = np.maximum((1./(1.-ibetap) * gpp/np.amax(gpp)) - const.tiny, 0.)
    isarr = np.size(np.shape(betas))
    if (isarr==0):
        ibetas = np.ones(nd) * betas
    else:
        ibetas = betas
        nn = np.size(ibetas)
        if (nn != nd):
            print 'CUNTZ_GLEIXNER: betas must be size 1 or size(idecdate)'
            return False
    # effective fractionation through leaf boundary layer, stomata and mesophyll
    isarr = np.size(np.shape(epsa))
    if (isarr==0):
        iepsa = np.ones(nd) * epsa
    else:
        iepsa = epsa
        nn = np.size(iepsa)
        if (nn != nd):
            print 'CUNTZ_GLEIXNER: epsa must be size 1 or size(idecdate)'
            return False
    # Vcyt
    if (nss & (~np.any(Vcyt != None))):
        print 'CUNTZ_GLEIXNER: Vcyt must be given if non-steady state'
        return False
    if nss:
        isarr = np.size(np.shape(Vcyt))
        if (isarr==0):
            iVcyt = np.ones(nd, dtype=np.int) * Vcyt
        else:
            iVcyt = Vcyt
            nn = np.size(iVcyt)
            if (nn != nd):
                print 'CUNTZ_GLEIXNER: Vcyt must be size 1 or size(idecdate)'
                return False
    #
    # derived variables
    iAss = np.abs(iGPP) - np.abs(iRd)  # Net assimilation
    iCc  = iCa - iAss/igtot             # CO2 @ chloroplast
    iVc   = np.abs(iGPP)/(1.-0.5*iPhi) # carboxylation
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
    for i in xrange(nd):
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
    # fraction of starch synthesis, same as T of Tchekez before but with different units,
    # also T of Tcherkez acts on A whereas bigT acts on GPP
    tmp    = 1./86400.
    ibigT  = 1. - idaylength*tmp # =0.058*6 @ Tcherkez
    ibetar = 1. - ibetas*(1.-ibetap)
    if not nocheck:
        if np.any(ibetar <= 0.):
            print 'CUNTZ_GLEIXNER: betar = 1-betas*(1-betap) <= 0.'
            return False
    iepsr  = np.where(ibetar > const.tiny, epss * (1.-ibetar)/ibetar, 1.) # limit to eps=1000 permil
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
        for i in xrange(nd):
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
                fphloem1       = (1.-ibigT[i]) * np.abs(iGPP[i]) - (ibetar[i]+ibetas[i]) * iRd[i]
                iRphloem[i]     = fphloem1 * iRcyt[i]
                fphloem2       = ibetap[i] * ibetas[i] * iRd[i]
                iRphloem[i]    += fphloem2 * (1.-epsp) * iRbio[i]
                iRphloem[i]    /= fphloem1 + fphloem2
                # Integrate starch
                fstarch    = ibigT[i] * np.abs(iGPP[i])
                Vstarchm1  = V0starch if i == 0 else iVstarch[i-1]
                Rstarchm1  = R0starch if i == 0 else iRstarch[i-1]
                iVstarch[i] = Vstarchm1 + fstarch*isecs[i]
                iRstarch[i] = (Rstarchm1*Vstarchm1 + fstarch*isecs[i]*iRnew_starch[i]) / iVstarch[i]
            else: # nss night
                Vstarchm1      = V0starch if i == 0 else iVstarch[i-1]
                fstarch        = -(Vstarchm1-s_resid)/nsecs[i]
                iVstarch[i]     = Vstarchm1 + fstarch*isecs[i]
                Rstarchm1      = R0starch if i == 0 else iRstarch[i-1]
                iRstarch[i]     = Rstarchm1
                iRnew_starch[i] = Rstarchm1
                iRchl[i]        = iRnew_starch[i]
                fs             = np.abs(fstarch)
                k1             = fs/iVcyt[i]*iRchl[i]
                k2             = (fs+np.abs(dVcytdt[i])) / iVcyt[i]
                Rcytm1         = R0cyt if i == 0 else iRcyt[i-1]
                iRcyt[i]        = k1/k2 + (Rcytm1-k1/k2) * np.exp(-k2*isecs[i])
                iRpyr[i]        = iRcyt[i]
                iRbio[i]        = iRpyr[i] if i == 0 else iRbio[i-1]
                iRm[i]          = iCa[i]/iCc[i]*iRa[i] + (1.-iCa[i]/iCc[i])*iRpyr[i]/(1.-iepsa[i])
                iRnew_cyt[i]    = iRchl[i]
                iRphloem[i]     = iRcyt[i]
    else: # steady-state
        for i in xrange(nd):
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
                fphloem1       = (1.-ibigT[i]) * np.abs(iGPP[i]) - (ibetar[i]+ibetas[i]) * iRd[i]
                iRphloem[i]     = fphloem1 * iRcyt[i]
                fphloem2       = ibetap[i] * ibetas[i] * iRd[i]
                iRphloem[i]    += fphloem2 * (1.-epsp) * iRbio[i]
                iRphloem[i]    /= fphloem1 + fphloem2
                # Integrate starch
                fstarch    = ibigT[i] * np.abs(iGPP[i])
                Vstarchm1  = V0starch if i == 0 else iVstarch[i-1]
                Rstarchm1  = R0starch if i == 0 else iRstarch[i-1]
                iVstarch[i] = Vstarchm1 + fstarch*isecs[i]
                iRstarch[i] = (Rstarchm1*Vstarchm1 + fstarch*isecs[i]*iRnew_starch[i]) / iVstarch[i]
            else: # ss night
                Vstarchm1      = V0starch if i == 0 else iVstarch[i-1]
                fstarch        = -(Vstarchm1-s_resid)/nsecs[i]
                iVstarch[i]     = Vstarchm1 + fstarch*isecs[i]
                Rstarchm1      = R0starch if i == 0 else iRstarch[i-1]
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

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod()

