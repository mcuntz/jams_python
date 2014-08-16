import numpy as np

def itc(H, zeta, ustar, varu, varw, vart, rho, lat, limit, outdir, plot=False):
    '''
    Calculation of integral turbulence characteristics after
    Thomas & Foken (2002). Plots itc's in outdir and returns flags if values
    deviate by limit from the itc model.
    
    
    Definition
    ----------
    itc(H, zeta, ustar, varu, varw, vart, rho, lat, limit, outdir, plot=False):

    
    Input
    ----- 
    H           np.array(N), sensible heat flux [W/m2]
    zeta        np.array(N), Monin-Obukhov stability parameter (z/L) [-]
    ustar       np.array(N), friction velocity  [m/s]
    varu        np.array(N), variance of u-component [m2/s2]
    varw        np.array(N), variance of w-component [m2/s2]
    vart        np.array(N), variance of tempreature [K2]
    rho         np.array(N), air density [kg/m3]
    lat         float, latitude of tower position [dec deg]
    limit       float, limit of exceedance from the model [-]
    outdir      str, path where the plots shall be saved
                  
                        
    Optional Input
    --------------
    plot        bool, if True plots are generated
    
    
    Output
    ------
    itcu          np.array(N), 0 where itcu did not exceed model limits, 1 where
                  model limit is exceeded, 2 where itcu exceeds stability limits
    itcw          np.array(N), 0 where itcw did not exceed model limits, 1 where
                  model limit is exceeded, 2 where itcw exceeds stability limits
    itct          np.array(N), 0 where itct did not exceed model limits, 1 where
                  model limit is exceeded, 2 where itct exceeds stability limits
    itc_u_pos.pdf plot with itcu of positive stability range
    itc_u_neg.pdf plot with itcu of negative stability range
    itc_w_pos.pdf plot with itcw of positive stability range
    itc_w_neg.pdf plot with itcw of negative stability range
    itc_t.pdf     plot with itct
    
    
    License
    -------
    This file is part of the UFZ Python library.

    The UFZ Python library is free software: you can redistribute it and/or 
    modify it under the terms of the GNU Lesser General Public License as 
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    The UFZ Python library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with The UFZ Python library.  If not,
    see <http://www.gnu.org/licenses/>.

    Copyright 2014 Arndt Piayda


    History
    -------
    Written,  AP, Aug 2014
    '''       
    
    
    ############################################################################
    # itc models
    def itc_flag(meas, mod, limit):
        return np.where(np.logical_or(meas>mod*(1.+limit), meas<mod*(1.-limit)), 1, 0)
    
    def u_pos(value, f):
        return 0.44*np.log((1. * f)/value)+6.3
    def u_neg(value):
        return 4.15*np.abs(value)**(1./8.)
    def w_pos(value, f):
        return 0.21*np.log((1. * f)/value)+3.1
    def w_neg(value):
        return 1.3*(1-2*value)**(1./3.)
    def t_all(value):
        return np.where(value>0.02, 1.4*np.abs(value)**(-1./4.),\
                   np.where((0.02>=value)&(value>-0.0625),\
                       0.5*np.abs(value)**(-1./2.),\
                           np.where((-0.0625>=value)&(value>-1),\
                               np.abs(value)**(-1./4.),\
                                   np.abs(value)**(-1./3.))))
        
    ############################################################################
    # calculate coriolis factor
    omega = 2.*np.pi/86164.1 # anglevelocity of earth, day of 23h56min4.1sec
    phi   = lat/360.*2.*np.pi
    f     = 2*omega*np.sin(phi) # coriolis factor

    ############################################################################
    # domains 
    zeta_neg = (-3.<zeta)  & (zeta<=-0.2)
    zeta_pos = (-0.2<zeta) & (zeta<0.4)
    
    ############################################################################
    # itcu
    itcu = np.where(np.isnan(varu), np.NaN, 0).astype(np.int)
    
    meas_u_pos      = np.sqrt(varu[zeta_pos])/ustar[zeta_pos]
    mod_u_pos       = u_pos(ustar[zeta_pos], f)
    itcu[zeta_pos] = itc_flag(meas_u_pos, mod_u_pos)
    
    meas_u_neg      = np.sqrt(varu[zeta_neg])/ ustar[zeta_neg]
    mod_u_neg       = u_neg(zeta[zeta_neg])
    itcu[zeta_neg] = itc_flag(meas_u_neg, mod_u_neg)
    
    itcu[(~np.isnan(itcu)) & (~zeta_neg) & (~zeta_pos)] = 2
    
    ############################################################################
    # itcw
    itcw = np.where(np.isnan(varw), np.NaN, 0).astype(np.int)
    
    meas_w_pos      = np.sqrt(varw[zeta_pos])/ustar[zeta_pos]
    mod_w_pos       = w_pos(ustar[zeta_pos], f)
    itcw[zeta_pos] = itc_flag(meas_w_pos, mod_w_pos)
    
    meas_w_neg      = np.sqrt(varw[zeta_neg])/ustar[zeta_neg]
    mod_w_neg       = w_neg(zeta[zeta_neg])
    itcw[zeta_neg] = itc_flag(meas_w_neg, mod_w_neg)
    
    itcw[(~np.isnan(itcw)) & (~zeta_neg) & (~zeta_pos)] = 2
    
    ############################################################################
    # itct
    itct = np.where(np.isnan(vart), np.NaN, 0).astype(np.int)
    
    meas_t = np.sqrt(vart) * 1004.67 * rho * ustar/np.abs(H)
    mod_t  = t_all(zeta)
    itct  = itc_flag(meas_t, mod_t)
    
    ############################################################################
    # plots
    if plot:
        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        
        x_pos     = np.log((1 * f)/ustar[zeta_pos])
        pos_range = np.arange(np.nanmin(ustar[zeta_pos]),
                              np.nanmax(ustar[zeta_pos]), 0.01)
        x_mod     = np.log((1 * f)/pos_range)

        x_neg     = zeta[zeta_neg]
        neg_range = np.arange(np.nanmin(x_neg),np.nanmax(x_neg), 0.01)
        
        ########################################################################
        # itcu: 0.44*ln((z_+*f)/u_star)+6.3 
        fig1 = plt.figure(1)
        sub1 = fig1.add_subplot(111)
        l1 = sub1.plot(x_pos, meas_u_pos, 'bo')
        l2 = sub1.plot(x_mod, u_pos(pos_range, f), 'g-',\
                       label='0.44*ln((z_+*f)/u_star)+6.3')
        l3 = sub1.plot(x_mod, u_pos(pos_range, f)*1.3, 'r-', label='+ 30%')
        l4 = sub1.plot(x_mod, u_pos(pos_range, f)*0.7, 'r-', label='- 30%')
        sub1.set_ylabel('sigma(u)/u_star')
        sub1.set_xlabel('ln(z_+*f/u_star)')
        sub1.legend(loc='best')
        plt.title('for -0.2 < zeta < 0.4')
        
        ########################################################################
        # itcu: 4.15*(np.abs(zeta))**(1/8) 
        fig2 = plt.figure(2)
        sub1 = fig2.add_subplot(111)
        l1 = sub1.plot(x_neg, meas_u_neg, 'bo')
        l2 = sub1.plot(neg_range, u_neg(neg_range), 'g-',\
                       label='4.15*(np.abs(zeta))**(1/8)')
        l3 = sub1.plot(neg_range, u_neg(neg_range)*1.3, 'r-', label='+30%')
        l4 = sub1.plot(neg_range, u_neg(neg_range)*0.7, 'r-', label='-30%')
        sub1.set_ylabel('sigma(u)/u_star')
        sub1.set_xlabel('zeta')
        sub1.legend(loc='best')
        plt.title('for -3 < zeta <= -0.2')
        
        ########################################################################
        # itcw: 0.21*ln((z_+*f)/u_star)+3.1
        fig3 = plt.figure(3)
        sub1 = fig3.add_subplot(111)
        l1 = sub1.plot(x_pos, meas_w_pos, 'bo')
        l2 = sub1.plot(x_mod, w_pos(pos_range, f), 'g-',\
                       label='0.21*ln((z_+*f)/u_star)+3.1')
        l3 = sub1.plot(x_mod, w_pos(pos_range, f)*1.3, 'r-', label='+ 30%')
        l4 = sub1.plot(x_mod, w_pos(pos_range, f)*0.7, 'r-', label='- 30%')
        sub1.set_ylabel('sigma(w)/u_star')
        sub1.set_xlabel('ln(z_+*f/u_star)')
        sub1.legend(loc='best')
        plt.title('for -0.2 < zeta < 0.4')
        
        ########################################################################
        # itcw: 1.3*(1-2*zeta)**(1/3) 
        fig4 = plt.figure(4)
        sub1 = fig4.add_subplot(111)
        l1 = sub1.plot(x_neg, meas_w_neg, 'bo')
        l2 = sub1.plot(neg_range, w_neg(neg_range), 'g-',\
                       label='1.3*(1-2*zeta)**(1/3)')
        l3 = sub1.plot(neg_range, w_neg(neg_range)*1.3, 'r-', label='+30%')
        l4 = sub1.plot(neg_range, w_neg(neg_range)*0.7, 'r-', label='-30%')
        sub1.set_ylabel('sigma(w)/u_star')
        sub1.set_xlabel('zeta')
        sub1.legend(loc='best')
        plt.title('for -3 < zeta <= -0.2')
        
        ########################################################################
        # itct: 1>zeta>0.02      => 1.4*np.abs(zeta)**(-1/4)
        #        0.02>zeta>-0.062 => 0.5*np.abs(zeta)**(-1/2)
        #        -0.062>zeta>-1   => 1.0*np.abs(zeta)**(-1/4)
        #        -1>zeta          => 1.0*np.abs(zeta)**(-1/3)
        fig5 = plt.figure(5)
        sub1 = fig5.add_subplot(111)
        x_mod = np.arange(-1,0.5,0.01)
        y_mod  = t_all(x_mod)
        l1 = sub1.plot(val[:,26], meas_t, 'bo')
        l2 = sub1.plot(x_mod, y_mod, 'g-', label='x*|zeta|**(-1/y)')
        l3 = sub1.plot(x_mod, y_mod*1.3, 'r-', label='+30%')
        l4 = sub1.plot(x_mod, y_mod*0.7, 'r-', label='-30%')
        sub1.set_ylabel('sigma(T)/T_star')
        sub1.set_xlabel('zeta')
        sub1.set_xlim(-1,0.5)
        sub1.set_ylim(0,10)
        sub1.legend(loc='best')
        plt.show()
        
        ########################################################################
        # save plots
        pp1  = pdf.PdfPages(outdir+'/itc_u_pos.pdf')
        pp2  = pdf.PdfPages(outdir+'/itc_u_neg.pdf')
        pp3  = pdf.PdfPages(outdir+'/itc_w_pos.pdf')
        pp4  = pdf.PdfPages(outdir+'/itc_w_neg.pdf')
        pp5  = pdf.PdfPages(outdir+'/itc_t.pdf')
        
        fig1.savefig(pp1, format='pdf')
        fig2.savefig(pp2, format='pdf')
        fig3.savefig(pp3, format='pdf')
        fig4.savefig(pp4, format='pdf')
        fig5.savefig(pp5, format='pdf')
        
        pp1.close()
        pp2.close()
        pp3.close()
        pp4.close()
        pp5.close()

    return itcu, itcw, itct