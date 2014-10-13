#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from ufz.brewer         import get_brewer
from ufz.signature2plot import signature2plot
from ufz.abc2plot       import abc2plot
from ufz.colours        import colours
from ufz.position       import position

def clockplot(sub, si, sti, stierr=None,
              usetex      = False,
              dxabc       = 2,              # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
              dyabc       = 0.5,            # % of (max-min) shift up from lower x-axis for a,b,c,... labels
              dxsig       = 1.23,           # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
              dysig       = -0.05,          # % of (max-min) shift up from lower x-axis for a,b,c,... labels
              elwidth     = 1.0,            # errorbar line width
              alwidth     = 1.0,            # axis line width
              glwidth     = 0.5,            # grid line width
              mcol        = '0.4',          # grid colour
              mcols       = ['0.0', '0.4', '0.4', '0.7', '0.7', '1.0'],      # stack colors
              lcols       = ['None', 'None', 'None', 'None', 'None', '0.0'], # stack bourder colours
              hatches     = [None, None, None, None, None, '//'],            # stack hatching
              llxbbox     =  0.0,           # x-anchor legend bounding box
              llybbox     = 1.15,           # y-anchor legend bounding box
              ymax        = 0.8,            # maximum of y-axis
              ntextsize   = 'medium',       # normal textsize
              bmod        = 0.5,            # fraction of ymax from center to start module colours
              alphamod    = 0.7,            # alpha channel for modules
              fwm         = 0.05,           # module width to remove at sides
              ylabel1     = 1.15,           # position of module names
              ylabel2     = 1.35,           # position of class names
              mtextsize   = 'large',        # 1.3*textsize # textsize of module labels
              bpar        = 0.4,            # fraction of ymax from center to start with parameter bars
              fwb         = [0.7,0.4,0.3],  # fractional width of bars
              plwidth     = 0.5,            # stack border line width
              bplabel     = 0.1,            # fractional distance of ymax of param numbers in centre from 0-line
              ptextsize   = 'medium',       # 'small' # 0.8*textsize # textsize of param numbers in centre
              space4yaxis = 2,              # space for y-axis (integer)
              ytextsize   = 'medium',       # 'small' # 0.8*textsize # textsize of y-axis
              dobw        = False,          # True: black & white
              docomp      = True,           # True: Print classification on top of modules
              dosig       = False,          # True: add signature to plot
              dolegend    = False,          # True: add legend to each subplot
              doabc       = False,          # True: add subpanel numbering
              do1legend   = False,          # True: add one legend in next subplot
              mod         = ['Interception',    'Snow',             'Soil moisture', 'Soil moisture',
                             'Direct\n runoff', 'Evapo-\n transp.', 'Interflow',     'Percolation',
                             'Routing',         'Geology'],      # module names
              comp        = ['P',               'P',                'S',             'ET',
                             'Q',               'ET',               'S',             'S',
                             'Q',               'S'],            # class names
              pmod        = [ 1,                 8,                  9,               8,
                              1,                 3,                  5,               3,
                              5,                 9],             # # of parameters per class
              saname      = ['Sobol', 'weighted Sobol', 'RMSE'], # stack names
              sig         = 'J Mai & M Cuntz'):                  # signature
    """
        The clock plot with modules and up to three index stacks.

        The plot currently defaults to mHM but it can be customized for any model output.


        Definition
        ----------
        def clockplot(sub, si, sti, stierr=None,
                      usetex=False,
                      dxabc       = 2,
                      dyabc       = 0.5,
                      dxsig       = 1.23,
                      dysig       = -0.05,
                      elwidth     = 1.0,
                      alwidth     = 1.0,
                      glwidth     = 0.5,
                      mcol        = '0.4',
                      mcols       = ['0.0', '0.4', '0.4', '0.7', '0.7', '1.0'],
                      lcols       = ['None', 'None', 'None', 'None', 'None', '0.0'],
                      hatches     = [None, None, None, None, None, '//'],
                      llxbbox     =  0.0,
                      llybbox     = 1.15,
                      ymax = 0.8,
                      ntextsize   = 'medium',
                      bmod        = 0.5,
                      alphamod    = 0.7,
                      fwm         = 0.05,
                      ylabel1     = 1.15,
                      ylabel2     = 1.35,
                      mtextsize   = 'large',
                      bpar        = 0.4,
                      fwb         = [0.7,0.4,0.3],
                      plwidth     = 0.5,
                      bplabel     = 0.1,
                      ptextsize   = 'medium',
                      space4yaxis = 2,
                      ytextsize   = 'medium',
                      dobw      = False,
                      docomp    = True,
                      dosig     = False,
                      dolegend  = False,
                      doabc     = False,
                      do1legend = False,
                      mod   = ['Interception',    'Snow',             'Soil moisture', 'Soil moisture',
                               'Direct\n runoff', 'Evapo-\n transp.', 'Interflow',     'Percolation',
                               'Routing',         'Geology'],
                      comp  = ['P',               'P',                'S',             'ET',
                               'Q',               'ET',               'S',             'S',
                               'Q',               'S'],
                      pmod  = [ 1,                 8,                  9,               8,
                                1,                 3,                  5,               3,
                                5,                 9],
                      saname = ['Sobol', 'weighted Sobol', 'RMSE'],
                      sig = 'J Mai & M Cuntz'):


        Input
        -----
        sub                          axes handle from for example
                                     sub = fig.add_axes(ufz.position(nrow,ncol,iplot), polar=True)
        si                           list, list of arrays, 1D-, or 2D-array si[nstacks, nparameters] of
                                     first-order Sobol' indexes
        sti                          list, list of arrays, 1D-, or 2D-array si[nstacks, nparameters] of
                                     total-order Sobol' indexes


        Optional Input
        --------------
        stierr = None                list, list of arrays, 1D-, or 2D-array si[nstacks, nparameters] of
                                     error bars of total-order Sobol' indexes
        usetex = False               True: use LaTeX rendering
        dxabc = 2                    % of (max-min) shift to the right from left y-axis for a,b,c,... labels
        dyabc = 0.5                  % of (max-min) shift up from lower x-axis for a,b,c,... labels
        dxsig = 1.23                 % of (max-min) shift to the right from left y-axis for signature
        dysig = -0.05                % of (max-min) shift up from lower x-axis for signature
        elwidth = 1.0                errorbar line width
        alwidth = 1.0                axis line width
        glwidth = 0.5                grid line width
        mcol = '0.4'                 grid colour
        mcols = ['0.0', '0.4', '0.4', '0.7', '0.7', '1.0']        stack colors
        lcols = ['None', 'None', 'None', 'None', 'None', '0.0']   stack border colours
        hatches = [None, None, None, None, None, '//']            stack hatching
        llxbbox =  0.0               x-anchor legend bounding box
        llybbox = 1.15,              y-anchor legend bounding box
        ymax = 0.8                   y-axis maximum
        ntextsize = 'medium'         normal textsize
        bmod = 0.5                   fraction of ymax to start module colours
        alphamod = 0.7               alpha channel for modules
        fwm = 0.05                   module width to remove at sides for space between modules
        ylabel1 = 1.15               fractional position of module names
        ylabel2 = 1.35               fractional position of class names
        mtextsize = 'large'          textsize of module labels
        bpar = 0.4                   fraction of ymax to start parameter bars
        fwb = [0.7,0.4,0.3]          fractional width of bars depending on number of index stacks
        plwidth = 0.5                stack border line width
        bplabel = 0.1                fractional distance of ymax of param numbers in centre from 0-line
        ptextsize = 'medium'         textsize of param numbers in centre
        space4yaxis = 2              space for y-axis top tickmark label (integer)
        ytextsize = 'medium'         textsize of y-axis tickmark labels
        dobw = False                 True: black & white; False: colour
        docomp = True                True: Print classification on top of modules
        dosig = False                True: add signature (sig) to plot
        dolegend = False             True: add legend to each subplot
        doabc = False                True: add subpanel numbering
        do1legend = False            True: add one legend in next subplot
        mod  = ['Interception',    'Snow',             'Soil moisture', 'Soil moisture',
                'Direct\n runoff', 'Evapo-\n transp.', 'Interflow',     'Percolation',
                'Routing',         'Geology'],                                           module names
        comp = ['P',               'P',                'S',             'ET',
                'Q',               'ET',               'S',             'S',
                'Q',               'S'],                                                 class names
        pmod = [ 1,                 8,                  9,               8,
                 1,                 3,                  5,               3,
                 5,                 9],                                                  number of parameters per class
        saname = ['Sobol', 'weighted Sobol', 'RMSE'],                                    stack names
        sig = 'J Mai & M Cuntz'                                                          signature


        Output
        ------
        clockplot on axes sub


        Restrictions
        ------------
        None


        Examples
        --------
        see below if __name__ == '__main__':


        License
        -------
        This file is part of the UFZ Python package.

        The UFZ Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2014 Matthias Cuntz


        History
        -------
        Written,  MC, JM, AP, Oct 2014
    """
    # Check si[nstacks, nparams]
    assert np.shape(si) == np.shape(sti), 'Si and STi must have same dimensions.'
    if stierr is not None:
        assert np.shape(sti) == np.shape(stierr), 'STi their errors stierr must have same dimensions.'
    si_shape = np.shape(si)
    if np.size(si_shape) == 1:
        nsi = 1
        isi  = np.array(si)
        isti = np.array(sti)
        isi  = isi[np.newaxis,:]
        isti = isti[np.newaxis,:]
        if stierr is not None:
            istierr = np.array(stierr)
            istierr = istierr[np.newaxis,:]
    elif np.size(si_shape) == 2:
        nsi = si_shape[0]
        if nsi > 3:
            raise ValueError('first data dimension must be <= 3, i.e. at most 3 stacks per parameter supported.')
        isi  = np.array(si)
        isti = np.array(sti)
        if stierr is not None:
            istierr = np.array(stierr)
    else:
        raise ValueError('input data must be 1D or 2D.')
    npar = isi.shape[1]
    idsi = isti-isi

    isaname = saname[:]
    # Prepare annotations
    if usetex:
        imod = []
        for i in mod:
            if '\n' in i:
                ss = i.split('\n')
            else:
                ss = [i]
            ss = [ r''+j.strip() for j in ss ]
            for j, s in enumerate(ss):
                ss[j] = r'$\mathrm{'+s+'}$'
                if '-' in ss[j]:
                    ss[j] = ss[j].replace('-}$', '}$-')
                if ' ' in ss[j]:
                    ss[j] = ss[j].replace(' ', '\ ')
            imod.append(ss)
        mod = imod
        comp = [ r'$\mathbf{'+i+'}$' for i in comp ]
        for j, s in enumerate(isaname):
            isaname[j] = r'$\mathrm{'+s+'}$'
            if '-' in isaname[j]:
                isaname[j] = isaname[j].replace('-}$', '}$-')
            if ' ' in isaname[j]:
                isaname[j] = isaname[j].replace(' ', '\ ')
    else:
        imod = []
        for i in mod:
            if '\n' in i:
                ss = i.split('\n')
            else:
                ss = [i]
            ss = [ r''+j.strip() for j in ss ]
            imod.append(ss)
        mod = imod

    # numbers
    nmod   = len(mod)  # number modules
    nparam = sum(pmod) # number of parameters
    param  = np.arange(nparam) + 1
    assert nparam == npar, 'si.shape[1] must be sum(pmod).'

    # colours
    if dobw:
        c = np.linspace(0.2, 0.85, nmod)
        c = np.ones(nmod)*0.7
        c = [ str(i) for i in c ]
    else:
        # c = [(165./255.,  0./255., 38./255.), # interception
        #      (215./255., 48./255., 39./255.), # snow
        #      (244./255.,109./255., 67./255.), # soil moisture
        #      (244./255.,109./255., 67./255.), # soil moisture
        #      (253./255.,174./255., 97./255.), # direct runoff
        #      (254./255.,224./255.,144./255.), # Evapotranspiration
        #      (171./255.,217./255.,233./255.), # interflow
        #      (116./255.,173./255.,209./255.), # percolation
        #      ( 69./255.,117./255.,180./255.), # routing
        #      ( 49./255., 54./255.,149./255.)] # geology
        c = get_brewer('rdylbu11', rgb=True)
        tmp = c.pop(5)   # rm yellow
        c.insert(2,c[2]) # same colour for both soil moistures

    # conversion factor from parameter number to radian
    n2rad = 2.*np.pi/(nparam+space4yaxis)

    # -------------------------------------------------------------------------
    # Plot
    #

    # Axes have to be defined externally
    # sub = fig.add_axes(ufz.position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), polar=True)
    sub.set_theta_zero_location('N') # 0 is North
    sub.set_theta_direction(-1)      # clockwise

    xlim = [0, 2.*np.pi]
    ylim = [-bpar*ymax, ymax]

    # coloured modules
    mleft   = (space4yaxis+np.cumsum([0]+pmod[:-1])+fwm)*n2rad # left start at space4yaxis
    mheight = np.ones(nmod)*ymax*(1.-bmod)                     # height from bmod*ymax to ymax
    mwidth  = (np.array(pmod)-2.*fwm)*n2rad                    # width is number of params per module
    bar1    = sub.bar(mleft, mheight, mwidth, bottom=bmod*ymax,
                      color=c, alpha=alphamod,
                      linewidth=0)

    # module and class labels
    xm = mleft+0.5*mwidth
    for i in range(nmod):
        # module
        for j, m in enumerate(mod[i]):
            mlabel12 = (ylabel2-ylabel1)*0.3
            if len(mod[i]) > 1:
                ylab = ylabel1 + (2*j-1)*mlabel12
            else:
                ylab = ylabel1
            label = sub.text(xm[i], ylab*ymax, m,
                             fontsize=ntextsize, horizontalalignment='center', verticalalignment='center')
            if (xm[i] < 0.5*np.pi) | (xm[i] > 1.5*np.pi):
                label.set_rotation(np.rad2deg(-xm[i]))
            else:
                label.set_rotation(np.rad2deg(-xm[i])+180.)
        # class
        if docomp:
            label = sub.text(xm[i], ylabel2*ymax, comp[i],
                             fontsize=mtextsize, fontweight='bold',
                             horizontalalignment='center', verticalalignment='center')
            if (xm[i] < 0.5*np.pi) | (xm[i] > 1.5*np.pi):
                label.set_rotation(np.rad2deg(-xm[i]))
            else:
                label.set_rotation(np.rad2deg(-xm[i])+180.)

    # y-axis
    # grid
    nyticks = 5
    dyy     = np.linspace(0,ymax,nyticks)
    gy      = np.delete(dyy[1:-1], nyticks//2-1)
    npoints = 100
    gxx     = np.linspace(0,2.*np.pi,npoints)
    for i in range(gy.size):
        grid = sub.plot(gxx, np.ones(npoints)*gy[i], linestyle='--', color=mcol, linewidth=glwidth)
    # in "axis normal coordinates" for rectangular ticks, etc.
    dyy    = np.array([0,ymax])
    yy     = ((1.+2.*bpar)*ymax+dyy)/((2.+2.*bpar)*ymax)
    yaxis  = sub.plot([0.5,0.5], yy, transform=sub.transAxes, linestyle='-', linewidth=alwidth, color='k')
    nyticks = 5
    xx      = np.ones(nyticks)*0.5
    dyy     = np.linspace(0,ymax,nyticks)
    yy      = ((1.+2.*bpar)*ymax+dyy)/((2.+2.*bpar)*ymax)
    ytickwidth = 0.015
    for i in range(nyticks):
        yticks = sub.plot([xx[i],xx[i]+ytickwidth], [yy[i],yy[i]], transform=sub.transAxes,
                          linestyle='-', linewidth=alwidth, color='k')
    # y-tickmarks at top and bottom
    tx = xx[0]+1.5*ytickwidth
    ty = yy[0]
    if usetex:
        tt = r'$\mathrm{0}$'
    else:
        tt = '0'
    label = sub.text(tx, ty, tt, transform=sub.transAxes,
                     fontsize=ytextsize, horizontalalignment='left', verticalalignment='bottom')
    bx = xx[-1]+1.5*ytickwidth
    by = yy[-1]
    if usetex:
        bt = r'$\mathrm{'+str(ymax)+'}$'
    else:
        bt = str(ymax)
    label = sub.text(bx, by, bt, transform=sub.transAxes,
                     fontsize=ytextsize, horizontalalignment='left', verticalalignment='top')

    # param numbers in center == x-tickmarks
    xticknames = [ str(i) for i in param[1::4] ]
    if usetex:
        xticknames = [ r'$\mathrm{'+i+'}$' for i in xticknames ]
    shiftx  = np.floor(1./(nsi+1)*10.)/10.
    pleft   = (param+(space4yaxis-1)+shiftx-0.5*fwb[nsi-1])*n2rad # center at 0.5 from param number
    pwidth  = np.ones(nparam)*fwb[nsi-1]*n2rad
    tx = pleft[1::4]+0.5*pwidth[1::4]
    ty = -np.ones(tx.size)*bplabel*ymax
    for i in range(tx.size):
        if (tx[i] < np.pi):
            rot = np.rad2deg(-tx[i])+90.
        else:
            rot = np.rad2deg(-tx[i])-90.
        label = sub.text(tx[i], ty[i], xticknames[i], rotation=rot, fontsize=ptextsize,
                         horizontalalignment='center', verticalalignment='center')

    # Index stacks
    for n in range(nsi):
        # params on bottom
        pleft = (param+(space4yaxis-1)+shiftx-0.5*fwb[nsi-1]+n*fwb[nsi-1])*n2rad # center at 0.5 from param number
        pheight = isi[n,:]
        bar2    = sub.bar(pleft, pheight, pwidth, bottom=0.,
                          facecolor=mcols[2*n], hatch=hatches[2*n], edgecolor=lcols[2*n],
                          linewidth=plwidth)
        # params on top
        pheight = idsi[n,:]
        pbottom = isi[n,:]
        bar2    = sub.bar(pleft, pheight, pwidth, bottom=pbottom,
                          facecolor=mcols[2*n+1], hatch=hatches[2*n+1], edgecolor=lcols[2*n+1],
                          linewidth=plwidth)
        if stierr is not None:
            # error bars
            xx      = pleft+0.5*pwidth
            yy      = isti[n,:]
            perr    = istierr[n,:]
            pewidth = 0.1*n2rad
            for i in range(xx.size):
                yerrm = sub.plot([xx[i],xx[i]], [yy[i]-perr[i],yy[i]+perr[i]], # middle line
                                 linestyle='-', linewidth=elwidth, color='k')
                yerrl = sub.plot([xx[i]-pewidth,xx[i]+pewidth], [yy[i]-perr[i],yy[i]-perr[i]], # lower bar
                                 linestyle='-', linewidth=elwidth, color='k')
                yerru = sub.plot([xx[i]-pewidth,xx[i]+pewidth], [yy[i]+perr[i],yy[i]+perr[i]], # upper bar
                                 linestyle='-', linewidth=elwidth, color='k')

    if dosig:
        signature2plot(sub, dxsig, dysig, sig, transform=sub.transAxes,
                       italic=True, small=True, mathrm=True, usetex=usetex)

    # General settings
    sub.set_xlim(xlim)
    sub.set_ylim(ylim)     # start at -bpar to get offset for bars
    sub.grid(False)
    sub.set_frame_on(False)
    sub.set_xticks([])
    sub.set_xticklabels([])
    sub.set_yticks([])
    sub.set_yticklabels([])

    # Fake subplot for legend and numbering
    spos = sub.get_position()
    lsub = fig.add_axes([spos.x0+llxbbox*spos.width, spos.y0+llybbox*spos.height, 0.5*spos.width, 0.1*spos.height])

    if dolegend:
        x1, y1 = lsub.transData.transform_affine(np.array([0,0]))
        x2, y2 = lsub.transData.transform_affine(np.array([1,1]))
        dpi = lsub.figure.dpi                         # pixels per inch
        ss  = mpl.rcParams['font.size']               # text size in points: 1 pt = 1/72 inch
        shifty = 1.0 * ss/72.*float(dpi)/float(y2-y1) # shift by 1.0 of textsize
        for n in range(nsi):
            dy1 = nsi*shifty - (n+1.)   * shifty
            dy2 = nsi*shifty - (n+0.01) * shifty
            lsub.fill_between([0,0.3],   [dy1,dy1], [dy2,dy2], linewidth=plwidth,
                              facecolor=mcols[2*n],   edgecolor=lcols[2*n], hatch=hatches[2*n])
            lsub.fill_between([0.3,0.6], [dy1,dy1], [dy2,dy2], linewidth=plwidth,
                              facecolor=mcols[2*n+1], edgecolor=lcols[2*n+1], hatch=hatches[2*n+1])
            lsub.text(0.65, 0.5*(dy1+dy2), isaname[n], fontsize=ntextsize, horizontalalignment='left', verticalalignment='center')
            if n == 0:
                t = r'$S_i$'
                lsub.text(0.15, dy2+0.05, t, fontsize=ntextsize, horizontalalignment='center', verticalalignment='bottom')
                t = r'$S_{Ti}-S_i$'
                lsub.text(0.45, dy2+0.05, t, fontsize=ntextsize, horizontalalignment='center', verticalalignment='bottom')

    lsub.set_xlim([0,1])
    lsub.set_ylim([0,1])

    # subplot numbering
    if doabc:
        abc2plot(lsub, dxabc, dyabc, iplot, lower=True, parenthesis='close',
                 bold=True, large=True,
                 mathrm=True, usetex=usetex,
                 horizontalalignment='right', verticalalignment='bottom')

    lsub.set_title('')
    lsub.set_xlabel('')
    lsub.set_ylabel('')
    lsub.set_xticks([])
    lsub.set_yticks([])
    lsub.set_axis_off()


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # nn = 52
    # ymax = 0.7
    # si  = np.arange(nn)/np.float(nn-1) * ymax * 0.5
    # sti = np.arange(nn)/np.float(nn-1) * ymax
    # stierr = 0.1*sti

    # pngbase = 'clockplot_'
    # pdffile = '' # 'clockplot.pdf'
    # usetex  = True

    # if (pdffile != '') & (pngbase != ''):
    #     print('\nError: PDF and PNG are mutually exclusive. Only either -p or -g possible.\n')
    #     parser.print_usage()
    #     import sys
    #     sys.exit()

    # if (pdffile == ''):
    #     if (pngbase == ''):
    #         outtype = 'x'
    #     else:
    #         outtype = 'png'
    # else:
    #     outtype = 'pdf'

    # # Main plot
    # nrow        = 3           # # of rows of subplots per figure
    # ncol        = 2           # # of columns of subplots per figure
    # hspace      = 0.05        # x-space between subplots
    # vspace      = 0.09        # y-space between subplots
    # right       = 0.9         # right space on page
    # textsize    = 6           # standard text size
    # dxabc       = 2           # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
    # dyabc       = 0.5           # % of (max-min) shift up from lower x-axis for a,b,c,... labels
    # dxsig       = 1.23        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
    # dysig       = -0.05       # % of (max-min) shift up from lower x-axis for a,b,c,... labels

    # lwidth      = 1.5         # linewidth
    # elwidth     = 1.0         # errorbar line width
    # alwidth     = 1.0         # axis line width
    # glwidth     = 0.5         # grid line width
    # msize       = 1.0         # marker size
    # mwidth      = 1.0         # marker edge width
    # mcol1       = '0.0'       # primary marker colour
    # mcol2       = '0.4'                     # secondary
    # mcol3       = '0.0' # third
    # mcols       = ['0.0', '0.4', '0.4', '0.7', '0.7', '1.0']
    # lcol1       = colours('blue')   # primary line colour
    # lcol2       = '0.4'
    # lcol3       = '0.0'
    # lcols       = ['None', 'None', 'None', 'None', 'None', '0.0']
    # hatches     = [None, None, None, None, None, '//']

    # # Legend
    # llxbbox     =  0.0        # x-anchor legend bounding box
    # llybbox     = 1.15       # y-anchor legend bounding box
    # llrspace    = 0.          # spacing between rows in legend
    # llcspace    = 1.0         # spacing between columns in legend
    # llhtextpad  = 0.4         # the pad between the legend handle and text
    # llhlength   = 1.5         # the length of the legend handles
    # frameon     = False       # if True, draw a frame around the legend. If None, use rc

    # # PNG
    # dpi           = 600
    # transparent   = False
    # bbox_inches   = 'tight'
    # pad_inches    = 0.035

    # # Clock options
    # ymax = 0.8
    # ntextsize   = 'medium'       # normal textsize
    # # modules
    # bmod        = 0.5            # fraction of ymax from center to start module colours
    # alphamod    = 0.7            # alpha channel for modules
    # fwm         = 0.05           # module width to remove at sides
    # ylabel1     = 1.15           # position of module names
    # ylabel2     = 1.35           # position of class names
    # mtextsize   = 'large'        # 1.3*textsize # textsize of module labels
    # # bars
    # bpar        = 0.4            # fraction of ymax from center to start with parameter bars
    # fwb         = [0.7,0.4,0.3]  # width of bars
    # plwidth     = 0.5
    # # parameters in centre
    # bplabel     = 0.1            # fractional distance of ymax of param numbers in centre from 0-line
    # ptextsize   = 'medium'       # 'small' # 0.8*textsize # textsize of param numbers in centre
    # # yaxis
    # space4yaxis = 2              # space for y-axis (integer)
    # ytextsize   = 'medium'       # 'small' # 0.8*textsize # textsize of y-axis

    # import matplotlib as mpl
    # if (outtype == 'pdf'):
    #     mpl.use('PDF') # set directly after import matplotlib
    #     import matplotlib.pyplot as plt
    #     from matplotlib.backends.backend_pdf import PdfPages
    #     # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    #     mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
    #     mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    #     if usetex:
    #         mpl.rc('text', usetex=True)
    #     else:
    #         #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #         mpl.rc('font',**{'family':'serif','serif':['times']})
    #     mpl.rc('text.latex', unicode=True)
    # elif (outtype == 'png'):
    #     mpl.use('Agg') # set directly after import matplotlib
    #     import matplotlib.pyplot as plt
    #     mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    #     if usetex:
    #         mpl.rc('text', usetex=True)
    #     else:
    #         #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #         mpl.rc('font',**{'family':'serif','serif':['times']})
    #     mpl.rc('text.latex', unicode=True)
    #     mpl.rc('savefig', dpi=dpi, format='png')
    # else:
    #     import matplotlib.pyplot as plt
    #     mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
    # mpl.rc('font', size=textsize)
    # mpl.rc('lines', linewidth=lwidth, color='black')
    # mpl.rc('axes', linewidth=alwidth, labelcolor='black')
    # mpl.rc('path', simplify=False) # do not remove

    # if (outtype == 'pdf'):
    #     print('Plot PDF ', pdffile)
    #     pdf_pages = PdfPages(pdffile)
    # elif (outtype == 'png'):
    #     print('Plot PNG ', pngbase)
    # else:
    #     print('Plot X')

    # ifig = 0
    # ifig += 1
    # iplot = 0
    # print('Plot - Fig ', ifig)
    # fig = plt.figure(ifig)

    # iplot += 1
    # sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), polar=True)
    # clockplot(sub, si, sti, usetex=usetex, dolegend=True)

    # iplot += 1
    # sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), polar=True)
    # clockplot(sub, si, sti, stierr, llybbox=1.00, usetex=usetex)

    # iplot += 1
    # sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace), polar=True)
    # clockplot(sub, [si,si[::-1]], [sti,sti[::-1]], [stierr,stierr[::-1]], llybbox=1.1, doabc=True, dolegend=True)

    # if (outtype == 'pdf'):
    #     pdf_pages.savefig(fig)
    #     plt.close(fig)
    # elif (outtype == 'png'):
    #     pngfile = pngbase+"{0:04d}".format(ifig)+".png"
    #     fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    #     plt.close(fig)

    # if (outtype == 'pdf'):
    #     pdf_pages.close()
    # elif (outtype == 'png'):
    #     pass
    # else:
    #     plt.show()
