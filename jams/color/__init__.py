#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Functions working with colours and producing colour tables.


    Provided colours
    ----------------
    Can be accessed also with colours(list of colours)
    ufzdarkblue, ufzblue, ufzlightblue,
    ufzred, ufzorange, ufzyellow,
    ufzdarkgreen, ufzgreen, ufzlightgreen,
    ufzgray1, ufzgray2, ufzgray3,
    ufzgrey1, ufzgrey2, ufzgrey3,
    ufzdarkgray, ufzgray, ufzlightgray,
    ufzdarkgrey, ufzgrey, ufzlightgrey,
    ufzblack, ufzwhite,
    darkblue=ufzdarkblue, blue=ufzblue, lightblue=ufzlightblue,
    red=ufzred, orange=ufzorange, yellow=ufzyellow,
    darkgreen=ufzdarkgreen, green=ufzgreen, lightgreen=ufzlightgreen,
    gray1=ufzgray1, gray2=ufzgray2, gray3=ufzgray3,
    grey1=ufzgrey1, grey2=ufzgrey2, grey3=ufzgrey3,
    darkgray=ufzdarkgray, gray=ufzgray, lightgray=ufzlightgray,
    darkgrey=ufzdarkgrey, grey=ufzgrey, lightgrey=ufzlightgrey,
    black=ufzblack, white=ufzwhite


    Provided colour dictionnaries
    -----------------------------
    Can be accessed also with get_brewer(colour palette name)
    brewer_sequential, brewer_diverging, brewer_qualitative, 
    oregon_sequential, oregon_diverging, oregon_qualitative, 
    ncl_large, ncl_small, ncl_meteo_swiss,
    mathematica
    chroma_brewer, chroma_x11


    Provided functions
    ------------------
    colors             Wrapper for colours.
    colours            Get JAMS colours.
    get_brewer         Register and return Brewer colormap.
    plot_brewer        Plot available Brewer color maps in pdf file.
    print_brewer       Print available Brewer colormap names.
    register_brewer    Register Brewer colormap with matplotlib.
    rgb_blend          Calculates colour between two given colors in rgb space.
    rgb_gradient       n interpolated colours in rgb space between several colours.
    rgb_range          n interpolated colours between two colours in rgb space.
    limit              Limit number between given min and max
    luminance          Relative brightness normalized to 0 for darkest black and 1 for lightest white
    rgb2hex,   hex2rgb,   hex2rgb01    Converter between different colour spaces:
    rgb2hsi,   hsi2rgb,   hsi2rgb01    hex string, HSI, HSV, HSL, Lab, LcH
    rgb2hsl,   hsl2rgb,   hsl2rgb01    RGB [0-255], RGB01 [0-1]
    rgb2hsv,   hsv2rgb,   hsv2rgb01    "
    rgb2lab,   lab2rgb,   lab2rgb01    "
    lab2lch,   lch2lab                 "
    rgb2lch,   lch2rgb,   lch2rgb01    "
    rgb2rgb01, rgb012rgb               "
    col2rgb,   col2rgb01               Convert any colour known by matplotlib (plus RGB [0-255]) to RGB [0-1] or [0-255]
    bezier                             Bezier interpolation between 2-5 colours in Lab colour space
    sron_colors        Distinct colour palettes of Paul Tol at SRON - Netherlands Institute for Space Research
    sron_maps          Colour maps of Paul Tol at SRON - Netherlands Institute for Space Research


    Example
    -------
    >>> import numpy as np
    >>> from autostring import astr
    >>> print(astr(np.array(jams.color.ufzdarkblue), 4))
    ['0.0000' '0.2431' '0.4314']

    >>> print(colours('ufzdarkblue', rgb256=True))
    (0, 62, 110)

    >>> print(colours(names=True)[0:3])
    ['ufzdarkblue', 'ufzblue', 'ufzlightblue']

    >>> print(astr(np.array(jams.color.colours('JAMSDARKBLUE')), 4))
    ['0.0000' '0.2431' '0.4314']

    >>> print(astr(np.array(jams.color.colours('DarkBlue')), 4))
    ['0.0000' '0.2431' '0.4314']

    >>> print(jams.color.colours(['orange','ufzdarkblue'], rgb256=True))
    [(207, 104, 0), (0, 62, 110)]


    # Print
    # jams.color.print_brewer('diverging')

    # Plot
    # jams.color.plot_brewer('brew_')

    # Get Names
    # names = jams.color.get_brewer(names='sequential')

    # Register colour map and get colour map handle
    # cc = jams.color.get_brewer('RdYlBu11')
    # plt.pcolormesh(np.outer(np.arange(cc.N), np.ones(cc.N)), cmap=cc)

    # Get RGB colours of colour map
    # cc = jams.color.get_brewer('blues4', rgb=True)
    # mark1 = sub.plot(x, y)
    # plt.setp(mark1, linestyle='None', marker='o', markeredgecolor=cc[0], markerfacecolor='None')

    # Register and use colour map
    # from scipy import misc
    # lena = misc.lena()
    # plt.imshow(l, cmap=mpl.cm.rainbow)
    # jams.color.register_brewer('ncl_meteo_swiss')
    # plt.imshow(l, cmap=mpl.cm.get_cmap('hotcold_18lev'))
    # plt.imshow(l, cmap=jams.get_brewer('hotcold_18lev'))


    >>> print(jams.color.brewer_sequential['blues4'])
    [(239, 243, 255), (189, 215, 231), (107, 174, 214), (33, 113, 181)]

    >>> jams.color.print_brewer('qualitative')[0:7]
    ['set33', 'set34', 'set35', 'set36', 'set37', 'set38', 'set39']

    >>> print(astr(np.array(jams.color.get_brewer('blues4', rgb=True)[0]), 4))
    ['0.9373' '0.9529' '1.0000']

    >>> print(jams.color.get_brewer('Blues4', rgb256=True)[0])
    (239, 243, 255)

    >>> cc = jams.color.get_brewer('bLuEs4', rgb256=True, reverse=True)
    >>> print(cc[-1])
    (239, 243, 255)
    >>> print(cc[0])
    (33, 113, 181)

    >>> print(astr(np.array(jams.color.get_brewer('blues4', rgb256=True, grey=True)[0]), 4))
    ['242.9897' '242.9897' '242.9897']


    >>> r = (1.0,0.0,0.0)
    >>> b = (0.0,0.0,1.0)
    >>> print(jams.color.rgb_blend(r,b,0.0), jams.color.rgb_blend(r,b,0.5), jams.color.rgb_blend(r,b,1.0))
    (1.0, 0.0, 0.0) (0.5, 0.0, 0.5) (0.0, 0.0, 1.0)

    >>> print(jams.color.rgb_range(r,b,3))
    [(1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0)]

    >>> print(jams.color.rgb_range(r,b,3,pow=2))
    [(1.0, 0.0, 0.0), (0.75, 0.0, 0.25), (0.0, 0.0, 1.0)]

    >>> print(jams.color.rgb_gradient([r,b],[0.0,1.0],3))
    [(1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0)]

    >>> print(jams.color.rgb_gradient([r,r,b,b],[0.0,0.25,0.75,1.0],5, cmap='MyGradient'))
    [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 0.0, 0.5), (0.0, 0.0, 1.0), (0.0, 0.0, 1.0)]


    >>> print(jams.color.chroma_brewer['OrRd'])
    ['#fff7ec', '#fee8c8', '#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f', '#b30000', '#7f0000']

    >>> print(jams.color.chroma_x11['indigo'])
    #4b0082

    >>> print(jams.color.rgb2hex(1, 101, 201))
    #0165c9


    >>> print(jams.color.limit(-1))
    0

    >>> print(jams.color.limit(2))
    1

    >>> print(jams.color.limit(267, mini=0, maxi=255))
    255

    >>> print(jams.color.luminance(*(0,0,0)))
    0.0

    >>> print(jams.color.luminance(*(255,255,255)))
    1.0

    >>> print(jams.color.luminance(*(255,0,0)))
    0.2126


    >>> print(jams.color.hex2rgb(jams.color.rgb2hex(1, 101, 201)))
    (1, 101, 201)

    >>> print(jams.color.hsi2rgb(*jams.color.rgb2hsi(1, 101, 201)))
    (1, 101, 201)

    >>> print(jams.color.hsl2rgb(*jams.color.rgb2hsl(1, 101, 201)))
    (1, 101, 201)

    >>> print(jams.color.hsv2rgb(*jams.color.rgb2hsv(1, 101, 201)))
    (1, 101, 201)

    >>> print(jams.color.lab2rgb(*jams.color.rgb2lab(1, 101, 201)))
    (1, 101, 201)

    >>> print(jams.color.lch2rgb(*jams.color.rgb2lch(1, 101, 201)))
    (1, 101, 201)


    import numpy as np
    from jams.color import chroma_brewer, hex2rgb01

    import matplotlib as mpl
    from matplotlib.pylab import *

    fig = figure(figsize=(8,6))
    ax1 = fig.add_axes([0.05, 0.90, 0.9, 0.10])
    ax2 = fig.add_axes([0.05, 0.75, 0.9, 0.10])
    ax3 = fig.add_axes([0.05, 0.60, 0.9, 0.10])
    ax4 = fig.add_axes([0.05, 0.45, 0.9, 0.10])
    ax5 = fig.add_axes([0.05, 0.30, 0.9, 0.10])
    ax6 = fig.add_axes([0.05, 0.15, 0.9, 0.10])
    ax7 = fig.add_axes([0.05, 0.00, 0.9, 0.10])

    # 7 sequential colours of increasing luminance
    nn     = 7
    cc     = jams.color.bezier(['black', 'red', 'yellow', 'white'], nn)
    cmap   = mpl.colors.ListedColormap(cc)
    norm   = mpl.colors.BoundaryNorm(np.arange(cmap.N+1), cmap.N)
    cb     = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')

    # 255 sequential colours of decreasing luminance
    cc     = jams.color.bezier(['white', 'yellow', 'red', 'black'], reverse=True)
    cmap   = mpl.colors.ListedColormap(cc)
    norm   = mpl.colors.BoundaryNorm(np.arange(cmap.N+1), cmap.N)
    cb     = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='horizontal')

    # 9 diverging colours from 5 given colours, first increasing then decreasing luminance
    nn     = 9
    cc     = jams.color.bezier(['darkred', 'deeppink', 'lightyellow', 'lightgreen', 'teal'], nn)
    cmap   = mpl.colors.ListedColormap(cc)
    norm   = mpl.colors.BoundaryNorm(np.arange(cmap.N+1), cmap.N)
    cb     = mpl.colorbar.ColorbarBase(ax3, cmap=cmap, norm=norm, orientation='horizontal')

    # 9 sequential colours of decreasing luminance
    nn     = 9
    cc     = jams.color.bezier([ jams.color.hex2rgb01(i) for i in jams.color.chroma_brewer['Oranges'][::3] ], nn)
    cmap   = mpl.colors.ListedColormap(cc)
    norm   = mpl.colors.BoundaryNorm(np.arange(cmap.N+1), cmap.N)
    cb     = mpl.colorbar.ColorbarBase(ax4, cmap=cmap, norm=norm, orientation='horizontal')

    # 5 sequential colours of decreasing luminance registered as MyBrewer colour map with matplotlib
    nn     = 5
    cc     = jams.color.bezier([ jams.color.rgb2rgb01(*i) for i in [(255,255,178), (253,141,60), (189,0,38)] ], nn, cmap='MyBrewer')
    norm   = mpl.colors.BoundaryNorm(np.arange(nn+1), nn)
    cb     = mpl.colorbar.ColorbarBase(ax5, cmap=mpl.cm.get_cmap('MyBrewer'), norm=norm, orientation='horizontal')

    # 9 diverging colours from 3 given colours, first increasing then decreasing luminance
    nn     = 9
    cc     = jams.color.bezier(['darkred', 'lightyellow', 'teal'], nn)
    cmap   = mpl.colors.ListedColormap(cc)
    norm   = mpl.colors.BoundaryNorm(np.arange(cmap.N+1), cmap.N)
    cb     = mpl.colorbar.ColorbarBase(ax6, cmap=cmap, norm=norm, orientation='horizontal')

    # 9 diverging colours from 3 given colours, first increasing then decreasing luminance
    # but interpolated in L*C*h space instead of L*a*b
    nn     = 9
    cc     = jams.color.bezier(['darkred', 'lightyellow', 'teal'], nn, lch=True)
    cmap   = mpl.colors.ListedColormap(cc)
    norm   = mpl.colors.BoundaryNorm(np.arange(cmap.N+1), cmap.N)
    cb     = mpl.colorbar.ColorbarBase(ax7, cmap=cmap, norm=norm, orientation='horizontal')

    show()


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2015 Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  MC, Mar 2015
"""

# UFZ colours
from .ufz_colours       import ufzdarkblue, ufzblue, ufzlightblue
from .ufz_colours       import ufzred, ufzorange, ufzyellow
from .ufz_colours       import ufzdarkgreen, ufzgreen, ufzlightgreen
from .ufz_colours       import ufzgray1, ufzgray2, ufzgray3
from .ufz_colours       import ufzgrey1, ufzgrey2, ufzgrey3
from .ufz_colours       import ufzdarkgray, ufzgray, ufzlightgray
from .ufz_colours       import ufzdarkgrey, ufzgrey, ufzlightgrey
from .ufz_colours       import ufzblack, ufzwhite
from .ufz_colours       import darkblue, blue, lightblue
from .ufz_colours       import red, orange, yellow
from .ufz_colours       import darkgreen, green, lightgreen
from .ufz_colours       import gray1, gray2, gray3
from .ufz_colours       import grey1, grey2, grey3
from .ufz_colours       import darkgray, gray, lightgray
from .ufz_colours       import darkgrey, grey, lightgrey
from .ufz_colours       import black, white
from .colours           import colours, colors

# colorbrewer
from .brewer_colors     import brewer_sequential, brewer_diverging, brewer_qualitative
from .brewer_colors     import oregon_sequential, oregon_diverging, oregon_qualitative
from .brewer_colors     import ncl_large, ncl_small, ncl_meteo_swiss
from .brewer_colors     import mathematica
from .brewer            import get_brewer, register_brewer, print_brewer, plot_brewer

# rgb interpolation
from .rgb               import rgb_blend, rgb_gradient, rgb_range

# chroma.js
from .chroma_colors     import chroma_brewer, chroma_x11
from .chroma_helper     import limit, luminance
from .chroma_converter  import rgb2hex,   hex2rgb,   hex2rgb01
from .chroma_converter  import rgb2hsi,   hsi2rgb,   hsi2rgb01
from .chroma_converter  import rgb2hsl,   hsl2rgb,   hsl2rgb01
from .chroma_converter  import rgb2hsv,   hsv2rgb,   hsv2rgb01
from .chroma_converter  import rgb2lab,   lab2rgb,   lab2rgb01
from .chroma_converter  import lab2lch,   lch2lab
from .chroma_converter  import rgb2lch,   lch2rgb,   lch2rgb01
from .chroma_converter  import rgb2rgb01, rgb012rgb
from .chroma_converter  import col2rgb,   col2rgb01
from .chroma_bezier     import bezier

# sron
from .sron              import sron_colors, sron_colours, sron_maps

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.5'
__revision__ = "Revision: 2419"
__date__     = 'Date: 16.05.2016'
