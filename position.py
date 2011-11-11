#!/usr/bin/env python
import numpy as np

def position(row=1, col=1, num=1,
             left=0.125, right=0.9, bottom=0.1, top=0.9,
             wspace=0.1, hspace=0.1,
             sortcol=False, golden=False, inversegolden=False,
             figsize=(1.,1.)):
    """
        Gives positions of subplots.
        To be used with add_axes instead of subplot.
        All dimensions are fractions of the figure width or height.
        Figure and subplot spaces are the same as for figure.subplotparams
        except for hspace and wspace (halved).
        If the figsize keyword is given, a rectangular section of the figure
        will be used.

        Definition
        ----------
        def position(row=1, col=1, num=1, 
                     left=0.125, right=0.9, bottom=0.1, top=0.9,
                     wspace=0.1, hspace=0.1,
                     sortcol=False, golden=False, inversegolden=False,
                     figsize(1.,1.)):

        Optional Input
        -----
        row            number of subplot rows (default 1)
        col            number of subplot columns (default 1)
        num            subplot number (default 1)
        left           left border of plot (default 0.125)
        right          right border of plot (default 0.9)
        bottom         bottom border of plot (default 0.1)
        top            top border of plot (default 0.9)
        hspace         space between columns (default 0.1)
        wspace         space between rows (default 0.1)
        sortcol        fill columns then rows (default False)
        golden         golden ratio of width/height = (1+sqrt(5))/2
                       (default False)
        inversegolden  golden ratio of height/width
                       (overwritten by golden) (default False)
        figsize        (width, height) of figure as given by e.g.
                       matplotlib.rcParams['figure.figsize'].
                       Scales everything to rectangular section
                       (default (1,1))

        Output
        ------
        position array with [left, bottom, width, height)
        to be used with fig.add_axes.
    
        Example
        -------
        # Use, for example, as follows
        # fig1 = figure(1)
        # sub1 = fig1.add_axes(position(2,2,1))
        # sub2 = fig1.add_axes(position(2,2,2))

        # if you want to have a true rectangle
        # figsize = matplotlib.rcParams['figure.figsize']
        # sub = fig1.add_axes(position(1,1,1,figsize=figsize,left=0.1))

        # if you want to have a true golden ratio
        # sub = fig1.add_axes(position(1,1,1,figsize=figsize,golden=True))

        # Doctest examples
        >>> print position(2,2,1)
        [ 0.125   0.55    0.3375  0.35  ]
        >>> print position(2,2,1,sortcol=True)
        [ 0.125   0.55    0.3375  0.35  ]
        >>> print position(2,2,1,golden=True)
        [ 0.125       0.40858647  0.3375      0.20858647]
        >>> print position(2,2,1,inversegolden=True)
        [ 0.125      0.55       0.2163119  0.35     ]
        >>> print position(2,2,1,golden=True,sortcol=True)
        [ 0.125       0.40858647  0.3375      0.20858647]
        >>> print position(2,2,1,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.)
        [ 0.   0.5  0.5  0.5]
        >>> print position(2,2,2,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.)
        [ 0.5  0.5  0.5  0.5]
        >>> print position(2,2,3,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.)
        [ 0.   0.   0.5  0.5]
        >>> print position(2,2,4,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.)
        [ 0.5  0.   0.5  0.5]
        >>> print position(2,2,1,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.,golden=True)
        [ 0.          0.30901699  0.5         0.30901699]
        >>> print position(2,2,2,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.,golden=True)
        [ 0.5         0.30901699  0.5         0.30901699]
        >>> print position(2,2,3,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.,golden=True)
        [ 0.          0.          0.5         0.30901699]
        >>> print position(2,2,4,top=1.,bottom=0.,left=0.,right=1.,hspace=0.,wspace=0.,golden=True)
        [ 0.5         0.          0.5         0.30901699]
        >>> figsize=[8,11]
        >>> print position(2,2,1,golden=True,sortcol=True,figsize=figsize)
        [ 0.125       0.32442652  0.3375      0.15169925]
        >>> print position(2,2,1,figsize=figsize,left=0.1)
        [ 0.1         0.42727273  0.35        0.25454545]
        >>> print position(2,2,1,figsize=figsize,left=0.1,golden=True)
        [ 0.1         0.33004502  0.35        0.15731774]

        
        History
        -------
        Written, MC, Aug. 2009
    """
    #
    # Check
    nplots = row*col
    if num > nplots:
        print 'POSITION: num > number of plots %s.' % (num, nplots)
        return None
    if right-left <= 0.:
        print 'POSITION: right %s < left %s.' % (right, left)
        return None
    if top-bottom <= 0.:
        print 'POSITION: top %s < bottom %s.' % (top, bottom)
        return None
    #
    # Scaling to figsize
    scalex = figsize[1]/float(max(figsize))
    scaley = figsize[0]/float(max(figsize))
    #
    # width, height
    dx = (right-left-(col-1)*hspace)/col
    dy = (top-bottom-(row-1)*wspace)/row
    #
    # golden ratio
    ratio = (1.+np.sqrt(5.))/2.
    if golden:
        width = dx
        height = dx / ratio
        checkheight = (top-bottom-row*height) - (row-1)*wspace
        if checkheight < 0.:
            height = dy
            width = dy * ratio
            checkwidth = (right-left-col*width) - (col-1)*hspace
            if checkwidth < 0.:
                print 'POSITION: golden ratio does not work. Have to recode.'
                return None
    else:
        if inversegolden:
            height = dy
            width = dy / ratio
            checkwidth = (right-left-col*width) - (col-1)*hspace
            if checkwidth < 0.:
                width = dx
                height = dx * ratio
                checkheight = (top-bottom-row*height) - (row-1)*wspace
                if checkheight < 0.:
                    print 'POSITION: inverse golden ratio does not work. Have to recode.'
                    return None
        else:
            width = dx
            height = dy
    #
    # order row/colmn, column/row
    if sortcol:
        irow = (num-1) % row
        icol = (num-1) // row
    else:
        irow = (num-1) // col
        icol = (num-1) % col
    #
    # position
    pos = np.empty(4)
    pos[0] = left   + icol*(width+hspace)          *scalex
    pos[1] = bottom + (row-1-irow)*(height+wspace) *scaley
    pos[2] = width  *scalex
    pos[3] = height *scaley
    #
    return pos

if __name__ == '__main__':
    import doctest
    doctest.testmod()
