#!/usr/bin/env python
from __future__ import print_function
#
# This script reads a DEM from a netcdf file and calculates the flow direction
# following ArcGis convention:
#
# flow direction is assumed like this
#           64             y-axis
#       32      128          |
#   16      -1       1       |
#        8       2          \|/
#            4               V
#      x-axis ------>
# ORIGIN is in the upper left corner
#
# sinks are marked by -1
#
# Nomenclature:
#
# co - channel orders
# fd - flow direction
# fa - flow accumulation
# yy - row index
# xx - column index
#
# author:  Stephan Thober, David Schaefer
# created: 01.07.2015
import numpy as np

# global variables for calculating upstream cells
# yy_offset, xx_offset = np.meshgrid(np.arange(-1, 2), np.arange(-1, 2))
yy_offset, xx_offset = np.meshgrid(np.arange(1, -2, -1), np.arange(-1, 2))
yy_offset = yy_offset.ravel() # yy_offset for neighboring cells
xx_offset = xx_offset.ravel() # xx_offset for neighboring cells
# local flow direction sorted in the order of yy_offset and xx_offset
local_flow_direction = np.array([8, 16, 32, 4, -1, 64, 2, 1, 128])
# same as local_flow_direction but reverted. This means these are the
# flow directions of neighboring cells flowing into the local cell.
inflow_direction = np.array([128, 1, 2, 64, -9999, 4, 32, 16, 8])


class stream_network(object):
    def __init__(self, dem=None, fdir=None, co=None, do_co=False, fa=None, do_fa=False):
        # initialize all arrays
        self.dem = None # digital elevation model
        self.fdir = None # flow direction
        self.sinks = None # sinks
        self.co = None # channel order
        self.fa = None # flow accumulation
        # consistency check
        if fdir == None and dem == None:
            raise ValueError('***ERROR: specify either dem or fdir to create a stream_network object')
        if fdir == None:
            # create flow direction if necessary
            self.dem = dem
            self.fdir = self.flow_direction()
        else:
            self.fdir = fdir
        # assign flow accumulation
        if not fa == None:
            self.fa = fa
        # assign channel order
        if not co == None:
            self.co = co
        # assign sinks
        self.sinks = self._get_sinks()
        # get channel order and flow accumulation
        if do_co and do_fa:
            self.co = np.zeros(self.fdir.shape)
            self.fa = np.zeros(self.fdir.shape)
            for ii in np.arange(self.sinks[0].shape[0]):
                co, fa = self.network_properties(self.fdir, self.sinks[0][ii], self.sinks[1][ii],
                                            do_co=do_co, co=self.co,
                                            do_fa=do_fa, fa=self.fa)
        elif do_co and not do_fa:
            self.co = np.zeros(self.fdir.shape)
            for ii in np.arange(self.sinks[0].shape[0]):
                self.co = self.network_properties(self.fdir, self.sinks[0][ii], self.sinks[1][ii],
                                                  do_co=do_co, co=self.co,
                                                  do_fa=do_fa)
        elif not do_co and do_fa:
            self.fa = np.zeros(self.fdir.shape)
            for ii in np.arange(self.sinks[0].shape[0]):
                self.fa = self.network_properties(self.fdir, self.sinks[0][ii], self.sinks[1][ii],
                                                  do_co=do_co,
                                                  do_fa=do_fa, fa=self.fa)
    
    def flow_direction(self):
        # global variable used: correct_direction
        fd = np.zeros(self.dem.shape)
        #
        for ii in np.arange(self.dem.shape[0]):
            for jj in np.arange(self.dem.shape[1]):
                if self.dem.mask[ii, jj]:
                    continue
                # get mask of neighbors and y and x locations
                neighbors, yy, xx = _get_neighbors(self.dem, ii, jj)
                # get position of cell with steepest gradient
                pos_min = np.ma.argmin(self.dem[yy, xx] - self.dem[ii, jj])
                fd[ii, jj] = local_flow_direction[neighbors][pos_min]
        return fd


    def network_properties(self, fd, yy, xx, print_info=False, do_co=True, co=None, do_fa=True, fa=None):
        """
            Calculates channel order number and flow accumulation starting from one sink in a flow direction map


            Definition
            ----------
            def network_properties(self, fd, yy, xx, print_info=False, do_co=True, co=None, do_fa=True, fa=None):


            Input
            -----
            file         self - stream_network object
            fd           flow direction field, basically stream_network.fd
            yy           row coordinate of sink
            xx           column coordinate of sink


            Optional Input Parameters
            -------------------------
            print_info   write additional info on std_out
            do_co        calculate channel order
            co           given channel order field
            do_fa        calculate flow accumulation
            fa           given flow accumulation

            Options
            -------
            print_info   True: write additional information
                         False: do not write additional information (default)
            do_co        True: calculate channel order (default)
                         False: do not channel order
            co           None: no channel order field specified, will be created (default)
            do_fa        True: calculate flow accumulation (default)
                         False: do not flow accumulation
            fa           None: no flow accumulation field specified, will be created (default)

            Output
            ------
            Depending on options:
                co, fa if do_co=True and do_fa=True
                co if do_co=True and not do_fa=True
                fa if not do_co=True and do_fa=True


            Restrictions
            ------------
                The origin of the flow direction field is assumed to be located in the upper left corner.
                Then flow directions are following this convention:

                flow direction is assumed like this
                          64             y-axis
                      32      128          |
                  16      -1       1       |
                       8       2          \|/
                           4               V
                     x-axis ------>

            Examples
            --------
            >>> # Create some data
            >>> fd = np.ma.array([[  2,   1,   1,  2,   4,   4,  16,  8,  8],
            ...          [  1,   2,   1,  1,   4,   4,   4,  4, 16],
            ...          [128,   1,   1,  1,   1,   2,   4,  4, 16],
            ...          [  1,   1,  64, 64, 128,   1,   1,  4,  8],
            ...          [  1,  64,  64, 64, 128,  64,   4,  4,  4],
            ...          [ 64, 128,  64, 32, 128,   1,   4,  2,  4],
            ...          [  1,  64,  64, 64,   1,   1,   1,  1, -1],
            ...          [  1,   2, 128, 64, 128,   1, 128, 64, 64],
            ...          [128,   1, 128, 64, 128, 128,  64, 64, 16]])
            >>> stream_network(fdir=fd, do_fa=False, do_co=True).co
            array([[ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
                   [ 1.,  2.,  1.,  1.,  2.,  1.,  1.,  2.,  1.],
                   [ 1.,  1.,  3.,  3.,  3.,  3.,  1.,  2.,  1.],
                   [ 1.,  2.,  3.,  1.,  1.,  2.,  3.,  3.,  1.],
                   [ 1.,  1.,  2.,  1.,  1.,  1.,  1.,  3.,  1.],
                   [ 1.,  1.,  1.,  2.,  1.,  1.,  2.,  3.,  1.],
                   [ 1.,  1.,  1.,  2.,  1.,  2.,  3.,  3.,  4.],
                   [ 1.,  2.,  1.,  2.,  1.,  1.,  2.,  1.,  1.],
                   [ 1.,  1.,  2.,  1.,  1.,  1.,  1.,  1.,  1.]])
            >>> stream_network(fdir=fd, do_fa=True, do_co=True).fa
            array([[  1.,   1.,   2.,   3.,   1.,   2.,   1.,   1.,   1.],
                   [  1.,   4.,   1.,   2.,   7.,   3.,   2.,   3.,   1.],
                   [  1.,   1.,  28.,  31.,  39.,  44.,   3.,   5.,   1.],
                   [  1.,   5.,  22.,   2.,   1.,   4.,  52.,  58.,   1.],
                   [  2.,   3.,  16.,   1.,   1.,   2.,   1.,  60.,   1.],
                   [  1.,   3.,   2.,  10.,   1.,   1.,   3.,  61.,   2.],
                   [  1.,   2.,   1.,   9.,   1.,   3.,   7.,  16.,  81.],
                   [  1.,   3.,   1.,   7.,   1.,   2.,   5.,   3.,   1.],
                   [  1.,   1.,   5.,   1.,   1.,   1.,   1.,   2.,   1.]])

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

            Copyright 2009-2015 Stephan Thober, David Schaefer


            History
            -------
            Written,  ST, Dec 2015
            Modified                
        """
        # calculate channel order according to Strahler 1952 and flow accumulation
        #
        # do_co - flag for calculating channel order
        # do_fa - flag for calculating flow accumulation
        #
        # channel order of headwater is ONE, if channels join, the channel order
        # of the resulting stream is the highest one of the inflowing streams,
        # if two or more than two inflowing streams have the highest channel order,
        # the channel order of the resulting stream is one higher than the highest
        # channel order of the inflowing streams
        # -------------------------------------------------------------------------
        if co == None:
            co = np.zeros(fd.shape)
        if fa == None and do_fa:
            fa = np.zeros(fd.shape)
        # flow direction stack to emulate recursion
        if not do_co and not do_fa:
            raise ValueERROR('***ERROR: neither fa nor co calculated')
        fd_stack = [[yy, xx]] # start at initial sink
        while fd_stack:        
            if print_info:
                print('current flow accumulation stack: ', fd_stack)
            upstream = self._get_upstream(fd, fd_stack[-1])
            if print_info:
                print('upstream locations: ', upstream)
            ext = [l for l in upstream if co[l[0],l[1]] == 0]
            if ext:
                fd_stack.extend(ext)
                continue
            # all upstream cells are available
            # note that headwaters dont have an upstream cell
            cell = fd_stack.pop() # save cell
            if do_co:
                co_upstream = [co[loc[0], loc[1]] for loc in upstream]
                co_max = np.amax(co_upstream + [1])      
                if len(np.where(co_upstream == co_max)[0]) > 1:
                    co_max += 1
                co[cell[0], cell[1]] = co_max
                if print_info:
                    print('co (channel order) of upstream: ', co_upstream)
            if do_fa:
                fa_upstream = [fa[loc[0], loc[1]] for loc in upstream]        
                fa[cell[0], cell[1]] = np.sum(fa_upstream) + 1
                if print_info:
                    print('sum of upstream: ', fa_upstream)
        if do_co and do_fa:
            return co, fa
        elif do_co and not do_fa:
            return co
        elif not do_co and do_fa:
            return fa

    def _get_sinks(self, sink=-1):
        return np.ma.where(self.fdir == sink)


    def _get_neighbors(self, arr, yy_loc, xx_loc):
        # global variables used: yy_offset, xx_offset
        #
        yy_ind = yy_offset + yy_loc
        xx_ind = xx_offset + xx_loc
        # create mask for valid neighbors
        neighbors = ((yy_ind >= 0) &
                     (yy_ind < arr.shape[0]) &
                     (xx_ind >= 0) &
                     (xx_ind < arr.shape[1]))
        return neighbors, yy_ind[neighbors], xx_ind[neighbors]


    def _get_upstream(self, fd, loc):
        # global variable used: inflow_direction
        #
        # get mask of neighbors and y and x locations
        neighbors, yy, xx = self._get_neighbors(fd, loc[0], loc[1])
        # mask inflowing cells
        upstream_mask = (fd.data[yy, xx] == inflow_direction[neighbors])
        yy_upstream = yy[upstream_mask]
        xx_upstream = xx[upstream_mask]
        return [[yy_upstream[ii], xx_upstream[ii]] for ii in np.arange(np.sum(upstream_mask))]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    print('Done')
