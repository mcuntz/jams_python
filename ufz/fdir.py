#!/usr/bin/env python
#
# This script reads a DEM from a netcdf file and calculates the flow direction
# following ArcGis convention:
#
#           64               N
#       32      128          A
#   16      -1       1      /|\
#        8       2           |
#            4               |
#
# CAUTION*** origin is in the lower left NOT upper right ***CAUTION
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
# author:  Stephan Thober
# created: 01.07.2015
import numpy as np

# global variables for calculating upstream cells
yy_offset, xx_offset = np.meshgrid(np.arange(-1, 2), np.arange(-1, 2))
yy_offset = yy_offset.ravel() # yy_offset for neighboring cells
xx_offset = xx_offset.ravel() # xx_offset for neighboring cells
# local flow direction sorted in the order of yy_offset and xx_offset
local_flow_direction = np.array([8, 16, 32, 4, -1, 64, 2, 1, 128])
# same as local_flow_direction but reverted. This means these are the
# flow directions of neighboring cells flowing into the local cell.
inflow_direction = np.array([128, 1, 2, 64, -9999, 4, 32, 16, 8])


def get_neighbors(arr, yy_loc, xx_loc):
    # global variables used: yy_offset, xx_offset
    yy_ind = yy_offset + yy_loc
    xx_ind = xx_offset + xx_loc
    # create mask for valid neighbors
    neighbors = ((yy_ind >= 0) &
                 (yy_ind < arr.shape[0]) &
                 (xx_ind >= 0) &
                 (xx_ind < arr.shape[1]))
    return neighbors, yy_ind[neighbors], xx_ind[neighbors]


def flow_direction(fd, dem):
    # global variable used: correct_direction
    #
    for ii in np.arange(dem.shape[0]):
        for jj in np.arange(dem.shape[1]):
            if dem.mask[ii, jj]:
                continue
            # get mask of neighbors and y and x locations
            neighbors, yy, xx = get_neighbors(dem, ii, jj)
            # get position of cell with steepest gradient
            pos_min = np.ma.argmin(dem[yy, xx] - dem[ii, jj])
            fd[ii, jj] = local_flow_direction[neighbors][pos_min]
    return fd


def get_sinks(fd, sink=-1):
    return np.ma.where(fd == sink)


def get_upstream(fd, loc):
    # global variable used: inflow_direction
    #
    # get mask of neighbors and y and x locations
    neighbors, yy, xx = get_neighbors(fd, loc[0], loc[1])
    # mask inflowing cells
    upstream_mask = (fd.data[yy, xx] == inflow_direction[neighbors])
    yy_upstream = yy[upstream_mask]
    xx_upstream = xx[upstream_mask]
    return [[yy_upstream[ii], xx_upstream[ii]] for ii in np.arange(np.sum(upstream_mask))]


def network_properties(fd, yy, xx, print_info=False, do_co=True, co=None, do_fa=True, fa=None):
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
    # flow direction stack to emulate recursion
    if not do_co and not do_fa:
        raise ValueERROR('***ERROR: neither fa nor co calculated')
    fd_stack = [[yy, xx]] # start at initial sink
    while fd_stack:        
        if print_info:
            print('current flow accumulation stack: ', fd_stack)
        upstream = get_upstream(fd, fd_stack[-1])
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

if __name__ == '__main__':
    print('Done!')
