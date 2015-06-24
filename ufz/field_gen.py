#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

import numpy as np
import numpy.random as r
import sys


"""Note: The different classes are published under different licenses."""


class Field(object):
    """Generates a 2d random hydraulic conductivity field.

        Use Field.K(x, y) for calculating the field.

        A Kraichnan algorithm to calculate the field is used.
        At the moment only a Gaussian auto-correlation function is used,
        but this could change in the future.
        The field can either be generated
        -on a regular grid,
        -for a given tuple of positions (e.g. for particles)
        -or for single positions.
        This behaviour can be changed by setting the "mode" to
        "grid", "tuple", or "single". "grid" is the standard setting.

        Args:
            K_G (float): geometric mean of hydraulic conductivity (exp(mu))
            variance (float): variance of the conductivity field
            corr_len (float or array_like): correlation length of the
                conductivity field, if a single value is passed,
                isotropy is assumed
            no_random_modes (int): number of random modes to approximate the
                 Fourier sum
            seed (int, optional): set the seed of the master RNG, if "None",
                a random seed is used

        Examples:
        >>> f = Field(1, 1, [1, 1], 100, 15011997)
        >>> f.mode = 'single'
        >>> round(f.K(0,0), 4)
        0.2104
        >>> x_tuple = [ 4, 0, 3]
        >>> y_tuple = [-1, 0, 1]
        >>> f.mode = 'tuple'
        >>> round(f.K(x_tuple, y_tuple)[1], 4)
        0.2104
        >>> x_grid = np.arange(0, 5, 0.5)
        >>> y_grid = np.arange(0, 5, 1)
        >>> f.mode = 'grid'
        >>> round(f.K(x_grid, y_grid)[0,0], 4)
        0.2104

        #testing reshaping in all combinations
        >>> f.mode = 'grid'
        >>> f.mode = 'single'
        >>> f.mode = 'grid'
        >>> f.mode = 'tuple'
        >>> f.mode = 'grid'
        >>> f.mode = 'single'
        >>> f.mode = 'tuple'
        >>> f.mode = 'single'

        License:
        This class is part of the UFZ Python package.

        The UFZ Python package is free software: you can redistribute it
        and/or modify it under the terms of the GNU Lesser General Public
        License as published by the Free Software Foundation, either version 3
        of the License, or (at your option) any later version.

        The UFZ Python package is distributed in the hope that it will be
        useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
        Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public
        License along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2015 Lennart Schüler
    """
    def __init__(self, K_G, variance, corr_len, no_random_modes,
                 seed=None):
        #geometric mean
        self.K_G = K_G
        #standard deviation
        self.sigma = np.sqrt(variance)

        #correlation length
        #assuming isotropy, if scalar correlation length is given
        if not hasattr(corr_len, '__getitem__'):
            self.l = np.array((corr_len, corr_len))
        else:
            self.l = np.array(corr_len)

        self.N = no_random_modes

        self._mode = 'grid'
        self.K = self.K_grid

        self.reset(seed)

    def reset(self, seed=None):
        """Resets the field with a new seed.

        Args:
            seed (int, optional): master seed for different RNG streams
        """
        r1, r2, r3 = self.init_RNG(seed)

        self.k1 = r1.normal(0., 1./self.l[0], self.N)
        self.k2 = r2.normal(0., 1./self.l[1], self.N)
        self.phi = r3.uniform(0., 2*np.pi, self.N)
        #reshape according to mode
        self.mode = self._mode

    def init_RNG(self, seed=None):
        """Initialises different independant RNG streams from a master RNG.

        Args:
            seed (int, optional): master seed for different RNG streams
        """
        self.__master_RNG = r.RandomState(seed)
        # MC - seed with maximum 4294967295 instead of sys.maxsize for Mac OS X
        # self.master_RNG = \
        #     lambda: self.__master_RNG.random_integers(sys.maxsize-1)
        self.master_RNG = \
            lambda: self.__master_RNG.random_integers(4294967295)
        r1 = r.RandomState(self.master_RNG())
        r2 = r.RandomState(self.master_RNG())
        r3 = r.RandomState(self.master_RNG())
        return (r1, r2, r3)

    @property
    def mode(self):
        """Returns the mode (grid, tuple, single)."""
        return self._mode

    @mode.setter
    def mode(self, m):
        """Sets the mode (grid, tuple, single).

        Args:
            m (string): mode of calculation
        """
        if m == 'grid':
            self._mode = m
            self.K = self.K_grid
            self.k1 = np.squeeze(self.k1)
            self.k2 = np.squeeze(self.k2)
            self.phi = np.squeeze(self.phi)
            self.k1 = np.reshape(self.k1, (1, 1, len(self.k1)))
            self.k2 = np.reshape(self.k2, (1, 1, len(self.k2)))
            self.phi = np.reshape(self.phi, (1, 1, len(self.phi)))
        elif m == 'tuple':
            self._mode = m
            self.K = self.K_tuple
            self.k1 = np.squeeze(self.k1)
            self.k2 = np.squeeze(self.k2)
            self.phi = np.squeeze(self.phi)
            self.k1 = np.reshape(self.k1, (1, len(self.k1)))
            self.k2 = np.reshape(self.k2, (1, len(self.k2)))
            self.phi = np.reshape(self.phi, (1, len(self.phi)))
        elif m == 'single':
            self._mode = m
            self.K = self.K_single
            try:
                self.k1 = np.squeeze(self.k1, axis=1)
                self.k2 = np.squeeze(self.k2, axis=1)
                self.phi = np.squeeze(self.phi, axis=1)
            except:
                pass
        else:
            raise ValueError('Unknown mode {}'.format(m))

    def reshape_axis_to_grid(self, x, y):
        """Reshape given positions for vectorisation of grid calculation.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            (x, y) (array_like, array_like): reshaped positions
        """
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        x = np.reshape(x, (len(x), 1, 1))
        y = np.reshape(y, (1, len(y), 1))
        return (x, y)

    def reshape_axis_to_tuple(self, x, y):
        """Reshape given positions for vectorisation of tuple calculation.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            (x, y) (array_like, array_like): reshaped positions
        """
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        x = np.reshape(x, (len(x), 1))
        y = np.reshape(y, (len(y), 1))
        return (x, y)

    def K_grid(self, x, y):
        """Calculate the hydraulic conductivity K for given grid (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (2d array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_grid(x, y)
        return self.K_G * np.exp(self.f(x, y))

    def K_single(self, x, y):
        """Calculate the hydraulic conductivity K for given point (x, y).

        Args:
            x (float): x position
            y (float): y position
        Returns:
            K (float): log-hydraulic conductivity
        """
        return float(self.K_G * np.exp(self.f(x, y)))

    def K_tuple(self, x, y):
        """Calculate the hydraulic conductivity K for given tuples x, y.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (2d array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_tuple(x, y)
        return self.K_G * np.exp(self.f(x, y))

    def f(self, x, y):
        """Calculate the Fourier sum for the random field.

        Args:
            x (float or array_like): x position
            y (float or array_like): y position
        """
        return (self.sigma * np.sqrt(2./self.N) *
               np.sum(np.cos(self.k1*x + self.k2*y + self.phi), axis=-1))


class Incompr_Field(Field):
    """Generates a 2d random incompressible velocity field.

        Use Field.U(x, y) to obtain the velocity vector(s) (u, v) or
        Field.u(x, y) and Field.v(x, y) to obtain the components individually.
        The second method should be slightly slower (though not tested).

        A Kraichnan algorithm together with a projector to ensure
        incompressiblity is used to calculate the field.
        At the moment only a Gaussian auto-correlation function is used,
        but this could change in the future.
        The field can either be generated
        -on a regular grid,
        -for a given tuple of positions (e.g. for particles)
        -or for single positions.
        This behaviour can be changed by setting the "mode" to
        "grid", "tuple", or "single". "grid" is the standard setting.

        If you use the K-methods for calculating the conductivities via this
        class, you have to set the value for K_G after instantiating. Otherwise
        a value of K_G = 1 is assumed.

        Args:
            mean_velocity (float): mean velocity in e1-direction
            variance (float): variance of the conductivity field
            corr_len (float or array_like): correlation length of the
                conductivity field, if a single value is passed,
                isotropy is assumed
            no_random_modes (int): number of random modes to approximate the
                 Fourier sum
            seed (int, optional): set the seed, if "None",
                a random seed is used

        Examples:
        >>> f = Incompr_Field(1, 1, [1, 1], 100, 16052001)
        >>> f.mode = 'single'
        >>> np.around(f.U(0, 0), 4)
        array([-0.1997,  0.3987])
        >>> round(f.u(0, 0), 4)
        -0.1997
        >>> round(f.v(0, 0), 4)
        0.3987
        >>> round(f.K(0,0), 4)
        4.6852
        >>> x_tuple = [ 4, 0, 3]
        >>> y_tuple = [-1, 0, 1]
        >>> f.mode = 'tuple'
        >>> np.around(f.U(x_tuple, y_tuple)[0], 4)
        array([ 0.0112, -0.1997,  0.9485])
        >>> round(f.u(x_tuple, y_tuple)[1], 4)
        -0.1997
        >>> round(f.v(x_tuple, y_tuple)[1], 4)
        0.3987
        >>> round(f.K(x_tuple, y_tuple)[1], 4)
        4.6852
        >>> x_grid = np.arange(0, 5, 0.5)
        >>> y_grid = np.arange(0, 5, 1)
        >>> f.mode = 'grid'
        >>> round(f.U(x_grid, y_grid)[0][0,0], 4)
        -0.1997
        >>> round(f.U(x_grid, y_grid)[1][0,0], 4)
        0.3987
        >>> round(f.u(x_grid, y_grid)[0,0], 4)
        -0.1997
        >>> round(f.v(x_grid, y_grid)[0,0], 4)
        0.3987
        >>> round(f.K(x_grid, y_grid)[0,0], 4)
        4.6852

        #testing reshaping in all combinations
        >>> f.mode = 'grid'
        >>> f.mode = 'single'
        >>> f.mode = 'grid'
        >>> f.mode = 'tuple'
        >>> f.mode = 'grid'
        >>> f.mode = 'single'
        >>> f.mode = 'tuple'
        >>> f.mode = 'single'

        License:
        This class is part of the UFZ Python package.

        The UFZ Python package is free software: you can redistribute it
        and/or modify it under the terms of the GNU Lesser General Public
        License as published by the Free Software Foundation, either version 3
        of the License, or (at your option) any later version.

        The UFZ Python package is distributed in the hope that it will be
        useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
        Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public
        License along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2015 Lennart Schüler
    """
    def __init__(self, mean_velocity, variance, corr_len, no_random_modes,
                 seed=None):
        self.u_bar = mean_velocity
        self.U = self.U_grid
        self.u = self.u_grid
        self.v = self.v_grid
        super(Incompr_Field, self).__init__(1., variance, corr_len,
                                            no_random_modes, seed)

    #TODO can super be used here!?
    @Field.mode.setter
    def mode(self, m):
        """Sets the mode (grid, tuple, single).

            Needs to be overridden, because the self.U, self.u and self.v
            attributes need to be adapted to the given mode.

        Args:
            m (string): mode of calculation
        """
        if m == 'grid':
            self.U = self.U_grid
            self.u = self.u_grid
            self.v = self.v_grid
            self.K = self.K_grid
            self.k1 = np.squeeze(self.k1)
            self.k2 = np.squeeze(self.k2)
            self.phi = np.squeeze(self.phi)
            self.k1 = np.reshape(self.k1, (1, 1, len(self.k1)))
            self.k2 = np.reshape(self.k2, (1, 1, len(self.k2)))
            self.phi = np.reshape(self.phi, (1, 1, len(self.phi)))
        elif m == 'tuple':
            self.U = self.U_tuple
            self.u = self.u_tuple
            self.v = self.v_tuple
            self.K = self.K_tuple
            self.k1 = np.squeeze(self.k1)
            self.k2 = np.squeeze(self.k2)
            self.phi = np.squeeze(self.phi)
            self.k1 = np.reshape(self.k1, (1, len(self.k1)))
            self.k2 = np.reshape(self.k2, (1, len(self.k2)))
            self.phi = np.reshape(self.phi, (1, len(self.phi)))
        elif m == 'single':
            self.U = self.U_single
            self.u = self.u_single
            self.v = self.v_single
            self.K = self.K_single
            try:
                self.k1 = np.squeeze(self.k1, axis=1)
                self.k2 = np.squeeze(self.k2, axis=1)
                self.phi = np.squeeze(self.phi, axis=1)
            except:
                pass
        else:
            raise ValueError('Unknown mode {}'.format(m))
        self._mode = m

    def U_grid(self, x, y):
        """Calculate the velocity vector for grid (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            (u, v) (2d array_like, 2d array_like): velocity vector
                on grid (x, y)
        """
        x, y = self.reshape_axis_to_grid(x, y)
        u = (self.u_bar - self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_x(x, y))
        v = (-self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_y(x, y))
        return (u, v)

    def u_grid(self, x, y):
        """Calculate the longitudinal velocity component for grid (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            u (2d array_like): velocity in e1-direction on grid (x, y)
        """
        x, y = self.reshape_axis_to_grid(x, y)
        return (self.u_bar - self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_x(x, y))

    def v_grid(self, x, y):
        """Calculate the transversal velocity component for grid (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            v (2d array_like): velocity in e2-direction on grid (x, y)
        """
        x, y = self.reshape_axis_to_grid(x, y)
        return (-self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_y(x, y))

    def U_single(self, x, y):
        """Calculate the velocity vector for point (x, y).

        Args:
            x (float): x position
            y (float): y position
        Returns:
            (u, v) (float, float): velocity vector at point (x, y)
        """
        u = float(self.u_bar - self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_x(x, y))
        v = float(-self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_y(x, y))
        return (u, v)

    def u_single(self, x, y):
        """Calculate the longitudinal velocity component for point (x, y).

        Args:
            x (float): x position
            y (float): y position
        Returns:
            u (float): velocity in e1-direction at point (x, y)
        """
        return float(self.u_bar - self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_x(x, y))

    def v_single(self, x, y):
        """Calculate the transversal velocity component for point (x, y).

        Args:
            x (float): x position
            y (float): y position
        Returns:
            v (float): velocity in e2-direction at point (x, y)
        """
        return float(-self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_y(x, y))

    def U_tuple(self, x, y):
        """Calculate the velocity vectors for tuples (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            (u, v) (array_like, array_like): velocity vectors for tuples (x, y)
        """
        x, y = self.reshape_axis_to_tuple(x, y)
        u = (self.u_bar - self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_x(x, y))
        v = (-self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_y(x, y))
        return (u, v)

    def u_tuple(self, x, y):
        """Calculate the longitudinal velocity component for tuples (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            u (array_like): velocity in e1-direction for tuples (x, y)
        """
        x, y = self.reshape_axis_to_tuple(x, y)
        return (self.u_bar - self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_x(x, y))

    def v_tuple(self, x, y):
        """Calculate the transversal velocity component for tuples (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            v (array_like): velocity in e2-direction for tuples (x, y)
        """
        x, y = self.reshape_axis_to_tuple(x, y)
        return (-self.sigma * self.u_bar *
            np.sqrt(2.0/self.N) * self.f_incompr_y(x, y))

    def f_incompr_x(self, x, y):
        """Calculate the (long.) Fourier sum for the incompr. random field.

        Args:
            x (float or array_like): x position
            y (float or array_like): y position
        """
        return np.sum(self.p_x() * np.cos(self.k1*x + self.k2*y + self.phi),
                      axis=-1)

    def f_incompr_y(self, x, y):
        """Calculate the (trans.) Fourier sum for the incompr. random field.

        Args:
            x (float or array_like): x position
            y (float or array_like): y position
        """
        return np.sum(self.p_y() * np.cos(self.k1*x + self.k2*y + self.phi),
                      axis=-1)

    def p_x(self):
        """Projector in long. dir. to ensure incompressibility."""
        return 1.0 - self.k1**2 / (self.k1**2 + self.k2**2)
    def p_y(self):
        """Projector in trans. dir. to ensure incompressibility."""
        return -self.k1 * self.k2 / (self.k1**2 + self.k2**2)


class Filtered_Incompr_Field(Incompr_Field):
    """Generates a 2d random filtered incompressible velocity field.

        Use Field.U_filtered(x, y) to obtain the velocity vector(s) (u, v) or
        Field.u_filtered(x, y) and Field.v_filtered(x, y) to obtain the
        components individually.
        The second method should be slightly slower (though not tested).

        A coarse grained Kraichnan algorithm together with a projector to
        ensure incompressiblity is used to calculate the field.
        At the moment only Gaussian auto-correlation functions are supported,
        but this could change in the future.
        The field can either be generated
        -on a regular grid,
        -for a given tuple of positions (e.g. for particles)
        -or for single positions.
        This behaviour can be changed by setting the "mode" to
        "grid", "tuple", or "single". "grid" is the standard setting.

        If you use the K-methods for calculating the conductivities via this
        class, you have to set the value for K_G after instantiating. Otherwise
        a value of K_G = 1 is assumed.

        Args:
            mean_velocity (float): mean velocity in e1-direction
            variance (float): variance of the conductivity field
            corr_len (float, array_like): correlation length of the
                conductivity field, if a single value is passed,
                isotropy is assumed
            filter_width (float): the size of the filter
            no_random_modes (int): number of random modes to approximate the
                 Fourier sum
            seed (int, optional): set the seed, if "None",
                a random seed is used

        Examples:
        >>> f = Filtered_Incompr_Field(1, 1, [1, 1], 4, 100, 609006)
        >>> f.mode = 'single'
        >>> np.around(f.U_filtered(0, 0), 4)
        array([ 1.6522,  0.1172])
        >>> round(f.u_filtered(0, 0), 4)
        1.6522
        >>> round(f.v_filtered(0, 0), 4)
        0.1172
        >>> np.around(f.U(0, 0), 4)
        array([ 2.4872, -0.118 ])
        >>> round(f.u(0, 0), 4)
        2.4872
        >>> round(f.v(0, 0), 4)
        -0.118
        >>> round(f.K(0,0), 4)
        0.04
        >>> x_tuple = [ 4, 0, 3]
        >>> y_tuple = [-1, 0, 1]
        >>> f.mode = 'tuple'
        >>> np.around(f.U_filtered(x_tuple, y_tuple)[0], 4)
        array([ 1.0285,  1.6522,  1.5197])
        >>> round(f.u_filtered(x_tuple, y_tuple)[1], 4)
        1.6522
        >>> round(f.v_filtered(x_tuple, y_tuple)[1], 4)
        0.1172
        >>> round(f.U(x_tuple, y_tuple)[0][1], 4)
        2.4872
        >>> round(f.U(x_tuple, y_tuple)[1][1], 4)
        -0.118
        >>> round(f.u(x_tuple, y_tuple)[1], 4)
        2.4872
        >>> round(f.v(x_tuple, y_tuple)[1], 4)
        -0.118
        >>> round(f.K(x_tuple, y_tuple)[1], 4)
        0.04
        >>> x_grid = np.arange(0, 5, 0.5)
        >>> y_grid = np.arange(0, 5, 1)
        >>> f.mode = 'grid'
        >>> round(f.U_filtered(x_grid, y_grid)[0][0,0], 4)
        1.6522
        >>> round(f.U_filtered(x_grid, y_grid)[1][0,0], 4)
        0.1172
        >>> round(f.u_filtered(x_grid, y_grid)[0,0], 4)
        1.6522
        >>> round(f.v_filtered(x_grid, y_grid)[0,0], 4)
        0.1172
        >>> round(f.U(x_grid, y_grid)[0][0,0], 4)
        2.4872
        >>> round(f.U(x_grid, y_grid)[1][0,0], 4)
        -0.118
        >>> round(f.u(x_grid, y_grid)[0,0], 4)
        2.4872
        >>> round(f.v(x_grid, y_grid)[0,0], 4)
        -0.118
        >>> round(f.K(x_grid, y_grid)[0,0], 4)
        0.04

        #testing reshaping in all combinations
        >>> f.mode = 'grid'
        >>> f.mode = 'single'
        >>> f.mode = 'grid'
        >>> f.mode = 'tuple'
        >>> f.mode = 'grid'
        >>> f.mode = 'single'
        >>> f.mode = 'tuple'
        >>> f.mode = 'single'

        License:
        This class is NOT released under the GNU Lesser General Public License,
        yet.

        If you use this class, please contact Lennart Schüler.

        Copyright 2015, Lennart Schüler
    """
    def __init__(self, mean_velocity, variance, corr_len, filter_width,
                 no_random_modes, seed=None):
        #only isotropic filter supported at the moment
        #TODO check derivation of filtered equations
        self.L = filter_width
        self.U_filtered = self.U_filtered_grid
        self.u_filtered = self.u_filtered_grid
        self.v_filtered = self.v_filtered_grid
        super(Filtered_Incompr_Field, self).__init__(mean_velocity, variance,
                                                     corr_len, no_random_modes,
                                                     seed)

    def reset(self, seed=None):
        """Resets the field with a new seed.

        Needs to be completely reimplemented, otherwise the base-class calls
        inherented mode-setter and tries to reshape self.gamma,
        which does not exist yet. On the other hand, self.gamma depends
        on self.k1 and self.k2, which are created by super().reset and thus do
        not exist at time of creation of self.gamma.

        Args:
            seed (int, optional): master seed for different RNG streams
        """
        r1, r2, r3 = self.init_RNG(seed)

        self.k1 = r1.normal(0., 1./self.l[0], self.N)
        self.k2 = r2.normal(0., 1./self.l[1], self.N)
        self.phi = r3.uniform(0., 2*np.pi, self.N)
        self.gamma = (np.sin(self.k1 * 0.5 * self.L) *
                      np.sin(self.k2 * 0.5 * self.L) /
                      (self.k1 * self.k2))
        #reshape according to mode
        self.mode = self._mode

    #TODO can super be used here!?
    @Field.mode.setter
    def mode(self, m):
        """Sets the mode (grid, tuple, single).

            Needs to be overridden, because the self.U_filtered,
            self.u_filtered and self.v_filtered attributes need to be adapted
            to the given mode.

        Args:
            m (string): mode of calculation
        """
        if m == 'grid':
            self.U = self.U_grid
            self.u = self.u_grid
            self.v = self.v_grid
            self.U_filtered = self.U_filtered_grid
            self.u_filtered = self.u_filtered_grid
            self.v_filtered = self.v_filtered_grid
            self.K = self.K_grid
            self.k1 = np.squeeze(self.k1)
            self.k2 = np.squeeze(self.k2)
            self.phi = np.squeeze(self.phi)
            self.gamma = np.squeeze(self.gamma)
            self.k1 = np.reshape(self.k1, (1, 1, len(self.k1)))
            self.k2 = np.reshape(self.k2, (1, 1, len(self.k2)))
            self.phi = np.reshape(self.phi, (1, 1, len(self.phi)))
            self.gamma = np.reshape(self.gamma, (1, 1, len(self.gamma)))
        elif m == 'tuple':
            self.U = self.U_tuple
            self.u = self.u_tuple
            self.v = self.v_tuple
            self.U_filtered = self.U_filtered_tuple
            self.u_filtered = self.u_filtered_tuple
            self.v_filtered = self.v_filtered_tuple
            self.K = self.K_tuple
            self.k1 = np.squeeze(self.k1)
            self.k2 = np.squeeze(self.k2)
            self.phi = np.squeeze(self.phi)
            self.gamma = np.squeeze(self.gamma)
            self.k1 = np.reshape(self.k1, (1, len(self.k1)))
            self.k2 = np.reshape(self.k2, (1, len(self.k2)))
            self.phi = np.reshape(self.phi, (1, len(self.phi)))
            self.gamma = np.reshape(self.gamma, (1, len(self.gamma)))
        elif m == 'single':
            self.U = self.U_single
            self.u = self.u_single
            self.v = self.v_single
            self.U_filtered = self.U_filtered_single
            self.u_filtered = self.u_filtered_single
            self.v_filtered = self.v_filtered_single
            self.K = self.K_single
            try:
                self.k1 = np.squeeze(self.k1, axis=1)
                self.k2 = np.squeeze(self.k2, axis=1)
                self.phi = np.squeeze(self.phi, axis=1)
                self.gamma = np.squeeze(self.gamma, axis=1)
            except:
                pass
        else:
            raise ValueError('Unknown mode {}'.format(m))
        self._mode = m

    def U_filtered_grid(self, x, y):
        """Calculate the filtered velocity vector for grid (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            (u, v) (2d array_like, 2d array_like): filtered velocity vector
                on grid (x, y)
        """
        x, y = self.reshape_axis_to_grid(x, y)
        u = (self.u_bar - self.sigma * self.u_bar / self.L**2 *
                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))
        v = (-self.sigma * self.u_bar / self.L**2 *
            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))
        return (u, v)

    def u_filtered_grid(self, x, y):
        """Calculate the filtered long. velocity component for grid (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            u (2d array_like): filtered velocity in e1-direction on grid (x, y)
        """
        x, y = self.reshape_axis_to_grid(x, y)
        return (self.u_bar - self.sigma * self.u_bar / self.L**2 *
                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))

    def v_filtered_grid(self, x, y):
        """Calculate the filtered trans. velocity component for grid (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            u (2d array_like): filtered velocity in e2-direction on grid (x, y)
        """
        x, y = self.reshape_axis_to_grid(x, y)
        return (-self.sigma * self.u_bar / self.L**2 *
            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))


    def U_filtered_single(self, x, y):
        """Calculate the filtered velocity vector for point (x, y).

        Args:
            x (float): x position
            y (float): y position
        Returns:
            (u, v) (float, float): filtered velocity vector at point (x, y)
        """
        u = float(self.u_bar - self.sigma * self.u_bar / self.L**2 *
                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))
        v = float(-self.sigma * self.u_bar / self.L**2 *
            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))
        return (u, v)

    def u_filtered_single(self, x, y):
        """Calculate the filtered long. velocity component for point (x, y).

        Args:
            x (float): x position
            y (float): y position
        Returns:
            u (float): velocity in e1-direction at point (x, y)
        """
        return float(self.u_bar - self.sigma * self.u_bar / self.L**2 *
                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))

    def v_filtered_single(self, x, y):
        """Calculate the filtered trans. velocity component for point (x, y).

        Args:
            x (float): x position
            y (float): y position
        Returns:
            u (float): velocity in e2-direction at point (x, y)
        """
        return float(-self.sigma * self.u_bar / self.L**2 *
            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))

    def U_filtered_tuple(self, x, y):
        """Calculate the filtered velocity vectors for tuples (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            (u, v) (array_like, array_like): filtered velocity vectors for
                tuples (x, y)
        """
        x, y = self.reshape_axis_to_tuple(x, y)
        u = (self.u_bar - self.sigma * self.u_bar / self.L**2 *
                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))
        v = (-self.sigma * self.u_bar / self.L**2 *
            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))
        return (u, v)

    def u_filtered_tuple(self, x, y):
        """Calculate the filtered long. velocity component for tuples (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            u (array_like): velocity in e1-direction for tuples (x, y)
        """
        x, y = self.reshape_axis_to_tuple(x, y)
        return (self.u_bar - self.sigma * self.u_bar / self.L**2 *
                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))

    def v_filtered_tuple(self, x, y):
        """Calculate the filtered trans. velocity component for tuples (x, y).

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            u (array_like): velocity in e2-direction for tuples (x, y)
        """
        x, y = self.reshape_axis_to_tuple(x, y)
        return (-self.sigma * self.u_bar / self.L**2 *
            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))

    def f_incompr_filtered_x(self, x, y):
        """Calculate the filtered (long.) Fourier sum for the incompr. field.

        Args:
            x (float or array_like): x position
            y (float or array_like): y position
        """
        return np.sum(self.p_x() * self.gamma *
                      np.cos(self.k1*x + self.k2*y + self.phi), axis=-1)

    def f_incompr_filtered_y(self, x, y):
        """Calculate the filtered (trans.) Fourier sum for the incompr. field.

        Args:
            x (float or array_like): x position
            y (float or array_like): y position
        """
        return np.sum(self.p_y() * self.gamma *
                      np.cos(self.k1*x + self.k2*y + self.phi), axis=-1)

if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # f = Field(1, 1, [1, 1], 100, 15011997)
    # f.mode = 'single'
    # round(f.K(0,0), 4)
