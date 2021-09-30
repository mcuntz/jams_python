#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

import numpy as np
import numpy.random as rand
# import sys


"""Note: The different classes are published under different licenses."""


class RNG(object):
    """Wrapper class for random number generators with a consistent interface.

        This class can be used as a base class to pass along a common
        interface for RNGs, if different distributions are needed.

         Args:
             seed (int, optional): set the seed of the master RNG, if "None",
                 a random seed is used

         Examples:
             >>> rng = RNG(19031977)

             np.around(rng('gau', [1, 1]), 4)
             array([ 0.8996,  0.9699, -0.3083, -1.7153])
             >>> rng.seed = 28091977

             np.around(rng('gau', [1, 1], 2), 4)
             array([[-1.7336,  0.3146],
                    [ 0.284 ,  0.3108],
                    [-0.449 ,  2.1735],
                    [ 0.1099, -1.5756]])

         License:
         --------
         This file is part of the JAMS Python package, distributed under the MIT
         License. The JAMS Python package originates from the former UFZ Python library,
         Department of Computational Hydrosystems, Helmholtz Centre for Environmental
         Research - UFZ, Leipzig, Germany.

         Copyright (c) 2015 Lennart Schueler

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
    """
    def __init__(self, seed=None):
        self.seed = seed

    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, new_seed=None):
        self._seed = new_seed
        self._master_RNG_fct = rand.RandomState(new_seed)
        # random_integers is deprecated
        # self._master_RNG = \
        #     lambda: self._master_RNG_fct.random_integers(2**32-1)
        self._master_RNG = \
            lambda: self._master_RNG_fct.randint(1, 4294967295 + 1)

    def gau(self, corr_len, size=None):
        r1 = rand.RandomState(self._master_RNG())
        r2 = rand.RandomState(self._master_RNG())
        k1 = r1.normal(0., 1/corr_len[0], size)
        k2 = r2.normal(0., 1/corr_len[0], size)
        return k1, k2

    def exp(self, corr_len, size=None):
        raise NotImplementedError('Exponential correlation function not '
                                  'yet implemented')

    def __call__(self, corr_fct, corr_len, size=None):
        try:
            k1, k2 = getattr(self, corr_fct)(corr_len, size)
        except AttributeError:
            raise ValueError('Unknown correlation function {0}'.format(corr_fct))
        r1 = rand.RandomState(self._master_RNG())
        r2 = rand.RandomState(self._master_RNG())
        Z1 = r1.normal(size=size)
        Z2 = r2.normal(size=size)
        return Z1, Z2, k1, k2


class RandMeth(object):
    """Randomization method for calculating standard random fields.

    Examples:
        >>> rm = RandMeth('gau', [1, 1], seed=12091986)
        >>> rm.mode = 'single'

        round(rm.Y(0,0), 4)
        3.0624
        >>> x_tuple = [ 4, 0, 3]
        >>> y_tuple = [-1, 0, 1]
        >>> rm.mode = 'unstructured'

        round(rm.Y(x_tuple, y_tuple)[1], 4)
        3.0624
        >>> rm.Y(x_tuple, y_tuple).shape
        (3,)
        >>> x_grid = np.arange(0, 5, 0.5)
        >>> y_grid = np.arange(0, 5, 1)
        >>> rm.mode = 'structured'

        round(rm.Y(x_grid, y_grid)[0,0], 4)
        3.0624
        >>> rm.Y(x_tuple, y_tuple).shape
        (3, 3)

        #testing reshaping in all combinations
        >>> rm.mode = 'structured'
        >>> rm.mode = 'single'
        >>> rm.mode = 'structured'
        >>> rm.mode = 'unstructured'
        >>> rm.mode = 'structured'
        >>> rm.mode = 'single'
        >>> rm.mode = 'unstructured'
        >>> rm.mode = 'single'

    """
    def __init__(self, corr_fct, corr_len, random_modes_no=100, seed=None,
                 mode='structured', rng_class=RNG):
        self.corr_fct = corr_fct

        #correlation length
        #assuming isotropy, if scalar correlation length is given
        if not hasattr(corr_len, '__getitem__'):
            self.l = np.array((corr_len, corr_len))
        else:
            self.l = np.array(corr_len)

        self.N = random_modes_no
        self._mode = mode
        self.rng = rng_class(seed)
        self.reset()

    def reset(self, seed=None):
        self.Z1, self.Z2, self.k1, self.k2 = self.rng(self.corr_fct, self.l,
                                                      size=self.N)
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
        # random_integers is deprecated
        # self.master_RNG = \
        #     lambda: self.__master_RNG.random_integers(4294967295)
        self._master_RNG = \
            lambda: self._master_RNG_fct.randint(1, 4294967295 + 1)
        r1 = r.RandomState(self.master_RNG())
        r2 = r.RandomState(self.master_RNG())
        r3 = r.RandomState(self.master_RNG())
        return (r1, r2, r3)

    @property
    def mode(self):
        """Returns the mode (grid, tuple, single)."""
        return self._mode

    def _set_mode(self, m):
        """Sets the mode (structured, unstructured, single).

        Args:
            m (string): mode of calculation
        """
        if m == 'structured':
            self.Y = self.Y_structured
            self.k1 = np.squeeze(self.k1)
            self.k2 = np.squeeze(self.k2)
            self.k1 = np.reshape(self.k1, (1, 1, len(self.k1)))
            self.k2 = np.reshape(self.k2, (1, 1, len(self.k2)))
        elif m == 'unstructured':
            self.Y = self.Y_unstructured
            self.k1 = np.squeeze(self.k1)
            self.k2 = np.squeeze(self.k2)
            self.k1 = np.reshape(self.k1, (1, len(self.k1)))
            self.k2 = np.reshape(self.k2, (1, len(self.k2)))
        elif m == 'single':
            self.Y = self.Y_single
            try:
                self.k1 = np.squeeze(self.k1, axis=1)
                self.k2 = np.squeeze(self.k2, axis=1)
            except:
                pass
        else:
            raise ValueError('Unknown mode {}'.format(m))
        self._mode = m

    def reshape_axis_to_structured(self, x, y):
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

    def reshape_axis_to_unstructured(self, x, y):
        """Reshape given positions for vectorisation on unstructured grid.

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

    def calc_Y(self, x, y):
        return (np.sqrt(2./self.N) *
                np.sum(self.Z1 * np.cos(2*np.pi*self.k1*x) +
                       self.Z2 * np.sin(2*np.pi*self.k2*y), axis=-1))

    def Y_single(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (float): log-hydraulic conductivity
        """
        return float(self.calc_Y(x, y))

    def Y_unstructured(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_unstructured(x, y)
        return self.calc_Y(x, y)

    def Y_structured(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (2d array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_structured(x, y)
        return self.calc_Y(x, y)

    mode = property(fget = lambda self:        self._get_mode(),
                    fset = lambda self, value: self._set_mode(value))


class IncomprRandMeth(RandMeth):
    """Randomization method for calculating standard random fields.

        Examples:
        >>> rm = IncomprRandMeth('gau', [1, 1], 100, 16052001)

        #>>> rm.mode = 'single'
        #>>> np.around(rm.incompr_Y(0, 0), 4)
        #array([-0.1997,  0.3987])

#        >>> round(f.u(0, 0), 4)
#        -0.1997
#        >>> round(f.v(0, 0), 4)
#        0.3987
#        >>> round(f.K(0,0), 4)
#        4.6852
#        >>> x_tuple = [ 4, 0, 3]
#        >>> y_tuple = [-1, 0, 1]
#        >>> f.mode = 'unstructured'
#        >>> np.around(f.U(x_tuple, y_tuple)[0], 4)
#        array([ 0.0112, -0.1997,  0.9485])
#        >>> round(f.u(x_tuple, y_tuple)[1], 4)
#        -0.1997
#        >>> round(f.v(x_tuple, y_tuple)[1], 4)
#        0.3987
#        >>> round(f.K(x_tuple, y_tuple)[1], 4)
#        4.6852
#        >>> x_grid = np.arange(0, 5, 0.5)
#        >>> y_grid = np.arange(0, 5, 1)
#        >>> f.mode = 'structured'
#        >>> round(f.U(x_grid, y_grid)[0][0,0], 4)
#        -0.1997
#        >>> round(f.U(x_grid, y_grid)[1][0,0], 4)
#        0.3987
#        >>> round(f.u(x_grid, y_grid)[0,0], 4)
#        -0.1997
#        >>> round(f.v(x_grid, y_grid)[0,0], 4)
#        0.3987
#        >>> round(f.K(x_grid, y_grid)[0,0], 4)
#        4.6852
#
#        #testing reshaping in all combinations
#        >>> f.mode = 'structured'
#        >>> f.mode = 'single'
#        >>> f.mode = 'structured'
#        >>> f.mode = 'unstructured'
#        >>> f.mode = 'structured'
#        >>> f.mode = 'single'
#        >>> f.mode = 'unstructured'
#        >>> f.mode = 'single'

    """
    def __init__(self, corr_fct, corr_len, random_modes_no=100, seed=None,
                 mode='structured', rng_class=RNG):
        super(IncomprRandMeth, self).__init__(corr_fct, corr_len,
                                              random_modes_no, seed, mode,
                                              rng_class)
        self.mode = mode

    def _set_mode(self, m):
        """Sets the mode (structured, unstructured, single).

        Args:
            m (string): mode of calculation
        """
        if m == 'structured':
            self.incompr_Y = self.incompr_Y_structured
            self.incompr_Y_x = self.incompr_Y_structured_x
            self.incompr_Y_y = self.incompr_Y_structured_y
        elif m == 'unstructured':
            self.incompr_Y = self.incompr_Y_unstructured
            self.incompr_Y_x = self.incompr_Y_unstructured_x
            self.incompr_Y_y = self.incompr_Y_unstructured_y
        elif m == 'single':
            self.incompr_Y = self.incompr_Y_single
            self.incompr_Y_x = self.incompr_Y_single_x
            self.incompr_Y_y = self.incompr_Y_single_y
        else:
            raise ValueError('Unknown mode {}'.format(m))
        super(IncomprRandMeth, self)._set_mode(m)
        #super(IncomprRandMeth, self.__class__).mode.fset(self, m)

    def calc_incompr_Y_x(self, x, y):
        """Calculate the (long.) Fourier sum for the incompr. random field.

        Args:
            x (float or array_like): x position
            y (float or array_like): y position
        """
        return (np.sqrt(2./self.N) * np.sum(self.p_x() *
                   (self.Z1 * np.cos(2*np.pi*self.k1*x) +
                    self.Z2 * np.sin(2*np.pi*self.k2*y)), axis=-1))

    def calc_incompr_Y_y(self, x, y):
        """Calculate the (trans.) Fourier sum for the incompr. random field.

        Args:
            x (float or array_like): x position
            y (float or array_like): y position
        """
        return (np.sqrt(2./self.N) * np.sum(self.p_y() *
                   (self.Z1 * np.cos(2*np.pi*self.k1*x) +
                    self.Z2 * np.sin(2*np.pi*self.k2*y)), axis=-1))

    def incompr_Y_single_x(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (float): log-hydraulic conductivity
        """
        return float(self.calc_incompr_Y_x(x, y))
    def incompr_Y_single_y(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (float): log-hydraulic conductivity
        """
        return float(self.calc_incompr_Y_y(x, y))
    def incompr_Y_single(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (float): log-hydraulic conductivity
        """
        return self.incompr_Y_single_x, self.incompr_Y_single_y

    def incompr_Y_unstructured(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_unstructured(x, y)
        return self.calc_incompr_Y_x(x, y), self.calc_incompr_Y_y(x, y)
    def incompr_Y_unstructured_x(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_unstructured(x, y)
        return self.calc_incompr_Y_x(x, y)
    def incompr_Y_unstructured_y(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_unstructured(x, y)
        return self.calc_incompr_Y_y(x, y)

    def incompr_Y_structured(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (2d array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_structured(x, y)
        return self.calc_incompr_Y_x(x, y), self.calc_incompr_Y_y(x, y)
    def incompr_Y_structured_x(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (2d array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_structured(x, y)
        return self.calc_incompr_Y_x(x, y)
    def incompr_Y_structured_y(self, x, y):
        """Calculate the standard random field.

        Args:
            x (array_like): x position
            y (array_like): y position
        Returns:
            K (2d array_like): log-hydraulic conductivity
        """
        x, y = self.reshape_axis_to_structured(x, y)
        return self.calc_incompr_Y_y(x, y)

    def p_x(self):
        """Projector in long. dir. to ensure incompressibility."""
        return 1.0 - self.k1**2 / (self.k1**2 + self.k2**2)
    def p_y(self):
        """Projector in trans. dir. to ensure incompressibility."""
        return -self.k1 * self.k2 / (self.k1**2 + self.k2**2)


class Field(object):
    """Generates a 2d random hydraulic conductivity field.

        Use Field.K(x, y) for calculating the field.

        A Randomization algorithm (Hesse et al. 2014, Environ. Model. Softw.) to
        calculate the field is used.
        At the moment only a Gaussian auto-correlation function is used,
        but this could change in the future.
        The field can either be generated
        -on a regular grid,
        -for a given tuple of positions (e.g. for particles)
        -or for single positions.
        This behaviour can be changed by setting the "mode" to
        "structured", "unstructured", or "single". "structured" is the standard setting.

        Args:
            K_G (float): geometric mean of hydraulic conductivity (exp(mu))
            variance (float): variance of the conductivity field
            corr_len (float or array_like): correlation length of the
                conductivity field, if a single value is passed,
                isotropy is assumed
            random_modes_no (int): number of random modes to approximate the
                 Fourier sum
            seed (int, optional): set the seed of the master RNG, if "None",
                a random seed is used

        Examples:
        >>> f = Field(1, 1, [1, 1], 100, 15011997)
        >>> f.mode = 'single'

        round(f.K(0,0), 4)
        2.1771
        >>> x_tuple = [ 4, 0, 3]
        >>> y_tuple = [-1, 0, 1]
        >>> f.mode = 'unstructured'

        round(f.K(x_tuple, y_tuple)[1], 4)
        2.1771
        >>> f.K(x_tuple, y_tuple).shape
        (3,)
        >>> x_grid = np.arange(0, 5, 0.5)
        >>> y_grid = np.arange(0, 5, 1)
        >>> f.mode = 'structured'

        round(f.K(x_grid, y_grid)[0,0], 4)
        2.1771
        >>> f.K(x_tuple, y_tuple).shape
        (3, 3)

        #testing reshaping in all combinations
        >>> f.mode = 'structured'
        >>> f.mode = 'single'
        >>> f.mode = 'structured'
        >>> f.mode = 'unstructured'
        >>> f.mode = 'structured'
        >>> f.mode = 'single'
        >>> f.mode = 'unstructured'
        >>> f.mode = 'single'

        License:
        --------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2015 Lennart Schueler

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
    """
    def __init__(self, K_G, variance, corr_len, random_modes_no,
                 seed=None):
        #geometric mean
        self.K_G = K_G
        #standard deviation
        self.sigma = np.sqrt(variance)

        self.rand_meth = RandMeth('gau', corr_len, random_modes_no, seed)

    def _get_mode(self):
        return self.rand_meth.mode

    def _set_mode(self, m):
        self.rand_meth.mode = m

    def reset(self, seed=None):
        self.rand_meth.reset(seed)

    def K(self, x, y):
        return self.K_G * np.exp(self.sigma * self.rand_meth.Y(x, y))

    mode = property(fget = lambda self:        self._get_mode(),
                    fset = lambda self, value: self._set_mode(value))

class IncomprField(Field):
    """Generates a 2d random incompressible velocity field.

        Use Field.U(x, y) to obtain the velocity vector(s) (u, v) or
        Field.u(x, y) and Field.v(x, y) to obtain the components individually.
        The second method should be slightly slower (though not tested).

        A Randomization algorithm (Hesse et al. 2014, Environ. Model. Softw.) to
        calculate the field is used.
        At the moment only a Gaussian auto-correlation function is used,
        but this could change in the future.
        The field can either be generated
        -on a regular grid,
        -for a given tuple of positions (e.g. for particles)
        -or for single positions.
        This behaviour can be changed by setting the "mode" to
        "structured", "unstructured", or "single". "structured" is the standard setting.

        If you use the K-methods for calculating the conductivities via this
        class, you have to set the value for K_G after instantiating. Otherwise
        a value of K_G = 1 is assumed.

        Args:
            mean_velocity (float): mean velocity in e1-direction
            variance (float): variance of the conductivity field
            corr_len (float or array_like): correlation length of the
                conductivity field, if a single value is passed,
                isotropy is assumed
            random_modes_no (int): number of random modes to approximate the
                 Fourier sum
            seed (int, optional): set the seed, if "None",
                a random seed is used

        Examples:
#        >>> f = IncomprField(1, 1, [1, 1], 100, 16052001)
#        >>> f.mode = 'single'
#        >>> np.round(f.U(0, 0), 4)
#        array([-0.1997,  0.3987])
#
#        >>> round(f.u(0, 0), 4)
#        -0.1997
#        >>> round(f.v(0, 0), 4)
#        0.3987
#        >>> round(f.K(0,0), 4)
#        4.6852
#        >>> x_tuple = [ 4, 0, 3]
#        >>> y_tuple = [-1, 0, 1]
#        >>> f.mode = 'unstructured'
#        >>> np.around(f.U(x_tuple, y_tuple)[0], 4)
#        array([ 0.0112, -0.1997,  0.9485])
#        >>> round(f.u(x_tuple, y_tuple)[1], 4)
#        -0.1997
#        >>> round(f.v(x_tuple, y_tuple)[1], 4)
#        0.3987
#        >>> round(f.K(x_tuple, y_tuple)[1], 4)
#        4.6852
#        >>> x_grid = np.arange(0, 5, 0.5)
#        >>> y_grid = np.arange(0, 5, 1)
#        >>> f.mode = 'structured'
#        >>> round(f.U(x_grid, y_grid)[0][0,0], 4)
#        -0.1997
#        >>> round(f.U(x_grid, y_grid)[1][0,0], 4)
#        0.3987
#        >>> round(f.u(x_grid, y_grid)[0,0], 4)
#        -0.1997
#        >>> round(f.v(x_grid, y_grid)[0,0], 4)
#        0.3987
#        >>> round(f.K(x_grid, y_grid)[0,0], 4)
#        4.6852
#
#        #testing reshaping in all combinations
#        >>> f.mode = 'structured'
#        >>> f.mode = 'single'
#        >>> f.mode = 'structured'
#        >>> f.mode = 'unstructured'
#        >>> f.mode = 'structured'
#        >>> f.mode = 'single'
#        >>> f.mode = 'unstructured'
#        >>> f.mode = 'single'

        License:
        --------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2015 Lennart Schueler

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
    """
    def __init__(self, mean_velocity, variance, corr_len, random_modes_no,
                 seed=None):
        self.u_bar = mean_velocity
        self.incompr_rand_meth = IncomprRandMeth('gau', corr_len,
                                                 random_modes_no, seed)
        super(IncomprField, self).__init__(mean_velocity, variance, corr_len,
                                           random_modes_no, seed)

    def _get_mode(self):
        #TODO call self.rand_meth.mode?
        return self.incompr_rand_meth.mode

    def _set_mode(self, m):
        #TODO call self.rand_meth.mode?
        self.incompr_rand_meth.mode = m

    def u(self, x, y):
        return (self.u_bar - self.sigma * self.u_bar *
                self.incompr_rand_meth.incompr_Y_x(x, y))
    def v(self, x, y):
        return (-self.sigma * self.u_bar *
                self.incompr_rand_meth.incompr_Y_y(x, y))
    def U(self, x, y):
        return self.u(x, y), self.v(x, y)


class FilteredIncomprField(IncomprField):
    def __init__(self, mean_velocity, variance, corr_len, filter_size,
                 random_modes_no, seed=None):
        self.filter_size = filter_size
        self.filtered_incompr_rand_meth = FilteredIncomprRandMeth('gau',
            corr_len, filter_size, random_modes_no, seed)
        super(FilteredIncomprField, self).__init__(mean_velocity, variance,
                                                   corr_len, random_modes_no,
                                                   seed)
    def _get_mode(self):
        #TODO call self.rand_meth.mode?
        return self.filtered_rand_meth.mode

    def _set_mode(self, m):
        #TODO call self.rand_meth.mode?
        self.filtered_rand_meth.mode = m

    def u_filtered(self, x, y):
        return (self.u_bar - self.sigma * self.u_bar *
                self.filtered_incompr_rand_meth.filtered_incompr_Y_x(x, y))
    def v_filtered(self, x, y):
        return (-self.sigma * self.u_bar *
                self.filtered_incompr_rand_meth.filtered_incompr_Y_y(x, y))
    def U_filtered(self, x, y):
        return self.u(x, y), self.v(x, y)


#class FilteredIncomprField(IncomprField):
#    """Generates a 2d random filtered incompressible velocity field.
#
#        Use Field.U_filtered(x, y) to obtain the velocity vector(s) (u, v) or
#        Field.u_filtered(x, y) and Field.v_filtered(x, y) to obtain the
#        components individually.
#        The second method should be slightly slower (though not tested).
#
#        A coarse grained Randomization algorithm
#        (Hesse et al. 2014, Environ. Model. Softw.) together with a projector to
#        ensure incompressiblity is used to calculate the field.
#        At the moment only Gaussian auto-correlation functions are supported,
#        but this could change in the future.
#        The field can either be generated
#        -on a regular grid,
#        -for a given tuple of positions (e.g. for particles)
#        -or for single positions.
#        This behaviour can be changed by setting the "mode" to
#        "structured", "unstructured", or "single". "structured" is the standard setting.
#
#        If you use the K-methods for calculating the conductivities via this
#        class, you have to set the value for K_G after instantiating. Otherwise
#        a value of K_G = 1 is assumed.
#
#        Args:
#            mean_velocity (float): mean velocity in e1-direction
#            variance (float): variance of the conductivity field
#            corr_len (float, array_like): correlation length of the
#                conductivity field, if a single value is passed,
#                isotropy is assumed
#            filter_width (float): the size of the filter
#            random_modes_no (int): number of random modes to approximate the
#                 Fourier sum
#            seed (int, optional): set the seed, if "None",
#                a random seed is used
#
#        Examples:
#        >>> f = FilteredIncomprField(1, 1, [1, 1], 4, 100, 609006)
#        >>> f.mode = 'single'
#        >>> np.around(f.U_filtered(0, 0), 4)
#        array([ 1.6522,  0.1172])
#        >>> round(f.u_filtered(0, 0), 4)
#        1.6522
#        >>> round(f.v_filtered(0, 0), 4)
#        0.1172
#        >>> np.around(f.U(0, 0), 4)
#        array([ 2.4872, -0.118 ])
#        >>> round(f.u(0, 0), 4)
#        2.4872
#        >>> round(f.v(0, 0), 4)
#        -0.118
#        >>> round(f.K(0,0), 4)
#        0.04
#        >>> x_tuple = [ 4, 0, 3]
#        >>> y_tuple = [-1, 0, 1]
#        >>> f.mode = 'tuple'
#        >>> np.around(f.U_filtered(x_tuple, y_tuple)[0], 4)
#        array([ 1.0285,  1.6522,  1.5197])
#        >>> round(f.u_filtered(x_tuple, y_tuple)[1], 4)
#        1.6522
#        >>> round(f.v_filtered(x_tuple, y_tuple)[1], 4)
#        0.1172
#        >>> round(f.U(x_tuple, y_tuple)[0][1], 4)
#        2.4872
#        >>> round(f.U(x_tuple, y_tuple)[1][1], 4)
#        -0.118
#        >>> round(f.u(x_tuple, y_tuple)[1], 4)
#        2.4872
#        >>> round(f.v(x_tuple, y_tuple)[1], 4)
#        -0.118
#        >>> round(f.K(x_tuple, y_tuple)[1], 4)
#        0.04
#        >>> x_grid = np.arange(0, 5, 0.5)
#        >>> y_grid = np.arange(0, 5, 1)
#        >>> f.mode = 'structured'
#        >>> round(f.U_filtered(x_grid, y_grid)[0][0,0], 4)
#        1.6522
#        >>> round(f.U_filtered(x_grid, y_grid)[1][0,0], 4)
#        0.1172
#        >>> round(f.u_filtered(x_grid, y_grid)[0,0], 4)
#        1.6522
#        >>> round(f.v_filtered(x_grid, y_grid)[0,0], 4)
#        0.1172
#        >>> round(f.U(x_grid, y_grid)[0][0,0], 4)
#        2.4872
#        >>> round(f.U(x_grid, y_grid)[1][0,0], 4)
#        -0.118
#        >>> round(f.u(x_grid, y_grid)[0,0], 4)
#        2.4872
#        >>> round(f.v(x_grid, y_grid)[0,0], 4)
#        -0.118
#        >>> round(f.K(x_grid, y_grid)[0,0], 4)
#        0.04
#
#        #testing reshaping in all combinations
#        >>> f.mode = 'structured'
#        >>> f.mode = 'single'
#        >>> f.mode = 'structured'
#        >>> f.mode = 'tuple'
#        >>> f.mode = 'structured'
#        >>> f.mode = 'single'
#        >>> f.mode = 'tuple'
#        >>> f.mode = 'single'
#
#        License:
#        This class is NOT released under the GNU Lesser General Public License,
#        yet.
#
#        If you use this class, please contact Lennart Schueler.
#
#        Copyright 2015, Lennart Schueler
#    """
#    def __init__(self, mean_velocity, variance, corr_len, filter_width,
#                 random_modes_no, seed=None):
#        #only isotropic filter supported at the moment
#        #TODO check derivation of filtered equations
#        self.L = filter_width
#        self.U_filtered = self.U_filtered_structured
#        self.u_filtered = self.u_filtered_structured
#        self.v_filtered = self.v_filtered_structured
#        super(FilteredIncomprField, self).__init__(mean_velocity, variance,
#                                                   corr_len, random_modes_no,
#                                                   seed)
#
#    def reset(self, seed=None):
#        """Resets the field with a new seed.
#
#        Needs to be completely reimplemented, otherwise the base-class calls
#        inherented mode-setter and tries to reshape self.gamma,
#        which does not exist yet. On the other hand, self.gamma depends
#        on self.k1 and self.k2, which are created by super().reset and thus do
#        not exist at time of creation of self.gamma.
#
#        Args:
#            seed (int, optional): master seed for different RNG streams
#        """
#        self.rng = NormalRNG(seed)
#
#        self.k1 = self.rng.normal(0., 1./self.l[0], self.N)
#        self.k2 = self.rng.normal(0., 1./self.l[1], self.N)
#        self.phi = self.rng.uniform(0., 2*np.pi, self.N)
#        self.gamma = (np.sin(self.k1 * 0.5 * self.L) *
#                      np.sin(self.k2 * 0.5 * self.L) /
#                      (self.k1 * self.k2))
#        #reshape according to mode
#        self.mode = self._mode
#
#    #TODO can super be used here!?
#    @Field.mode.setter
#    def mode(self, m):
#        """Sets the mode (structured, unstructured, single).
#
#            Needs to be overridden, because the self.U_filtered,
#            self.u_filtered and self.v_filtered attributes need to be adapted
#            to the given mode.
#
#        Args:
#            m (string): mode of calculation
#        """
#        if m == 'structured':
#            self.U = self.U_structured
#            self.u = self.u_structured
#            self.v = self.v_structured
#            self.U_filtered = self.U_filtered_structured
#            self.u_filtered = self.u_filtered_structured
#            self.v_filtered = self.v_filtered_structured
#            self.K = self.K_structured
#            self.k1 = np.squeeze(self.k1)
#            self.k2 = np.squeeze(self.k2)
#            self.phi = np.squeeze(self.phi)
#            self.gamma = np.squeeze(self.gamma)
#            self.k1 = np.reshape(self.k1, (1, 1, len(self.k1)))
#            self.k2 = np.reshape(self.k2, (1, 1, len(self.k2)))
#            self.phi = np.reshape(self.phi, (1, 1, len(self.phi)))
#            self.gamma = np.reshape(self.gamma, (1, 1, len(self.gamma)))
#        elif m == 'tuple':
#            self.U = self.U_tuple
#            self.u = self.u_tuple
#            self.v = self.v_tuple
#            self.U_filtered = self.U_filtered_tuple
#            self.u_filtered = self.u_filtered_tuple
#            self.v_filtered = self.v_filtered_tuple
#            self.K = self.K_tuple
#            self.k1 = np.squeeze(self.k1)
#            self.k2 = np.squeeze(self.k2)
#            self.phi = np.squeeze(self.phi)
#            self.gamma = np.squeeze(self.gamma)
#            self.k1 = np.reshape(self.k1, (1, len(self.k1)))
#            self.k2 = np.reshape(self.k2, (1, len(self.k2)))
#            self.phi = np.reshape(self.phi, (1, len(self.phi)))
#            self.gamma = np.reshape(self.gamma, (1, len(self.gamma)))
#        elif m == 'single':
#            self.U = self.U_single
#            self.u = self.u_single
#            self.v = self.v_single
#            self.U_filtered = self.U_filtered_single
#            self.u_filtered = self.u_filtered_single
#            self.v_filtered = self.v_filtered_single
#            self.K = self.K_single
#            try:
#                self.k1 = np.squeeze(self.k1, axis=1)
#                self.k2 = np.squeeze(self.k2, axis=1)
#                self.phi = np.squeeze(self.phi, axis=1)
#                self.gamma = np.squeeze(self.gamma, axis=1)
#            except:
#                pass
#        else:
#            raise ValueError('Unknown mode {}'.format(m))
#        self._mode = m
#
#    def U_filtered_structured(self, x, y):
#        """Calculate the filtered velocity vector for grid (x, y).
#
#        Args:
#            x (array_like): x position
#            y (array_like): y position
#        Returns:
#            (u, v) (2d array_like, 2d array_like): filtered velocity vector
#                on grid (x, y)
#        """
#        x, y = self.reshape_axis_to_structured(x, y)
#        u = (self.u_bar - self.sigma * self.u_bar / self.L**2 *
#                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))
#        v = (-self.sigma * self.u_bar / self.L**2 *
#            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))
#        return (u, v)
#
#    def u_filtered_structured(self, x, y):
#        """Calculate the filtered long. velocity component for grid (x, y).
#
#        Args:
#            x (array_like): x position
#            y (array_like): y position
#        Returns:
#            u (2d array_like): filtered velocity in e1-direction on grid (x, y)
#        """
#        x, y = self.reshape_axis_to_structured(x, y)
#        return (self.u_bar - self.sigma * self.u_bar / self.L**2 *
#                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))
#
#    def v_filtered_structured(self, x, y):
#        """Calculate the filtered trans. velocity component for grid (x, y).
#
#        Args:
#            x (array_like): x position
#            y (array_like): y position
#        Returns:
#            u (2d array_like): filtered velocity in e2-direction on grid (x, y)
#        """
#        x, y = self.reshape_axis_to_structured(x, y)
#        return (-self.sigma * self.u_bar / self.L**2 *
#            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))
#
#
#    def U_filtered_single(self, x, y):
#        """Calculate the filtered velocity vector for point (x, y).
#
#        Args:
#            x (float): x position
#            y (float): y position
#        Returns:
#            (u, v) (float, float): filtered velocity vector at point (x, y)
#        """
#        u = float(self.u_bar - self.sigma * self.u_bar / self.L**2 *
#                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))
#        v = float(-self.sigma * self.u_bar / self.L**2 *
#            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))
#        return (u, v)
#
#    def u_filtered_single(self, x, y):
#        """Calculate the filtered long. velocity component for point (x, y).
#
#        Args:
#            x (float): x position
#            y (float): y position
#        Returns:
#            u (float): velocity in e1-direction at point (x, y)
#        """
#        return float(self.u_bar - self.sigma * self.u_bar / self.L**2 *
#                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))
#
#    def v_filtered_single(self, x, y):
#        """Calculate the filtered trans. velocity component for point (x, y).
#
#        Args:
#            x (float): x position
#            y (float): y position
#        Returns:
#            u (float): velocity in e2-direction at point (x, y)
#        """
#        return float(-self.sigma * self.u_bar / self.L**2 *
#            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))
#
#    def U_filtered_tuple(self, x, y):
#        """Calculate the filtered velocity vectors for tuples (x, y).
#
#        Args:
#            x (array_like): x position
#            y (array_like): y position
#        Returns:
#            (u, v) (array_like, array_like): filtered velocity vectors for
#                tuples (x, y)
#        """
#        x, y = self.reshape_axis_to_tuple(x, y)
#        u = (self.u_bar - self.sigma * self.u_bar / self.L**2 *
#                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))
#        v = (-self.sigma * self.u_bar / self.L**2 *
#            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))
#        return (u, v)
#
#    def u_filtered_tuple(self, x, y):
#        """Calculate the filtered long. velocity component for tuples (x, y).
#
#        Args:
#            x (array_like): x position
#            y (array_like): y position
#        Returns:
#            u (array_like): velocity in e1-direction for tuples (x, y)
#        """
#        x, y = self.reshape_axis_to_tuple(x, y)
#        return (self.u_bar - self.sigma * self.u_bar / self.L**2 *
#                  np.sqrt(32.0/self.N) * self.f_incompr_filtered_x(x, y))
#
#    def v_filtered_tuple(self, x, y):
#        """Calculate the filtered trans. velocity component for tuples (x, y).
#
#        Args:
#            x (array_like): x position
#            y (array_like): y position
#        Returns:
#            u (array_like): velocity in e2-direction for tuples (x, y)
#        """
#        x, y = self.reshape_axis_to_tuple(x, y)
#        return (-self.sigma * self.u_bar / self.L**2 *
#            np.sqrt(32.0/self.N) * self.f_incompr_filtered_y(x, y))
#
#    def f_incompr_filtered_x(self, x, y):
#        """Calculate the filtered (long.) Fourier sum for the incompr. field.
#
#        Args:
#            x (float or array_like): x position
#            y (float or array_like): y position
#        """
#        return np.sum(self.p_x() * self.gamma *
#                      np.cos(self.k1*x + self.k2*y + self.phi), axis=-1)
#
#    def f_incompr_filtered_y(self, x, y):
#        """Calculate the filtered (trans.) Fourier sum for the incompr. field.
#
#        Args:
#            x (float or array_like): x position
#            y (float or array_like): y position
#        """
#        return np.sum(self.p_y() * self.gamma *
#                      np.cos(self.k1*x + self.k2*y + self.phi), axis=-1)

if __name__ == '__main__':
    import doctest
    doctest.testmod()

    # f = Field(1, 1, [1, 1], 100, 15011997)
    # f.mode = 'single'
    # round(f.K(0,0), 4)
