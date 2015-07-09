#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import xlrd

__all__ = ['get_value_excel']

# --------------------------------------------------------------------

def get_value_excel(excelfile, sheet, variable, column):
    """
        Get a cell in a sheet of an Excel file
        where the column is identified by the name in the first row
        and the row is identified by the value in the column "headerout (final)".

        Both, variable and column can be iterables.

        If variable and/or column is None then the list of variables and columns is returned.

        
        Output
        ------
        content of the desired cell
        if column or variable is iterable, output is list
        if column and variable is iterable, output is list (variable) of list (columns)
        if column and/or variable is not given, then a list of all variable and/or columns is returned


        Examples
        --------
        mini = get_value_excel('CHS-measurements.xlsx', 'Logger Hohes Holz', 'BC1_Ptemp [degC]', 'Min')
        -> -20.0
        mini, maxi = get_value_excel('CHS-measurements.xlsx', 'Logger Hohes Holz', 'BC1_Ptemp [degC]', ['Min', 'Max'])
        -> -20.0, 50.0
        min1, min2 = get_value_excel('CHS-measurements.xlsx', 'Logger Hohes Holz', ['BC1_Ptemp [degC]', 'Tair50HMP [degC]'], 'Min')
        -> -20.0, -20.0
        minmax1, minmax2 = get_value_excel('CHS-measurements.xlsx', 'Logger Hohes Holz', ['BC1_Ptemp [degC]', 'Tair50HMP [degC]'], ['Min', 'Max'])
        -> [-20.0, 50.0], [-20.0, 50.0]
        vars = get_value_excel('CHS-measurements.xlsx', 'Logger Hohes Holz', None, '')
        columns = get_value_excel('CHS-measurements.xlsx', 'Logger Hohes Holz', '', None)
        vars, columns = get_value_excel('CHS-measurements.xlsx', 'Logger Hohes Holz', None, None)


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

        Copyright 2015 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2015
        Modified, MC, Jul 2015 - column/variable can be iterable
                               - column/variable can None, returning names
    """
    # open Excel file
    wb = xlrd.open_workbook(excelfile)
    # check sheet exists
    if sheet not in wb.sheet_names(): raise ValueError('No sheet '+sheet+' in file '+excelfile)
    sh = wb.sheet_by_name(sheet)
    # Return varible and column names
    columns   = sh.row_values(0,0)
    variables = sh.col_values(columns.index('headerout (final)'), start_rowx=1, end_rowx=sh.nrows)
    if (variable is None) and (column is None):
        return [variables, columns]
    elif (variable is None):
        return variables
    elif (column is None):
        return columns
    # check columns exist
    hh = ['headerout (final)']
    if isinstance(column, (list, tuple, np.ndarray)):
        hh.extend(column)
    else:
        hh.append(column)
    for i in hh:
        if i not in columns:
            raise ValueError('Column '+i+' not found in sheet '+sheet+' of file '+excelfile)
    # check columns
    if isinstance(variable, (list, tuple, np.ndarray)):
        for ivar in variable:
            if ivar not in variables:
                raise ValueError(ivar+' not in column "headerout (final)" in sheet '+sheet+' of file '+excelfile)
    else:
        if variable not in variables:
            raise ValueError(variable+' not in column "headerout (final)" in sheet '+sheet+' of file '+excelfile)
    # extract cells
    if isinstance(column, (list, tuple, np.ndarray)) and isinstance(variable, (list, tuple, np.ndarray)):
        out = list()
        for ivar in variable:
            out.append([ sh.col_values(columns.index(icol),
                                       start_rowx=1, end_rowx=sh.nrows)[variables.index(ivar)] for icol in column ])
    else:
        if isinstance(column, (list, tuple, np.ndarray)):     # columns
            out = [ sh.col_values(columns.index(icol),
                                  start_rowx=1, end_rowx=sh.nrows)[variables.index(variable)] for icol in column ]
        elif isinstance(variable, (list, tuple, np.ndarray)): # variables
            out = [ sh.col_values(columns.index(column),
                                  start_rowx=1, end_rowx=sh.nrows)[variables.index(ivar)] for ivar in variable ]
        else:                                                 # scalars
            parameter = sh.col_values(columns.index(column), start_rowx=1, end_rowx=sh.nrows)
            out = parameter[variables.index(variable)]
    # relase Excel file
    del wb
    
    return out

# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # base_folder = '/Volumes/Gruppen/tereno/CHS-Data'
    # setup_chs = base_folder+'/CHS-measurements.xlsx'
    # sheet = 'Logger Hohes Holz'
    # v = 'Tair50HMP [degC]'
    # print(get_value_excel(setup_chs, sheet, v, 'Min'))
    # print(get_value_excel(setup_chs, sheet, v, ['Min','Max']))
    # v = ['Tair50HMP [degC]', 'BC1_Ptemp [degC]']
    # print(get_value_excel(setup_chs, sheet, v, 'Min'))
    # print(get_value_excel(setup_chs, sheet, v, ['Min','Max']))
    # # print(get_value_excel(setup_chs, 'Logger Hohes Holz', None, ''))
    # print(get_value_excel(setup_chs, 'Logger Hohes Holz', '', None))
    # # print(get_value_excel(setup_chs, 'Logger Hohes Holz', None, None))

