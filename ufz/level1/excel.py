#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import xlrd

__all__ = ['get_value_excel']

# --------------------------------------------------------------------

def get_value_excel(excelfile, sheet, variable, column):
    """
        Get the value in a column of a sheet in an excelfile
        given a variable name that is looked-up in the column
        "headerout (final)" in the same sheet.


        Definition
        ----------
        def get_value_excel(excelfile, sheet, variable, column):


        Input
        -----
        excelfile   filename of Excel file
        sheet       Name of sheet in Excel file
        variable    variable name as given in column "headerout (final)" in sheet
        column      column name of the parameter wanted

        
        Output
        ------
        content of the desired cell


        Examples
        --------
        mini = get_value_excel('CHS-measurements.xlsx', 'Forest Hohes Holz', 'BC1_Ptemp [degC]', 'Min')
        -> -20.0


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
    """
    # open Excel file
    wb = xlrd.open_workbook(excelfile)
    # check sheet exists
    if sheet not in wb.sheet_names(): raise ValueError('No sheet'+sheet+' in file '+excelfile)
    sh = wb.sheet_by_name(sheet)
    # check columns exist
    hh = ['headerout (final)', column]
    for i in hh:
        if i not in sh.row_values(0,0):
            raise ValueError('Column '+i+' not found in sheet '+sheet+' of file '+excelfile)
    # read the columns
    variables = sh.col_values(sh.row_values(0,0).index('headerout (final)'), start_rowx=1, end_rowx=sh.nrows)
    parameter = sh.col_values(sh.row_values(0,0).index(column), start_rowx=1, end_rowx=sh.nrows)
    if variable not in variables:
        raise ValueError(variable+' not in column "headerout (final)" in sheet '+sheet+' of file '+excelfile)
    # extract cell
    out = parameter[variables.index(variable)]
    # relase Excel file
    del wb
    
    return out

# --------------------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
