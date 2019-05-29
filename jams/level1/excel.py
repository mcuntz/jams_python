#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import xlrd

__all__ = ['get_header_excel', 'get_value_excel']

# --------------------------------------------------------------------

def get_header_excel(excelfile, sheet):
    """
        Get the header row of an Excel sheet.

        Same as get_value_excel(excelfile, sheet, '', None).


        Input
        -----
        excelfile   Filename of Excel file
        sheet       Name of Sheet in Excel file


        Output
        ------
        list with entries of first row of Excel sheet


        Examples
        --------
        --> see __init__.py for full example of workflow

        hh = get_header_excel('CHS-measurements.xlsx', 'Logger Hohes Holz')
        -> ['variable description', 'date start', 'date end', ... ]


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2016 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2016
    """
    # open Excel file
    wb = xlrd.open_workbook(excelfile)
    # check sheet exists
    if sheet not in wb.sheet_names(): raise ValueError('No sheet '+sheet+' in file '+excelfile)
    sh = wb.sheet_by_name(sheet)
    # Return variable and column names
    columns   = sh.row_values(0,0)

    # relase Excel file
    del wb

    return columns


# --------------------------------------------------------------------

def get_value_excel(excelfile, sheet, variable, column, all_rows=False):

    """
        Get a cell in a sheet of an Excel file
        where the column is identified by the name in the first row
        and the row is identified by the value in the column "headerout (final)".

        Both, variable and column can be iterables.

        If variable and/or column is None then the list of variables and columns is returned.


        Input
        -----
        excelfile   Filename of Excel file
        sheet       Name of Sheet in Excel file
        variable    Variable name in column "headerout (final)"
        column      Name in header line


        Optional Input
        --------------
        all_rows    option to return whole column, i.e. all rows (default: False).
                    Should be chosen if there are/may be duplicates in vars/the "headerout (final)" column;
                    the result of get_value_excel(file,sheet,vars,col)
                    if vars = get_value_excel(file,sheet,None,'') is NOT equal to
                    get_value_excel(file,sheet,None,col,all_rows=True)
                    as the former cannot deal properly with duplicates in vars.


        Output
        ------
        content of the desired cell
        if column or variable is iterable, output is list
        if column and variable is iterable, output is list (variable) of list (columns)
        if column and/or variable is not given, then a list of all variable and/or columns is returned
        if column is iterable and all_rows==True, output is list of list (all rows of all columns)


        Examples
        --------
        --> see __init__.py for full example of workflow

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
        This file is part of the JAMS Python package.

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

        Copyright 2015 Matthias Cuntz


        History
        -------
        Written,  MC, Mar 2015
        Modified, MC, Jul 2015 - column/variable can be iterable
                               - column/variable can None, returning names
                      Aug 2015 - added docu of input variables
                  BD, Nov 2016 - all_rows option to get the whole column (i.e. all rows)
                                 in cases where the varnames contain duplicates
    """
    # open Excel file
    wb = xlrd.open_workbook(excelfile)
    # check sheet exists
    if sheet not in wb.sheet_names(): raise ValueError('No sheet '+sheet+' in file '+excelfile)
    sh = wb.sheet_by_name(sheet)
    # Return variable and column names
    columns   = sh.row_values(0,0)
    variables = sh.col_values(columns.index('headerout (final)'), start_rowx=1, end_rowx=sh.nrows)
    if (variable is None) and (column is None):
        return [variables, columns]
    elif (variable is None) and (not all_rows):
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
    if isinstance(variable, (list, tuple, np.ndarray)) and (not all_rows):
        for ivar in variable:
            if ivar not in variables:
                raise ValueError(ivar+' not in column "headerout (final)" in sheet '+sheet+' of file '+excelfile)
    elif not all_rows:
        if variable not in variables:
            raise ValueError(variable+' not in column "headerout (final)" in sheet '+sheet+' of file '+excelfile)
    # extract cells
    if isinstance(column, (list, tuple, np.ndarray)) and isinstance(variable, (list, tuple, np.ndarray)) and (not all_rows):
        out = list()
        for ivar in variable:
            out.append([ sh.col_values(columns.index(icol),
                                       start_rowx=1, end_rowx=sh.nrows)[variables.index(ivar)] for icol in column ])
    elif all_rows: # gives all rows for the column if no variable is given and all_rows==True
        out=[ sh.col_values(columns.index(icol),
                        start_rowx=1, end_rowx=sh.nrows) for icol in column ]

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

    # base_folder = '/Volumes/Gruppen/CHS-Data'
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
    # print(get_value_excel(setup_chs, 'Logger Hohes Holz', None, ''))

    # print(get_header_excel(setup_chs, sheet))

    # v = 'BC1_Voltage [V]'
    # print(get_value_excel(setup_chs, sheet, v, ['Flag_0', 'Flag_1', 'Flag_2', 'Flag_3', 'Flag_4', 'Flag_5']))
