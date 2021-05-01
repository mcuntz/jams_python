#!/usr/bin/env python
from __future__ import absolute_import, division, print_function
'''
Pandas DataFrame GUI

A minimalistic GUI for analyzing Pandas DataFrames based on wxPython.


Usage
-----
import numpy as np
import pandas as pd
filename = 'Hesse_DB2_1997-2011.csv'
# filename = '1.csv'
# not needed for US format
parser = lambda date: pd.datetime.strptime(date, '%Y-%m-%d %H:%M')
names = ['Date', 'FCO2 H1 (umol/m2.s)', 'FCO2_cor H1 (umol/m2.s)',
         'FH2O H1 (mmol/m2.s)', 'H H1 (W/m2)',
         'Rg_Kipp H1 (W/m2)', 'alb H1 (-)', 'Rnet_CNR1 H1 (W/m2)',
         'NDVI H1 (-)', 'Patm (hPa)',
         'T_Vais H1 (degC)', 'WS_EC H1 (m/s)', 'Prec H1 (mm/30min)']

# pandas
df = pd.read_csv(filename, ';', parse_dates=[0], date_parser=parser,
                 index_col=0, header=0, usecols=names)

# jams
from jams import fread, fsread, ascii2eng
names = names[1:]
head = fread(filename, skip=1, cname=names, header=True)
dat, date = fsread(filename, skip=1, cname=names, sname=['Date'],
                   squeeze=True, strarr=True, strip=False)
df = pd.DataFrame(dat, pd.to_datetime(ascii2eng(date)), head)

# This must be before any other call to matplotlib
# because it use the wxAgg backend.
# This means, do not use --pylab with ipython.
df.replace(-9999., np.nan, inplace=True)
from jams import dfgui
dfgui.show(df)


Features
--------
- Tabular view of data frame
- Columns are sortable (by clicking column header)
- Columns can be enabled/disabled (left click on 'Columns' tab)
- Columns can be rearranged (right click drag on 'Columns' tab)
- Generic filtering: Write arbitrary Python expression to filter rows.
  *Warning:* Uses Python's `eval` -- use with care.
  Select column. Use filter like:
      _ > 0.
      _ > datetime.datetime(1998,1,1)
- Histogram plots
- Scatter plots


Dependencies
------------
Needs wxpython.


License
-------
This file is part of the JAMS Python package,
distributed under the MIT License.

Copyright (c) 2019 Matthias Cuntz - mc (at) macu (dot) de

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
Written, Fabian Keller (fabian.keller@blue-yonder.com), 2014,
         http://github.com/bluenote10/PandasDataFrameGUI
Modified, Matthias Cuntz, Jan 2019 - exclude NaN in histograms
                                   - add index as 'Date' column automatically
                                   - plot controls
          Matthias Cuntz, Feb 2019 - added 2nd line/scatter in scatter plots
                                     -> Bug: time axis does not work with first
                                             line but with second line.
                                   - label x-axis on histograms
          Matthias Cuntz, Jun 2019 - make application
                                   - Removed bug that x-axis in scatter plots
                                     is only correct if line2 is chosen.
          Matthias Cuntz, Jul 2019 - added example to filter panel
          Matthias Cuntz, Aug 2019 - added checkbox for linking the two y-axes
                                     in scatter plots
          Matthias Cuntz, Sep 2020 - added checkbox to invert second axis
          Matthias Cuntz, Sep 2020 - changed wx.ListItemAttr() to wx.ItemAttr()
          Matthias Cuntz, Apr 2021 - undef on command line
                                   - flake8 compatible
          Matthias Cuntz, May 2021 - sort columns with command line option -c
'''
# --------------------------------------------------------------------
# import
#

# pip install wxpython
import numpy as np
import pandas as pd
# unused import required to allow 'eval' of date filters
import datetime
from bisect import bisect

try:
    import wx
except ImportError:
    import sys
    sys.path += ["/usr/lib/python2.7/dist-packages/wx-2.8-gtk2-unicode",
                 "/usr/lib/python2.7/dist-packages"]
    import wx

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

# try to get nicer plotting styles
try:
    import seaborn
    seaborn.set()
except ImportError:
    try:
        from matplotlib import pyplot as plt
        plt.style.use('ggplot')
    except AttributeError:
        pass


# --------------------------------------------------------------------
# Basic List Class
#
class ListCtrlDataFrame(wx.ListCtrl):

    # TODO: we could do something more sophisticated to come
    #       up with a reasonable column width...
    DEFAULT_COLUMN_WIDTH = 100
    TMP_SELECTION_COLUMN = 'tmp_selection_column'

    def __init__(self, parent, df, status_bar_callback):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT |
                             wx.LC_VIRTUAL | wx.LC_HRULES | wx.LC_VRULES |
                             wx.LB_MULTIPLE)
        self.status_bar_callback = status_bar_callback

        self.df_orig = df
        self.original_columns = self.df_orig.columns[:]
        self.current_columns = self.df_orig.columns[:]

        self.sort_by_column = None

        self._reset_mask()

        # prepare attribute for alternating colors of rows
        # self.attr_light_blue = wx.ListItemAttr()
        self.attr_light_blue = wx.ItemAttr()
        self.attr_light_blue.SetBackgroundColour("#D6EBFF")

        self.Bind(wx.EVT_LIST_COL_CLICK, self._on_col_click)
        self.Bind(wx.EVT_RIGHT_DOWN, self._on_right_click)

        self.df = pd.DataFrame({})  # init empty to force initial update
        self._update_rows()
        self._update_columns(self.original_columns)

    def _reset_mask(self):
        # self.mask = [True] * self.df_orig.shape[0]
        self.mask = pd.Series([True] * self.df_orig.shape[0],
                              index=self.df_orig.index)

    def _update_columns(self, columns):
        self.ClearAll()
        for i, col in enumerate(columns):
            self.InsertColumn(i, col)
            self.SetColumnWidth(i, self.DEFAULT_COLUMN_WIDTH)
        # Note that we have to reset the count as well because ClearAll()
        # not only deletes columns but also the count...
        self.SetItemCount(len(self.df))

    def set_columns(self, columns_to_use):
        """
        External interface to set the column projections.
        """
        self.current_columns = columns_to_use
        self._update_rows()
        self._update_columns(columns_to_use)

    def _update_rows(self):
        old_len = len(self.df)
        self.df = self.df_orig.loc[self.mask.values, self.current_columns]
        new_len = len(self.df)
        if old_len != new_len:
            self.SetItemCount(new_len)
            self.status_bar_callback(0, "Number of rows: {}".format(new_len))

    def apply_filter(self, conditions):
        """
        External interface to set a filter.
        """
        old_mask = self.mask.copy()

        if len(conditions) == 0:
            self._reset_mask()

        else:
            self._reset_mask()  # set all to True for destructive conjunction

            no_error = True
            for column, condition in conditions:
                if condition.strip() == '':
                    continue
                condition = condition.replace(
                    "_", "self.df_orig['{}']".format(column))
                print("Evaluating condition:", condition)
                try:
                    tmp_mask = eval(condition)
                    if (isinstance(tmp_mask, pd.Series) and
                        (tmp_mask.dtype == np.bool)):
                        self.mask &= tmp_mask
                except Exception as e:
                    print("Failed with:", e)
                    no_error = False
                    self.status_bar_callback(
                        1,
                        "Evaluating '{}' failed with: {}".format(condition, e)
                    )

            if no_error:
                self.status_bar_callback(1, "")

        has_changed = any(old_mask != self.mask)
        if has_changed:
            self._update_rows()

        return len(self.df), has_changed

    def get_selected_items(self):
        """
        Gets the selected items for the list control.
        Selection is returned as a list of selected indices,
        low to high.
        """
        selection = []
        current = -1    # start at -1 to get the first selected item
        while True:
            next = self.GetNextItem(current, wx.LIST_NEXT_ALL,
                                    wx.LIST_STATE_SELECTED)
            if next == -1:
                return selection
            else:
                selection.append(next)
                current = next

    def get_filtered_df(self):
        return self.df_orig.loc[self.mask, :]

    def _on_col_click(self, event):
        """
        Sort data frame by selected column.
        """
        # get currently selected items
        selected = self.get_selected_items()

        # append a temporary column to store the currently selected items
        self.df[self.TMP_SELECTION_COLUMN] = False
        self.df.iloc[selected, -1] = True

        # get column name to use for sorting
        col = event.GetColumn()

        # determine if ascending or descending
        if self.sort_by_column is None or self.sort_by_column[0] != col:
            ascending = True
        else:
            ascending = not self.sort_by_column[1]

        # store sort column and sort direction
        self.sort_by_column = (col, ascending)

        try:
            # pandas 0.17
            self.df.sort_values(self.df.columns[col], inplace=True,
                                ascending=ascending)
        except AttributeError:
            # pandas 0.16 compatibility
            self.df.sort(self.df.columns[col], inplace=True,
                         ascending=ascending)

        # deselect all previously selected
        for i in selected:
            self.Select(i, on=False)

        # determine indices of selection after sorting
        selected_bool = self.df.iloc[:, -1] == True
        selected = self.df.reset_index().index[selected_bool]

        # select corresponding rows
        for i in selected:
            self.Select(i, on=True)

        # delete temporary column
        del self.df[self.TMP_SELECTION_COLUMN]

    def _on_right_click(self, event):
        """
        Copies a cell into clipboard on right click. Unfortunately,
        determining the clicked column is not straightforward. This
        appraoch is inspired by the TextEditMixin in:
        /usr/lib/python2.7/dist-packages/wx-2.8-gtk2-unicode/wx/lib/mixins/listctrl.py
        More references:
        http://wxpython-users.1045709.n5.nabble.com/Getting-row-col-of-selected-cell-in-ListCtrl-td2360831.html
        https://groups.google.com/forum/#!topic/wxpython-users/7BNl9TA5Y5U
        https://groups.google.com/forum/#!topic/wxpython-users/wyayJIARG8c
        """
        if self.HitTest(event.GetPosition()) != wx.NOT_FOUND:
            x, y = event.GetPosition()
            row, flags = self.HitTest((x, y))

            col_locs = [0]
            loc = 0
            for n in range(self.GetColumnCount()):
                loc = loc + self.GetColumnWidth(n)
                col_locs.append(loc)

            scroll_pos = self.GetScrollPos(wx.HORIZONTAL)
            # this is crucial step to get the scroll pixel units
            unit_x, unit_y = self.GetMainWindow().GetScrollPixelsPerUnit()

            col = bisect(col_locs, x + scroll_pos * unit_x) - 1

            value = self.df.iloc[row, col]
            # print(row, col, scroll_pos, value)

            clipdata = wx.TextDataObject()
            clipdata.SetText(str(value))
            wx.TheClipboard.Open()
            wx.TheClipboard.SetData(clipdata)
            wx.TheClipboard.Close()

    def OnGetItemText(self, item, col):
        """
        Implements the item getter for a "virtual" ListCtrl.
        """
        value = self.df.iloc[item, col]
        # print("retrieving %d %d %s" % (item, col, value))
        return str(value)

    def OnGetItemAttr(self, item):
        """
        Implements the attribute getter for a "virtual" ListCtrl.
        """
        if item % 2 == 0:
            return self.attr_light_blue
        else:
            return None


# --------------------------------------------------------------------
# Date table panel
#
class DataframePanel(wx.Panel):
    """
    Panel providing the main data frame table view.
    """
    def __init__(self, parent, df, status_bar_callback):
        wx.Panel.__init__(self, parent)

        self.df_list_ctrl = ListCtrlDataFrame(self, df, status_bar_callback)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.df_list_ctrl, 1, wx.ALL | wx.EXPAND | wx.GROW, 5)
        self.SetSizer(sizer)
        self.Show()


class ListBoxDraggable(wx.ListBox):
    """
    Helper class to provide ListBox with extended behavior.
    """
    def __init__(self, parent, size, data, *args, **kwargs):

        wx.ListBox.__init__(self, parent, size, **kwargs)

        self.data = data

        self.InsertItems(data, 0)

        self.Bind(wx.EVT_LISTBOX, self.on_selection_changed)

        self.Bind(wx.EVT_LEFT_DOWN, self.on_left_down)

        self.Bind(wx.EVT_RIGHT_DOWN, self.on_right_down)
        self.Bind(wx.EVT_RIGHT_UP, self.on_right_up)
        self.Bind(wx.EVT_MOTION, self.on_move)

        self.index_iter = range(len(self.data))

        self.selected_items = [True] * len(self.data)
        self.index_mapping = list(range(len(self.data)))

        self.drag_start_index = None

        self.update_selection()
        self.SetFocus()

    def on_left_down(self, event):
        if self.HitTest(event.GetPosition()) != wx.NOT_FOUND:
            index = self.HitTest(event.GetPosition())
            self.selected_items[index] = not self.selected_items[index]
            # doesn't really work to update selection direclty (focus issues)
            # instead we wait for the EVT_LISTBOX event and fix the selection
            # there...
            # self.update_selection()
            # TODO: we could probably use wx.CallAfter
        event.Skip()

    def update_selection(self):
        # self.SetFocus()
        # print(self.selected_items)
        for i in self.index_iter:
            if self.IsSelected(i) and not self.selected_items[i]:
                # print("Deselecting", i)
                self.Deselect(i)
            elif not self.IsSelected(i) and self.selected_items[i]:
                # print("Selecting", i)
                self.Select(i)

    def on_selection_changed(self, evt):
        self.update_selection()
        evt.Skip()

    def on_right_down(self, event):
        if self.HitTest(event.GetPosition()) != wx.NOT_FOUND:
            index = self.HitTest(event.GetPosition())
            self.drag_start_index = index

    def on_right_up(self, event):
        self.drag_start_index = None
        event.Skip()

    def on_move(self, event):
        if self.drag_start_index is not None:
            if self.HitTest(event.GetPosition()) != wx.NOT_FOUND:
                index = self.HitTest(event.GetPosition())
                if self.drag_start_index != index:
                    self.swap(self.drag_start_index, index)
                    self.drag_start_index = index

    def swap(self, i, j):
        self.index_mapping[i], self.index_mapping[j] = self.index_mapping[j], self.index_mapping[i]
        self.SetString(i, self.data[self.index_mapping[i]])
        self.SetString(j, self.data[self.index_mapping[j]])
        self.selected_items[i], self.selected_items[j] = self.selected_items[j], self.selected_items[i]
        # self.update_selection()
        # print("Updated mapping:", self.index_mapping)
        new_event = wx.PyCommandEvent(wx.EVT_LISTBOX.typeId, self.GetId())
        self.GetEventHandler().ProcessEvent(new_event)

    def get_selected_data(self):
        selected = []
        for i, col in enumerate(self.data):
            if self.IsSelected(i):
                index = self.index_mapping[i]
                value = self.data[index]
                selected.append(value)
        # print("Selected data:", selected)
        return selected


# --------------------------------------------------------------------
# Column selection panel
#
class ColumnSelectionPanel(wx.Panel):
    """
    Panel for selecting and re-arranging columns.
    """
    def __init__(self, parent, columns, df_list_ctrl):
        wx.Panel.__init__(self, parent)

        self.columns = columns
        self.df_list_ctrl = df_list_ctrl

        self.list_box = ListBoxDraggable(self, -1, columns,
                                         style=wx.LB_EXTENDED)
        self.Bind(wx.EVT_LISTBOX, self.update_selected_columns)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.list_box, 1, wx.ALL | wx.EXPAND | wx.GROW, 5)
        self.SetSizer(sizer)
        self.list_box.SetFocus()

    def update_selected_columns(self, evt):
        selected = self.list_box.get_selected_data()
        self.df_list_ctrl.set_columns(selected)


# --------------------------------------------------------------------
# Filter panel
#
class FilterPanel(wx.Panel):
    """
    Panel for defining filter expressions.
    """
    def __init__(self, parent, columns, df_list_ctrl, change_callback):
        wx.Panel.__init__(self, parent)

        columns_with_neutral_selection = [''] + list(columns)
        self.columns = columns
        self.df_list_ctrl = df_list_ctrl
        self.change_callback = change_callback

        self.num_filters = 10

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        self.example = wx.StaticText(self, wx.ID_ANY,
                                     label="Select column. Use filter like:   _ > 0.   or   _ > datetime.datetime(1998,1,1)",
                                     size=wx.Size(560, 20))
        row_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # row_sizer.Add(self.example, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        row_sizer.Add(self.example, 0, wx.ALL, 1)
        self.main_sizer.Add(row_sizer, 0, wx.EXPAND)

        self.combo_boxes = []
        self.text_controls = []

        for i in range(self.num_filters):
            combo_box = wx.ComboBox(self,
                                    choices=columns_with_neutral_selection,
                                    style=wx.CB_READONLY)
            text_ctrl = wx.TextCtrl(self, wx.ID_ANY, '')

            self.Bind(wx.EVT_COMBOBOX, self.on_combo_box_select)
            self.Bind(wx.EVT_TEXT, self.on_text_change)

            row_sizer = wx.BoxSizer(wx.HORIZONTAL)
            row_sizer.Add(combo_box, 0, wx.ALL, 5)
            # row_sizer.Add(text_ctrl, 1, wx.ALL | wx.EXPAND |
            #               wx.ALIGN_RIGHT, 5)
            row_sizer.Add(text_ctrl, 1, wx.ALL | wx.EXPAND, 5)

            self.combo_boxes.append(combo_box)
            self.text_controls.append(text_ctrl)
            self.main_sizer.Add(row_sizer, 0, wx.EXPAND)

        self.SetSizer(self.main_sizer)

    def on_combo_box_select(self, event):
        self.update_conditions()

    def on_text_change(self, event):
        self.update_conditions()

    def update_conditions(self):
        # print("Updating conditions")
        conditions = []
        for i in range(self.num_filters):
            column_index = self.combo_boxes[i].GetSelection()
            condition = self.text_controls[i].GetValue()
            if column_index != wx.NOT_FOUND and column_index != 0:
                # since we have added a dummy column for "deselect",
                # we have to subtract one
                column = self.columns[column_index - 1]
                conditions += [(column, condition)]
        num_matching, has_changed = self.df_list_ctrl.apply_filter(conditions)
        if has_changed:
            self.change_callback()
        # print("Num matching:", num_matching)


# --------------------------------------------------------------------
# Histogram plot panel
#
class HistogramPlot(wx.Panel):
    """
    Panel providing a histogram plot.
    """
    def __init__(self, parent, columns, df_list_ctrl):
        wx.Panel.__init__(self, parent)

        columns_with_neutral_selection = [''] + list(columns)
        self.columns = columns
        self.df_list_ctrl = df_list_ctrl

        self.figure = Figure(facecolor="white", figsize=(1, 1))
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)

        chart_toolbar = NavigationToolbar2Wx(self.canvas)

        self.combo_box1 = wx.ComboBox(self,
                                      choices=columns_with_neutral_selection,
                                      style=wx.CB_READONLY)

        self.Bind(wx.EVT_COMBOBOX, self.on_combo_box_select)

        row_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # row_sizer.Add(self.combo_box1, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        row_sizer.Add(self.combo_box1, 0, wx.ALL, 5)
        row_sizer.Add(chart_toolbar, 0, wx.ALL, 5)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, flag=wx.EXPAND, border=5)
        sizer.Add(row_sizer)
        self.SetSizer(sizer)

    def on_combo_box_select(self, event):
        self.redraw()

    def redraw(self):
        column_index1 = self.combo_box1.GetSelection()
        if column_index1 != wx.NOT_FOUND and column_index1 != 0:
            # subtract one to remove the neutral selection index
            column_index1 -= 1

            df = self.df_list_ctrl.get_filtered_df()

            if len(df) > 0:
                self.axes.clear()

                column = df.iloc[:, column_index1]
                is_string_col = ((column.dtype == np.object) and
                                 isinstance(column.values[0], str))
                if is_string_col:
                    value_counts = column.value_counts().sort_index()
                    value_counts.plot(kind='bar', ax=self.axes)
                else:
                    self.axes.hist(
                        np.ma.array(column.values,
                                    mask=~np.isfinite(column.values)),
                        bins=100)
                    self.axes.xaxis.set_label_text(self.columns[column_index1])

                self.canvas.draw()


# --------------------------------------------------------------------
# Scatter plot panel
#
class ScatterPlot(wx.Panel):
    """
    Panel providing a scatter plot.
    """
    def __init__(self, parent, columns, df_list_ctrl):
        wx.Panel.__init__(self, parent)

        columns_with_neutral_selection = [''] + list(columns)
        self.columns = columns
        self.df_list_ctrl = df_list_ctrl

        self.figure = Figure(facecolor="white", figsize=(1, 1))
        self.axes  = self.figure.add_subplot(111)
        self.axes2 = self.axes.twinx()
        self.canvas = FigureCanvas(self, -1, self.figure)

        chart_toolbar = NavigationToolbar2Wx(self.canvas)

        # first line
        self.combo_box1 = wx.ComboBox(self,
                                      choices=columns_with_neutral_selection,
                                      style=wx.CB_READONLY)
        self.combo_box2 = wx.ComboBox(self,
                                      choices=columns_with_neutral_selection,
                                      style=wx.CB_READONLY)

        self.linestyle_label = wx.StaticText(self, wx.ID_ANY,
                                             label="linestyle",
                                             size=wx.Size(60, 20))
        self.linestyle = wx.TextCtrl(self, wx.ID_ANY, value="-",
                                     size=wx.Size(40, 20))
        self.linewidth_label = wx.StaticText(self, wx.ID_ANY,
                                             label="linewidth",
                                             size=wx.Size(60, 20))
        self.linewidth = wx.TextCtrl(self, wx.ID_ANY, value="1",
                                     size=wx.Size(25, 20))
        self.linecolor_label = wx.StaticText(self, wx.ID_ANY,
                                             label="linecolor",
                                             size=wx.Size(60, 20))
        self.linecolor = wx.TextCtrl(self, wx.ID_ANY, value="b",
                                     size=wx.Size(40, 20))
        self.marker_label = wx.StaticText(self, wx.ID_ANY, label="marker",
                                          size=wx.Size(50, 20))
        self.marker = wx.TextCtrl(self, wx.ID_ANY, value="None",
                                  size=wx.Size(40, 20))
        self.markersize_label = wx.StaticText(self, wx.ID_ANY,
                                              label="markersize",
                                              size=wx.Size(70, 20))
        self.markersize = wx.TextCtrl(self, wx.ID_ANY, value="1",
                                      size=wx.Size(25, 20))
        self.markerfacecolor_label = wx.StaticText(self, wx.ID_ANY,
                                                   label="markerfacecolor",
                                                   size=wx.Size(100, 20))
        self.markerfacecolor = wx.TextCtrl(self, wx.ID_ANY, value="b",
                                           size=wx.Size(40, 20))
        self.markeredgecolor_label = wx.StaticText(self, wx.ID_ANY,
                                                   label="markeredgecolor",
                                                   size=wx.Size(110, 20))
        self.markeredgecolor = wx.TextCtrl(self, wx.ID_ANY, value="b",
                                           size=wx.Size(40, 20))
        self.markeredgewidth_label = wx.StaticText(self, wx.ID_ANY,
                                                   label="markeredgewidth",
                                                   size=wx.Size(110, 20))
        self.markeredgewidth = wx.TextCtrl(self, wx.ID_ANY, value="1",
                                           size=wx.Size(25, 20))

        row_sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        # row_sizer1.Add(self.combo_box1, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer1.Add(self.combo_box2, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        row_sizer1.Add(self.combo_box1, 0, wx.ALL, 5)
        row_sizer1.Add(self.combo_box2, 0, wx.ALL, 5)
        row_sizer1.Add(chart_toolbar, 0, wx.ALL, 5)

        row_sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        # row_sizer2.Add(self.linestyle_label, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer2.Add(self.linestyle, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer2.Add(self.linewidth_label, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer2.Add(self.linewidth, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer2.Add(self.linecolor_label, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer2.Add(self.linecolor, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        row_sizer2.Add(self.linestyle_label, 0, wx.ALL, 1)
        row_sizer2.Add(self.linestyle, 0, wx.ALL, 5)
        row_sizer2.Add(self.linewidth_label, 0, wx.ALL, 1)
        row_sizer2.Add(self.linewidth, 0, wx.ALL, 5)
        row_sizer2.Add(self.linecolor_label, 0, wx.ALL, 1)
        row_sizer2.Add(self.linecolor, 0, wx.ALL, 5)

        row_sizer3 = wx.BoxSizer(wx.HORIZONTAL)
        # row_sizer3.Add(self.marker_label, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer3.Add(self.marker, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer3.Add(self.markersize_label, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer3.Add(self.markersize, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer3.Add(self.markerfacecolor_label, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer3.Add(self.markerfacecolor, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer3.Add(self.markeredgecolor_label, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer3.Add(self.markeredgecolor, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer3.Add(self.markeredgewidth_label, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer3.Add(self.markeredgewidth, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 5)
        row_sizer3.Add(self.marker_label, 0, wx.ALL, 1)
        row_sizer3.Add(self.marker, 0, wx.ALL, 5)
        row_sizer3.Add(self.markersize_label, 0, wx.ALL, 1)
        row_sizer3.Add(self.markersize, 0, wx.ALL, 5)
        row_sizer3.Add(self.markerfacecolor_label, 0, wx.ALL, 1)
        row_sizer3.Add(self.markerfacecolor, 0, wx.ALL, 5)
        row_sizer3.Add(self.markeredgecolor_label, 0, wx.ALL, 1)
        row_sizer3.Add(self.markeredgecolor, 0, wx.ALL, 5)
        row_sizer3.Add(self.markeredgewidth_label, 0, wx.ALL, 1)
        row_sizer3.Add(self.markeredgewidth, 0, wx.ALL, 5)

        # second line
        self.combo_box3 = wx.ComboBox(self,
                                      choices=columns_with_neutral_selection,
                                      style=wx.CB_READONLY)
        self.checkbox1 = wx.CheckBox(self, label="Same y-axis")
        self.checkbox2 = wx.CheckBox(self, label="Invert y-axis")

        self.linestyle_label2 = wx.StaticText(self, wx.ID_ANY,
                                              label="linestyle",
                                              size=wx.Size(60, 20))
        self.linestyle2 = wx.TextCtrl(self, wx.ID_ANY, value="--",
                                      size=wx.Size(40, 20))
        self.linewidth_label2 = wx.StaticText(self, wx.ID_ANY,
                                              label="linewidth",
                                              size=wx.Size(60, 20))
        self.linewidth2 = wx.TextCtrl(self, wx.ID_ANY, value="1",
                                      size=wx.Size(25, 20))
        self.linecolor_label2 = wx.StaticText(self, wx.ID_ANY,
                                              label="linecolor",
                                              size=wx.Size(60, 20))
        self.linecolor2 = wx.TextCtrl(self, wx.ID_ANY, value="r",
                                      size=wx.Size(40, 20))
        self.marker_label2 = wx.StaticText(self, wx.ID_ANY, label="marker",
                                           size=wx.Size(50, 20))
        self.marker2 = wx.TextCtrl(self, wx.ID_ANY, value="None",
                                   size=wx.Size(40, 20))
        self.markersize_label2 = wx.StaticText(self, wx.ID_ANY,
                                               label="markersize",
                                               size=wx.Size(70, 20))
        self.markersize2 = wx.TextCtrl(self, wx.ID_ANY, value="1",
                                       size=wx.Size(25, 20))
        self.markerfacecolor_label2 = wx.StaticText(self, wx.ID_ANY,
                                                    label="markerfacecolor",
                                                    size=wx.Size(100, 20))
        self.markerfacecolor2 = wx.TextCtrl(self, wx.ID_ANY, value="r",
                                            size=wx.Size(40, 20))
        self.markeredgecolor_label2 = wx.StaticText(self, wx.ID_ANY,
                                                    label="markeredgecolor",
                                                    size=wx.Size(110, 20))
        self.markeredgecolor2 = wx.TextCtrl(self, wx.ID_ANY, value="r",
                                            size=wx.Size(40, 20))
        self.markeredgewidth_label2 = wx.StaticText(self, wx.ID_ANY,
                                                    label="markeredgewidth",
                                                    size=wx.Size(110, 20))
        self.markeredgewidth2 = wx.TextCtrl(self, wx.ID_ANY, value="1",
                                            size=wx.Size(25, 20))

        row_sizer4 = wx.BoxSizer(wx.HORIZONTAL)
        # row_sizer4.Add(self.combo_box3, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer4.Add(self.checkbox1, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        row_sizer4.Add(self.combo_box3, 0, wx.ALL, 5)
        row_sizer4.Add(self.checkbox1, 0, wx.ALL, 5)
        row_sizer4.Add(self.checkbox2, 0, wx.ALL, 5)

        row_sizer5 = wx.BoxSizer(wx.HORIZONTAL)
        # row_sizer5.Add(self.linestyle_label2, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer5.Add(self.linestyle2, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer5.Add(self.linewidth_label2, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer5.Add(self.linewidth2, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer5.Add(self.linecolor_label2, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer5.Add(self.linecolor2, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        row_sizer5.Add(self.linestyle_label2, 0, wx.ALL, 1)
        row_sizer5.Add(self.linestyle2, 0, wx.ALL, 5)
        row_sizer5.Add(self.linewidth_label2, 0, wx.ALL, 1)
        row_sizer5.Add(self.linewidth2, 0, wx.ALL, 5)
        row_sizer5.Add(self.linecolor_label2, 0, wx.ALL, 1)
        row_sizer5.Add(self.linecolor2, 0, wx.ALL, 5)

        row_sizer6 = wx.BoxSizer(wx.HORIZONTAL)
        # row_sizer6.Add(self.marker_label2, 0, wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer6.Add(self.marker2, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer6.Add(self.markersize_label2, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer6.Add(self.markersize2, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer6.Add(self.markerfacecolor_label2, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer6.Add(self.markerfacecolor2, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer6.Add(self.markeredgecolor_label2, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer6.Add(self.markeredgecolor2, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # row_sizer6.Add(self.markeredgewidth_label2, 0,
        #                wx.ALL | wx.ALIGN_CENTER, 1)
        # row_sizer6.Add(self.markeredgewidth2, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        row_sizer6.Add(self.marker_label2, 0, wx.ALL, 1)
        row_sizer6.Add(self.marker2, 0, wx.ALL, 5)
        row_sizer6.Add(self.markersize_label2, 0, wx.ALL, 1)
        row_sizer6.Add(self.markersize2, 0, wx.ALL, 5)
        row_sizer6.Add(self.markerfacecolor_label2, 0, wx.ALL, 1)
        row_sizer6.Add(self.markerfacecolor2, 0, wx.ALL, 5)
        row_sizer6.Add(self.markeredgecolor_label2, 0, wx.ALL, 1)
        row_sizer6.Add(self.markeredgecolor2, 0, wx.ALL, 5)
        row_sizer6.Add(self.markeredgewidth_label2, 0, wx.ALL, 1)
        row_sizer6.Add(self.markeredgewidth2, 0, wx.ALL, 5)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, flag=wx.EXPAND, border=5)
        sizer.Add(row_sizer1)
        sizer.Add(row_sizer2)
        sizer.Add(row_sizer3)
        sizer.Add(row_sizer4)
        sizer.Add(row_sizer5)
        sizer.Add(row_sizer6)
        self.SetSizer(sizer)

        self.Bind(wx.EVT_COMBOBOX, self.on_combo_box_select)
        self.Bind(wx.EVT_TEXT, self.on_text_change)
        self.Bind(wx.EVT_CHECKBOX, self.on_checkbox)

    def on_combo_box_select(self, event):
        self.redraw()

    def on_text_change(self, event):
        self.redraw()

    def on_checkbox(self, event):
        self.redraw()

    def redraw(self, event=None):
        column_index1 = self.combo_box1.GetSelection()
        column_index2 = self.combo_box2.GetSelection()
        linestyle       = self.linestyle.GetLineText(0)
        linewidth       = float(self.linewidth.GetLineText(0))
        linecolor       = self.linecolor.GetLineText(0)
        marker          = self.marker.GetLineText(0)
        markersize      = float(self.markersize.GetLineText(0))
        markerfacecolor = self.markerfacecolor.GetLineText(0)
        markeredgecolor = self.markeredgecolor.GetLineText(0)
        markeredgewidth = float(self.markeredgewidth.GetLineText(0))
        # line 2
        column_index3 = self.combo_box3.GetSelection()
        linestyle2       = self.linestyle2.GetLineText(0)
        linewidth2       = float(self.linewidth2.GetLineText(0))
        linecolor2       = self.linecolor2.GetLineText(0)
        marker2          = self.marker2.GetLineText(0)
        markersize2      = float(self.markersize2.GetLineText(0))
        markerfacecolor2 = self.markerfacecolor2.GetLineText(0)
        markeredgecolor2 = self.markeredgecolor2.GetLineText(0)
        markeredgewidth2 = float(self.markeredgewidth2.GetLineText(0))
        sameyaxes   = self.checkbox1.GetValue()
        invertyaxis = self.checkbox2.GetValue()
        if ((column_index1 != wx.NOT_FOUND and column_index2 != wx.NOT_FOUND
             and column_index1 != 0 and column_index2 != 0) or
            (column_index1 != wx.NOT_FOUND and column_index3 != wx.NOT_FOUND
             and column_index1 != 0 and column_index3 != 0)):
            df = self.df_list_ctrl.get_filtered_df()
            if len(df) > 0:
                # subtract one to remove the neutral selection index
                column_index1 -= 1
                x  = df.iloc[:, column_index1].values
                column_index2 -= 1
                y  = df.iloc[:, column_index2].values
                column_index3 -= 1
                y2 = df.iloc[:, column_index3].values
                # It looks like using pandas dataframe.plot causes something
                # weird to crash in wx internally. Therefore we use plain
                # axes.plot functionality.
                # column_name1 = self.columns[column_index1]
                # column_name2 = self.columns[column_index2]
                # df.plot(kind='scatter', x=column_name1, y=column_name2)
                # Clear both axes first, otherwise x-axis only shows if line2
                # is chosen
                self.axes.clear()
                self.axes2.clear()
                ylim1 = [None, None]
                ylim2 = [None, None]
                if (column_index1 >= 0 and column_index2 >= 0):
                    self.axes.plot(x, y, linestyle=linestyle,
                                   linewidth=linewidth, color=linecolor,
                                   marker=marker, markersize=markersize,
                                   markerfacecolor=markerfacecolor,
                                   markeredgecolor=markeredgecolor,
                                   markeredgewidth=markeredgewidth)
                    self.axes.spines['left'].set_color(linecolor)
                    self.axes.tick_params(axis='y', colors=linecolor)
                    self.axes.xaxis.set_label_text(self.columns[column_index1])
                    self.axes.yaxis.label.set_color(linecolor)
                    self.axes.yaxis.set_label_text(self.columns[column_index2])
                    ylim1 = self.axes.get_ylim()
                if (column_index1 >= 0 and column_index3 >= 0):
                    self.axes2.plot(x, y2, linestyle=linestyle2,
                                    linewidth=linewidth2, color=linecolor2,
                                    marker=marker2, markersize=markersize2,
                                    markerfacecolor=markerfacecolor2,
                                    markeredgecolor=markeredgecolor2,
                                    markeredgewidth=markeredgewidth2)
                    self.axes2.spines['right'].set_color(linecolor2)
                    self.axes2.tick_params(axis='y', colors=linecolor2)
                    self.axes2.xaxis.set_label_text(
                        self.columns[column_index1])
                    self.axes2.yaxis.label.set_color(linecolor2)
                    self.axes2.yaxis.set_label_text(
                        self.columns[column_index3])
                    ylim2 = self.axes2.get_ylim()
                if sameyaxes:
                    if (ylim1[0] is not None) and (ylim2[0] is not None):
                        ymin = min(ylim1[0], ylim2[0])
                    else:
                        if (ylim1[0] is not None):
                            ymin = ylim1[0]
                        else:
                            ymin = ylim2[0]
                    if (ylim1[1] is not None) and (ylim2[1] is not None):
                        ymax = max(ylim1[1], ylim2[1])
                    else:
                        if (ylim1[1] is not None):
                            ymax = ylim1[1]
                        else:
                            ymax = ylim2[1]
                    if (ymin is not None) and (ymax is not None):
                        self.axes.set_ylim([ymin, ymax])
                        self.axes2.set_ylim([ymin, ymax])
                elif invertyaxis:
                    if (ylim2[0] is not None):
                        ylim2 = ylim2[::-1]
                        self.axes2.set_ylim(ylim2)
                self.canvas.draw()


# --------------------------------------------------------------------
# Main window
#
class MainFrame(wx.Frame):
    """
    The main GUI window.
    """
    def __init__(self, df):
        wx.Frame.__init__(self, None, -1, "Pandas DataFrame GUI")

        # Here we create a panel and a notebook on the panel
        p = wx.Panel(self)
        nb = wx.Notebook(p)
        self.nb = nb

        # add index as Date column
        df.insert(0, 'Date', pd.Series(pd.to_datetime(df.index),
                                       index=df.index))

        columns = df.columns[:]

        self.CreateStatusBar(2, style=0)
        self.SetStatusWidths([200, -1])

        # create the page windows as children of the notebook
        self.page1 = DataframePanel(nb, df, self.status_bar_callback)
        self.page2 = ColumnSelectionPanel(nb, columns, self.page1.df_list_ctrl)
        self.page3 = FilterPanel(nb, columns, self.page1.df_list_ctrl,
                                 self.selection_change_callback)
        self.page4 = HistogramPlot(nb, columns, self.page1.df_list_ctrl)
        self.page5 = ScatterPlot(nb, columns, self.page1.df_list_ctrl)

        # add the pages to the notebook with the label to show on the tab
        nb.AddPage(self.page1, "Data Frame")
        nb.AddPage(self.page2, "Columns")
        nb.AddPage(self.page3, "Filters")
        nb.AddPage(self.page4, "Histogram")
        nb.AddPage(self.page5, "Scatter Plot")

        nb.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.on_tab_change)

        # finally, put the notebook in a sizer for the panel to manage
        # the layout
        sizer = wx.BoxSizer()
        sizer.Add(nb, 1, wx.EXPAND)
        p.SetSizer(sizer)

        self.SetSize((1000, 800))
        self.Center()

    def on_tab_change(self, event):
        self.page2.list_box.SetFocus()
        page_to_select = event.GetSelection()
        wx.CallAfter(self.fix_focus, page_to_select)
        event.Skip(True)

    def fix_focus(self, page_to_select):
        page = self.nb.GetPage(page_to_select)
        page.SetFocus()
        if isinstance(page, DataframePanel):
            self.page1.df_list_ctrl.SetFocus()
        elif isinstance(page, ColumnSelectionPanel):
            self.page2.list_box.SetFocus()

    def status_bar_callback(self, i, new_text):
        self.SetStatusText(new_text, i)

    def selection_change_callback(self):
        self.page4.redraw()
        self.page5.redraw()


# --------------------------------------------------------------------
# Calling function
#
def show(df):
    """
    The main function to start the data frame GUI.
    """
    app = wx.App(False)
    frame = MainFrame(df)
    frame.Show()
    app.MainLoop()


# --------------------------------------------------------------------
# Script
#

if __name__ == "__main__":

    import argparse

    csort = False
    form  = '%Y-%m-%d %H:%M:%S'
    sep   = None
    undef = None
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
A minimalistic GUI for analyzing text files.
Date/Time must be in the first data column. Its format can be given on the
command line, otherwise assumed %Y-%m-%d %H:%M:%S.''')
    hstr = 'Sort columns by name.'
    parser.add_argument('-c', '--columnsort', action='store_true',
                        default=csort, dest='csort', help=hstr)
    hstr  = "Format codes of the platform's C library strftime(),"
    hstr += " used by Python."
    parser.add_argument('-f', '--format', action='store', default=form,
                        dest='form', metavar='C-format', help=hstr)
    hstr = "Field delimiter (Default: try , ; None)"
    parser.add_argument('-s', '--sep', action='store', default=sep, dest='sep',
                        metavar='field_delimiter', help=hstr)
    hstr  = 'Missing value. Will be set to NaN for plotting.'
    hstr += ' Use long version if negative number, e.g. --undef="-9999".'
    parser.add_argument('-u', '--undef', action='store', type=float,
                        default=undef, dest='undef',
                        metavar='undef', help=hstr)
    hstr = 'Comma (,) or semi-colon (;) delimited text file(s).'
    parser.add_argument('files', nargs='*', default=None, metavar='file(s)',
                        help=hstr)

    args    = parser.parse_args()
    csort   = args.csort
    form    = args.form
    sep     = args.sep
    undef   = args.undef
    infiles = args.files

    del parser, args

    import numpy as np
    import pandas as pd

    parser = lambda date: pd.datetime.strptime(date, form)

    infile = infiles[0]
    if sep is not None:
        df = pd.read_csv(infile, sep, parse_dates=[0], date_parser=parser,
                         index_col=0, header=0)
    else:
        try:
            sep = ','
            df = pd.read_csv(infile, sep, parse_dates=[0], date_parser=parser,
                             index_col=0, header=0)
        except:
            try:
                sep = ';'
                df = pd.read_csv(infile, sep, parse_dates=[0],
                                 date_parser=parser, index_col=0, header=0)
            except:
                try:
                    # Pandas try to detect delimiter by using csv.Sniffer.
                    sep = None
                    df = pd.read_csv(infile, sep, parse_dates=[0],
                                     date_parser=parser, index_col=0, header=0)
                except:
                    raise IOError('Pandas could not read input file: '+infile)

    if len(infiles) > 1:
        for infile in infiles[1:]:
            df1 = pd.read_csv(infile, sep, parse_dates=[0], date_parser=parser,
                              index_col=0, header=0)
            df = df.append(df1, sort=True)

    # # to concatenate files
    # df.to_csv('all.csv', na_rep=-9999., line_terminator='\r\n')

    if undef is not None:
        df.replace(undef, np.nan, inplace=True)

    if csort:
        df = df.reindex(sorted(df.columns), axis=1)

    # This must be before any other call to matplotlib
    # because it use the wxAgg backend.
    # This means, do not use --pylab with ipython.
    show(df)
