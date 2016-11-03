#!/usr/bin/env python

import sys               as sys
import os                as os
import matplotlib.pyplot as plt
import numpy             as np
from jams    import date2dec as d2d
from PyQt4   import QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg    as FC
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NT

class mainwin(QtGui.QWidget):
    """
        Graphical user interface for inspecting diver data (ground water table
        and temperature) stored in *.MON files exported from Schlumberger
        DiverOffice. User can select folder containing files. Interactive plot
        is generated for one or multiple files.
        
        
        Definition
        ----------
        class mainwin(QtGui.QWidget):


        Input
        -----
        

        Optional Input
        --------------


        Output
        ------


        Restrictions
        ------------
        Works only with *.MON files exported from Schlumberger DiverOffice. 
        Typical name pattern of *.MON files is necessary to extract information:
        e.g. do_baro_4401_160509130113_H4401.MON indicating 
             site_baro or diver_serial number_crypticnumbers_serial number.MON


        Examples
        --------


        License
        -------
        This file is part of the JAMS Python package.

        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2012-2014 Matthias Cuntz


        History
        -------
        Written,  AP, Jun 2016
        Modified, AP, Sep 2016 - added multiple selection
                  MC, Nov 2016 - 00, 01, etc. for integers not accepted by Python3
    """
    def __init__(self):
        # constructor of main window
        super(mainwin, self).__init__()
        self.init_gui()
    
    def init_gui(self):
        # define window properties
        self.setWindowTitle('DiverView')
        self.setWindowIcon(QtGui.QIcon('doicon.png'))
        
        # create boxes for main window        
        self.create_file_select_box()
        self.create_plot_box()
        
        # define main window layout and include widgets
        mainLayout = QtGui.QGridLayout()
        mainLayout.addWidget(self.filebox, 1, 0)
        mainLayout.addWidget(self.plotbox, 1, 1)
        self.setLayout(mainLayout)
        
        # set size policy for runtime scaling
        self.filebox.setSizePolicy(QtGui.QSizePolicy.Fixed,     QtGui.QSizePolicy.Expanding)
        self.plotbox.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        
        # define zero date and empty data dictionary
        self.date01 = d2d(yr=1, mo=1, dy=2, hr=0, mi=0, sc=0)
        self.data   = {}
        
    def create_file_select_box(self):
        # define box
        self.filebox = QtGui.QGroupBox("Select files")
  
        # create folder and open button
        self.folder = QtGui.QLineEdit()
        button1     = QtGui.QPushButton('Open')
        button1.clicked.connect(self.selectFolder)      
          
        # create file list
        self.files     = QtGui.QListWidget()
        self.files.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        selectionModel = self.files.selectionModel()
        selectionModel.selectionChanged.connect(self.change_selection)        
  
        # define layout for file box and include widgets
        layout = QtGui.QGridLayout()
        layout.addWidget(self.folder, 0, 0)
        layout.addWidget(button1, 0, 1)
        layout.addWidget(self.files, 1, 0)
        self.filebox.setLayout(layout)
    
    def create_plot_box(self):
        # define box
        self.plotbox = QtGui.QGroupBox("Plot")

        # define plot area with toolbar and button
        self.figure  = plt.figure()
        self.canvas  = FC(self.figure)
        self.toolbar = NT(self.canvas, self)
        button1      = QtGui.QPushButton('Prev')
        button2      = QtGui.QPushButton('Next')
        button1.clicked.connect(self.plot_prev)
        button2.clicked.connect(self.plot_next)

        # define layout for plot box
        layout = QtGui.QGridLayout()
        layout.addWidget(self.toolbar,0,0,1,2)
        layout.addWidget(self.canvas,1,0,1,2)
        layout.addWidget(button1,2,0)
        layout.addWidget(button2,2,1)
        self.plotbox.setLayout(layout)

    def selectFolder(self):
        # select folder
        folder = QtGui.QFileDialog.getExistingDirectory(self, 'Open Folder')
        self.folder.setText(folder)
        
        # get files
        self.filelist = os.listdir(folder)
        self.filelist = [file for file in self.filelist if file.split('.')[-1].lower()=='mon']
        self.filelist = sorted([file for file in self.filelist if '~' not in file])
        
        for f in self.filelist:
            self.files.addItem(f)
        
        # load data
        self.data.update(load_MON(str(folder)))
        self.files.setItemSelected(self.files.item(0), True)
    
    def diver_plot(self):
        # choose selected diver data
        keys = [u'_'.join(str(x.text()).split('_')[:3]) for x in self.files.selectedItems()]
        
        # create an axis
        self.figure.clf()
        ax  = self.figure.add_subplot(111)
        ax2 = ax.twinx()

        # plot data
        for item in keys:
            ax.plot_date(self.data[item][:,0]-self.date01, self.data[item][:,1], marker=None, ls='-', label=item)
            ax2.plot_date(self.data[item][:,0]-self.date01, self.data[item][:,2], marker=None, ls='-', label=item)

        # modify plot properties
        ax.set_ylabel('Pressure (cmH2O)')
        ax2.set_ylabel('Temperature (Celsius)')
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2)
        self.figure.autofmt_xdate(rotation=45)
   
        # refresh canvas
        self.canvas.draw()
            
    def plot_next(self):
        # select next diver
        current = self.files.row(self.files.selectedItems()[-1])
        index = current + 1 if current < self.files.count() - 1 else current
        self.clear_selection()
        self.files.setItemSelected(self.files.item(index), True)
        
    def plot_prev(self):
        # select previous diver
        current = self.files.row(self.files.selectedItems()[0])
        index = current - 1 if current > 0 else current
        self.clear_selection()
        self.files.setItemSelected(self.files.item(index), True)

    def change_selection(self):
        # update plot upon selection of diver
        self.diver_plot()
    
    def clear_selection(self):
        # clear selection in diver list
        for i in range(self.files.count()):
            self.files.setItemSelected(self.files.item(i), False)    

def load_MON(path):
    """
        Loads diver data from exported *.MON files from Schlumberger DiverOffice. 
        Skips the header of 53 lines. Reads first part of file name as plot id.
        
        
        Input
        -----
        path         path ot the folder where the files are. loads all *.MON files
                     from this folder.
        
        
        Output
        ------
        data         dict with names of loggers as dict.keys and np.ma.arrays with
                     timestamp and values of respective loggers as dict.values
        
        
        Restrictions
        ------------
        Works only with *.MON files exported from Schlumberger DiverOffice. 
        Typical name pattern of *.MON files is necessary to extract information:
        e.g. do_baro_4401_160509130113_H4401.MON indicating 
             site_baro or diver_serial number_crypticnumbers_serial number.MON
    
    
        Examples
        --------
    
    
        License
        -------
        This file is part of the JAMS Python package.
    
        The JAMS Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.
    
        The JAMS Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.
    
        You should have received a copy of the GNU Lesser General Public License
        along with the JAMS Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.
    
        Copyright 2012-2014 Matthias Cuntz
    
    
        History
        -------
        Written,  AP, Jun 2016
    """
    # get files
    filelist = os.listdir(path)    
    filelist = [file for file in filelist if file.split('.')[-1].lower()=='mon']
    filelist = sorted([file for file in filelist if '~' not in file])
    
    # load data
    data = {}
    for i, file in enumerate(filelist):
        timevalues = np.genfromtxt(path+'/'+file, dtype='|S100', skip_header=53, skip_footer=1)
        asciidate  = np.array([' '.join(x).split('.')[0].replace('/','-') for x in timevalues[:,:2]])
        jd         = d2d(eng=asciidate)        
        values     = timevalues[:,2:].astype(float)
        values     = np.ma.array(values, mask=values<-9998.)
        data['_'.join(file.split('_')[:-2])] = np.ma.concatenate((jd[np.newaxis].T,values),1)
    return data

# define event loop for running the app
def main():
    app = QtGui.QApplication(sys.argv)
    mw  = mainwin()
    mw.show()
    sys.exit(app.exec_())

# run event loop
if __name__ == '__main__':
    main()
