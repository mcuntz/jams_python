#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import os as os
import re
import shutil
from time import localtime
from scipy.stats import mode
                                    
def sltclean(indir, pat = '[a-zA-Z0-9]*.slt|[a-zA-Z0-9]*.SLT'):
    """       
        Moves *.slt files to a "deleted" folder to exclude them from further
        processing if they have a file size smaller than half of the regular
        file size. Regular file size is determined by mode(all file sizes in the
        folder). *.slt files are raw eddy covariance files (binary) recorded
        with EddyMeas (Kolle & Rebmann, 2007)
        
        
        Definition
        ----------
        sltclean(indir, pat = '[a-zA-Z0-9]*.slt|[a-zA-Z0-9]*.SLT'):
        
        
        Input
        ----- 
        indir       str, path of the folder containing the *.slt files 
        
        
        Optional Input
        --------------
        pat         str, regular expression, describing the name pattern of
                    the *.slt files in the indir folder
                    
        
        Output
        ------
        sltclean_X_X.log log file of the cleaning process
        
        
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
            
        Copyright 2014 Arndt Piayda
    
    
        History
        -------
        Written,  AP, Jul 2014
        
    """
    
    ###########################################################################
    # reading input directory
    dirlist = os.listdir(indir)
    sizelist, filelist = np.array([]), np.array([])
    if 'deleted' not in dirlist:
        os.mkdir('%s/deleted'%indir)
    
    ###########################################################################
    # remove all files and folders from list which are not *.slt files and get size
    pat = re.compile(pat)
    for item in dirlist:
        if re.search(pat, item):
            sizelist = np.append(sizelist, os.path.getsize('%s/%s' %(indir, item))/1000.)
            filelist = np.append(filelist, item)
    filesize = mode(sizelist)[0][0]
    
    ###########################################################################
    # move files to deleted which are too small and write log file
    delfiles = filelist[sizelist<filesize/2.]

    log = open('%s/deleted/sltclean%04i%02i%02i_%02i%02i%02i.log'\
               %((indir,)+localtime()[:6]), 'w')
    log.write('Regular file size: %i\n'%filesize)
    log.write('Moved to deleted:\n')
    for item in delfiles:
        shutil.move('%s/%s'%(indir,item), '%s/deleted'%indir)
        log.write('%s\n'%item)
    log.close()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
