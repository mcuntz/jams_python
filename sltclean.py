import os as os
import re
import shutil
import numpy as np
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
        This file is part of the UFZ Python library.
    
        The UFZ Python library is free software: you can redistribute it and/or 
        modify it under the terms of the GNU Lesser General Public License as 
        published by the Free Software Foundation, either version 3 of the License,
        or (at your option) any later version.
    
        The UFZ Python library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.
    
        You should have received a copy of the GNU Lesser General Public License
        along with The UFZ Python library.  If not,
        see <http://www.gnu.org/licenses/>.
    
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