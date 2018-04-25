#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

__all__ = ['directory_from_gui', 'directories_from_gui', 'file_from_gui', 'files_from_gui']

# -------------------------------------------------------------------------
# Choose directories in GUI
#

def directories_from_gui(initialdir='.', title='Choose one of several directories. Press Cancel when finished.'):
    """
        Opens consecutive directory selection dialogs, returns consecutiveley selected directories


        Definition
        ----------
        def directories_from_gui(initialdir='.', title='Choose one of several directories'):


        Optional Input
        --------------
        initialdir   Initial directory GUI opens in (default: '.')
        title        Title of GUI (default: 'Choose one of several directories. Press Cancel when finished.'


        Output
        ------
        Consecutiveley selected directories


        Examples
        --------
        if not dirs:
            dirs = directories_from_gui()
            if not dirs:
                raise ValueError('Error: no directories given.')


        History
        -------
        Written,  MC, Oct 2015
    """
    try:                # Python 3
        import tkinter as Tkinter
        import tkinter.filedialog as tkFileDialog
    except ImportError: # Python 2
        import Tkinter, tkFileDialog
    root = Tkinter.Tk()
    root.withdraw()                                      # hide root window, i.e. white square

    # always on top
    root.tk.call('wm', 'attributes', '.', '-topmost', 1) # focus on (hidden) window so that child is on top

    # # Make it almost invisible - no decorations, 0 size, top left corner.
    # # Then show window again and lift it to top so it can get focus
    # root.overrideredirect(True)
    # root.geometry('0x0+0+0')
    # root.deiconify()
    # root.lift()
    # root.focus_force()

    idir = initialdir
    alldirs = []
    while True:
        dirs = tkFileDialog.askdirectory(parent=root,
                                         title=title,
                                         initialdir=idir)
        if not dirs: break
        alldirs.append(dirs)
        idir = dirs

    root.destroy()

    return alldirs

# -------------------------------------------------------------------------
# Choose one directory in GUI
#

def directory_from_gui(initialdir='.', title='Choose directory'):
    """
        Opens directory selection dialog, returns selected directory


        Definition
        ----------
        def directory_from_gui(initialdir='.', title='Choose one directory'):


        Optional Input
        --------------
        initialdir   Initial directory GUI opens in (default: '.')
        title        Title of GUI (default: 'Choose directory'


        Output
        ------
        Selected directory


        Examples
        --------
        if not dirs:
            dirs = directory_from_gui()
            if not dirs:
                raise ValueError('Error: no directory given.')


        History
        -------
        Written,  MC, Jun 2014
    """
    try:                # Python 3
        import tkinter as Tkinter
        import tkinter.filedialog as tkFileDialog
    except ImportError: # Python 2
        import Tkinter, tkFileDialog
    root = Tkinter.Tk()
    root.withdraw()                                      # hide root window, i.e. white square

    # always on top
    root.tk.call('wm', 'attributes', '.', '-topmost', 1) # focus on (hidden) window so that child is on top

    # # Make it almost invisible - no decorations, 0 size, top left corner.
    # # Then show window again and lift it to top so it can get focus
    # root.overrideredirect(True)
    # root.geometry('0x0+0+0')
    # root.deiconify()
    # root.lift()
    # root.focus_force()

    dirs = tkFileDialog.askdirectory(parent=root,
                                     title=title,
                                     initialdir=initialdir)

    root.destroy()

    return dirs



# -------------------------------------------------------------------------
# Choose one file in GUI
#

def file_from_gui(initialdir='.', title='Choose file', multiple=False):
    """
        Wrapper for files_from_gui with multiple=False as default


        Definition
        ----------
        def file_from_gui(initialdir='.', title='Choose file', multiple=False):


        Optional Input
        --------------
        initialdir   Initial directory GUI opens in (default: '.')
        title        Title of GUI (default: 'Choose file'
        multiple     True:  allow selection of multiple files
                     False: only one single file to select (default)


        Output
        ------
        List of selected files


        Restrictions
        ------------
        Always returns a list even with multiple=False.


        Examples
        --------
        if not file:
            file = file_from_gui()
            if not file:
                raise ValueError('Error: no input file given.')


        History
        -------
        Written,  MC, Jun 2014
    """
    return files_from_gui(initialdir=initialdir, title=title, multiple=multiple)


# -------------------------------------------------------------------------
# Choose files in GUI
#

def files_from_gui(initialdir='.', title='Choose file(s)', multiple=True):
    """
        Opens file selection dialog, returns selected files


        Definition
        ----------
        def files_from_gui(initialdir='.', title='Choose file(s)', multiple=True):


        Optional Input
        --------------
        initialdir   Initial directory GUI opens in (default: '.')
        title        Title of GUI (default: 'Choose file(s)'
        multiple     True:  allow selection of multiple files (default)
                     False: only one single file to select


        Output
        ------
        List of selected files


        Restrictions
        ------------
        Always returns a list even with multiple=False.


        Examples
        --------
        if not files:
            files = files_from_gui()
            if not files:
                raise ValueError('Error: no input file(s) given.')


        History
        -------
        Written,  MC, Jun 2014
    """
    try:                # Python 3
        import tkinter as Tkinter
        import tkinter.filedialog as tkFileDialog
    except ImportError: # Python 2
        import Tkinter, tkFileDialog
    root = Tkinter.Tk()
    root.withdraw()                                      # hide root window, i.e. white square

    # always on top
    root.tk.call('wm', 'attributes', '.', '-topmost', 1) # focus on (hidden) window so that child is on top

    # # Make it almost invisible - no decorations, 0 size, top left corner.
    # # Then show window again and lift it to top so it can get focus
    # root.overrideredirect(True)
    # root.geometry('0x0+0+0')
    # root.deiconify()
    # root.lift()
    # root.focus_force()

    files = tkFileDialog.askopenfilename(parent=root,
                                         title=title,
                                         multiple=multiple,
                                         initialdir=initialdir)
    files = list(root.tk.splitlist(files))
    
    root.destroy()

    return files


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
