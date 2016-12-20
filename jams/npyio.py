#!/usr/bin/env python
# from __future__ import absolute_import
import numpy as np
import numpy.lib.format as format
from numpy.compat import basestring
import os

__all__ = ['savez', 'savez_compressed']


def savez(file, *args, **kwds):
    """
    Save several arrays into a single file in uncompressed ``.npz`` format.

    If arguments are passed in with no keywords, the corresponding variable
    names, in the ``.npz`` file, are 'arr_0', 'arr_1', etc. If keyword
    arguments are given, the corresponding variable names, in the ``.npz``
    file will match the keyword names.

    This is an extension of numpy.savez with a possible keyword arguments ``append``,
    which appends arrays in an already existing npz-file, ``update``, which also updates
    existing members of the zip archive, and ``compress`` which produces compressed npz
    files and is the same as savez_compressed.

    Note that ``update`` copies existing archive members to a new zip (just as the zip utility).
    It therefore needs doubel disk space temporarily.

    Also note that file must be a filename and no file handle, contrary to numpy.savez.

    Parameters
    ----------
    file : str
        File name (string) where the data will be saved. The ``.npz``
        extension will be appended to the file name if it is not already there.
    args : Arguments, optional
        Arrays to save to the file. Since it is not possible for Python to
        know the names of the arrays outside `savez`, the arrays will be saved
        with names "arr_0", "arr_1", and so on. These arguments can be any
        expression.
    kwds : Keyword arguments, optional
        Arrays to save to the file. Arrays will be saved in the file with the
        keyword names.
    append : Keyword argument, optional
        True = append to existing ``.npz`` file
        False = overwrite possibly existing ``.npz`` file (default).
    update : Keyword argument, optional
        True = update existing members of a zip archive, i.e. ``.npz`` file.
        False = raises ValueError if archive member already exists in ``.npz`` file (default).
    compress : Keyword argument, optional
        True = produce compressed ``.npz`` file; same as savez_compressed.
        False = produce uncompressed ``.npz``  (default).

    Returns
    -------
    None

    See Also
    --------
    savez_compressed : Save several arrays into a compressed ``.npz`` archive

    Notes
    -----
    The ``.npz`` file format is a zipped archive of files named after the
    variables they contain.  The archive is not compressed and each file
    in the archive contains one variable in ``.npy`` format. For a
    description of the ``.npy`` format, see `format`.

    When opening the saved ``.npz`` file with `load` a `NpzFile` object is
    returned. This is a dictionary-like object which can be queried for
    its list of arrays (with the ``.files`` attribute), and for the arrays
    themselves.

    Examples
    --------
    >>> import numpy as np
    >>> from tempfile import mkstemp
    >>> fd, outfile = mkstemp('')
    >>> os.close(fd)
    >>> outfile = outfile + '.npz'
    >>> x = np.arange(10)
    >>> y = np.sin(x)

    Using `savez` with \\*args, the arrays are saved with default names.

    >>> savez(outfile, x, y)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['arr_0', 'arr_1']
    >>> npzfile['arr_0']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> npzfile.close()

    Using `savez` with \\**kwds, the arrays are saved with the keyword names.

    >>> savez(outfile, x=x, y=y)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['x', 'y']
    >>> npzfile['x']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> npzfile.close()

    Using `savez` with compress, append, and update keywords.

    >>> savez(outfile, x=x, y=y, compress=True)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['x', 'y']
    >>> npzfile['x']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> npzfile.close()
    >>> x2 = 2*x
    >>> y2 = 2*y
    >>> savez(outfile, x2=x2, y2=y2, append=True)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['x', 'x2', 'y', 'y2']
    >>> npzfile['x2']
    array([ 0,  2,  4,  6,  8, 10, 12, 14, 16, 18])
    >>> npzfile.close()
    >>> savez(outfile, x=x2, update=True)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['x', 'x2', 'y', 'y2']
    >>> npzfile['x']
    array([ 0,  2,  4,  6,  8, 10, 12, 14, 16, 18])
    >>> npzfile.close()

    >>> import os
    >>> os.remove(outfile)

    """
    _savez(file, args, kwds, False)


def savez_compressed(file, *args, **kwds):
    """
    Save several arrays into a single file in compressed ``.npz`` format.

    If keyword arguments are given, then filenames are taken from the keywords.
    If arguments are passed in with no keywords, then stored file names are
    arr_0, arr_1, etc.

    This is an extension of numpy.savez_compressed with a possible keyword arguments ``append``,
    which appends arrays in an already existing npz-file, and ``update``, which also updates
    existing members of the zip archive.

    Note that ``update`` copies existing archive members to a new zip (just as the zip utility).
    It therefore needs doubel disk space temporarily.

    Also note that file must be a filename and no file handle, contrary to numpy.savez.

    Parameters
    ----------
    file : str
        File name (string) where the data will be saved. The ``.npz``
        extension will be appended to the file name if it is not already there.
    args : Arguments, optional
        Arrays to save to the file. Since it is not possible for Python to
        know the names of the arrays outside `savez`, the arrays will be saved
        with names "arr_0", "arr_1", and so on. These arguments can be any
        expression.
    kwds : Keyword arguments, optional
        Arrays to save to the file. Arrays will be saved in the file with the
        keyword names.
    append : Keyword argument, optional
        True = append to existing ``.npz`` file
        False = overwrite possibly existing ``.npz`` file (default).
    update : Keyword argument, optional
        True = update existing members of a zip archive, i.e. ``.npz`` file.
        False = raises ValueError if archive member already exists in ``.npz`` file (default).

    See Also
    --------
    savez : Save several arrays into an uncompressed ``.npz`` file format

    Notes
    -----
    The ``.npz`` file format is a zipped archive of files named after the
    variables they contain. The archive is compressed and each file
    in the archive contains one variable in ``.npy`` format. For a
    description of the ``.npy`` format, see `format`.

    When opening the saved ``.npz`` file with `load` a `NpzFile` object is
    returned. This is a dictionary-like object which can be queried for
    its list of arrays (with the ``.files`` attribute), and for the arrays
    themselves.

    Examples
    --------
    >>> import numpy as np
    >>> from tempfile import mkstemp
    >>> fd, outfile = mkstemp('')
    >>> os.close(fd)
    >>> outfile = outfile + '.npz'
    >>> x = np.arange(10)
    >>> y = np.sin(x)

    Using `savez` with \\*args, the arrays are saved with default names.

    >>> savez_compressed(outfile, x, y)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['arr_0', 'arr_1']
    >>> npzfile['arr_0']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> npzfile.close()

    Using `savez` with \\**kwds, the arrays are saved with the keyword names.

    >>> savez_compressed(outfile, x=x, y=y)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['x', 'y']
    >>> npzfile['x']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> npzfile.close()

    Using `savez` with append and update keywords.

    >>> savez_compressed(outfile, x=x, y=y)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['x', 'y']
    >>> npzfile['x']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> npzfile.close()
    >>> x2 = 2*x
    >>> y2 = 2*y
    >>> savez_compressed(outfile, x2=x2, y2=y2, append=True)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['x', 'x2', 'y', 'y2']
    >>> npzfile['x2']
    array([ 0,  2,  4,  6,  8, 10, 12, 14, 16, 18])
    >>> npzfile.close()
    >>> savez_compressed(outfile, x=x2, update=True)
    >>> npzfile = np.load(outfile)
    >>> list(np.sort(npzfile.files))
    ['x', 'x2', 'y', 'y2']
    >>> npzfile['x']
    array([ 0,  2,  4,  6,  8, 10, 12, 14, 16, 18])
    >>> npzfile.close()

    >>> import os
    >>> os.remove(outfile)
    """
    _savez(file, args, kwds, True)


def _savez(file, args, kwds, compress):
    # Import is postponed to here since zipfile depends on gzip, an optional
    # component of the so-called standard library.
    import zipfile
    import tempfile
    import shutil

    if isinstance(file, basestring):
        if not file.endswith('.npz'):
            file = file + '.npz'

    namedict = kwds
    for i, val in enumerate(args):
        key = 'arr_%d' % i
        if key in namedict.keys():
            raise ValueError("Cannot use un-named variables and keyword %s" % key)
        namedict[key] = val

    if compress:
        compression = zipfile.ZIP_DEFLATED
    else:
        compression = zipfile.ZIP_STORED

    # use compress keyword if given; only active in savez
    if 'compress' in namedict.keys():
        if namedict['compress']:
            compression = zipfile.ZIP_DEFLATED
        del namedict['compress']

    # append or update
    if 'append' in namedict.keys():
        appendit = namedict['append']
        del namedict['append']
    else:
        appendit = False
    if 'update' in namedict.keys():
        updateit = namedict['update']
        del namedict['update']
    else:
        updateit = False
    if appendit and updateit:
        raise KeyError("append and update mutually exclusive.")

    # check if file exists, otherwise it will be a simple write
    if not os.path.isfile(file):
        appendit = False
        updateit = False
        inzipf   = []
    else:
        zipf   = zipfile.ZipFile(file, mode="r")
        inzipf = zipf.namelist()
        inzipf = [ i[:-4] for i in inzipf ]
        zipf.close()
    allkeys = set(namedict.keys())
    allkeys.update(inzipf)

    # append if keyword is True
    if appendit:
        mode = "a"
        # check if new arrays already exist in zip file
        if len(inzipf) != 0:
            for key in namedict:
                if key in inzipf:
                    raise ValueError("array name already in npz-file: %s" % key)
    else:
        mode = "w"

    if not updateit:
        # Just add arrays to existing or non-existing file; duplicates were checked before
        zipf = zipfile.ZipFile(file, mode=mode, compression=compression, allowZip64=True)
        # Stage arrays in a temporary file on disk, before writing to zip.
        fd, tmpfile = tempfile.mkstemp(suffix='-numpy.npy')
        os.close(fd)
        try:
            for key, val in namedict.items():
                fname = key + '.npy'
                fid = open(tmpfile, 'wb')
                try:
                    format.write_array(fid, np.asanyarray(val))
                    fid.close()
                    fid = None
                    zipf.write(tmpfile, arcname=fname)
                finally:
                    if fid:
                        fid.close()
        finally:
            os.remove(tmpfile)
        zipf.close()
    else:
        # open existing zip file in read mode
        zipr = zipfile.ZipFile(file, mode="r")
        # open temporary zip file in write mode
        tempdir = tempfile.mkdtemp()
        try:
            tempname = os.path.join(tempdir, 'new.zip')
            zipw = zipfile.ZipFile(tempname, mode="w", compression=compression, allowZip64=True)
            for key in allkeys:
                # if in namedict then write new, else extract it from zipfile
                if key in namedict.keys():
                    # Stage arrays in a temporary file on disk, before writing to zip.
                    fd, tmpfile = tempfile.mkstemp(suffix='-numpy.npy')
                    os.close(fd)
                    try:
                        fname = key + '.npy'
                        fid = open(tmpfile, 'wb')
                        try:
                            format.write_array(fid, np.asanyarray(namedict[key]))
                            fid.close()
                            fid = None
                            zipw.write(tmpfile, arcname=fname)
                        finally:
                            if fid:
                                fid.close()
                    finally:
                        os.remove(tmpfile)
                else:
                    fname = key + '.npy'
                    zipr.extract(fname, tempdir)
                    tmpfile = os.path.join(tempdir, fname)
                    zipw.write(tmpfile, arcname=fname)
                    os.remove(tmpfile)
            # close both files and move new to old
            zipr.close()
            zipw.close()
            shutil.move(tempname, file)
        finally:
            shutil.rmtree(tempdir)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
