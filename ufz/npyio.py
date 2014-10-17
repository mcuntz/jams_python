#!/usr/bin/env python
# from __future__ import absolute_import

import numpy as np
import numpy.lib.format as format
from numpy.compat import basestring
import os

__all__ = ['save', 'savez', 'savez_compressed']


def zipfile_factory(*args, **kwargs):
    import zipfile
    kwargs['allowZip64'] = True
    return zipfile.ZipFile(*args, **kwargs)


def savez(file, *args, **kwds):
    """
    Save several arrays into a single file in uncompressed ``.npz`` format.

    If arguments are passed in with no keywords, the corresponding variable
    names, in the ``.npz`` file, are 'arr_0', 'arr_1', etc. If keyword
    arguments are given, the corresponding variable names, in the ``.npz``
    file will match the keyword names.

    This is a copy of numpy.savez but with a possible keyword argument ``append``,
    which appends arrays in an already existing npz-file. Also the keyword ``append``,
    which produces compressed instead of uncompressed npz files (same as savez+compressed).
    With this mechanism, compressed and uncompressed arrays can be stored in a then
    mixed npz file.

    Parameters
    ----------
    file : str or file
        Either the file name (string) or an open file (file-like object)
        where the data will be saved. If file is a string, the ``.npz``
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
        False = overwrite possible existing ``.npz`` file (default).
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
    >>> from tempfile import TemporaryFile, mkstemp
    >>> outfile = TemporaryFile()
    >>> x = np.arange(10)
    >>> y = np.sin(x)

    Using `savez` with \\*args, the arrays are saved with default names.

    >>> savez(outfile, x, y)
    >>> outfile.seek(0) # Only needed here to simulate closing & reopening file
    >>> npzfile = np.load(outfile)
    >>> npzfile.files
    ['arr_1', 'arr_0']
    >>> npzfile['arr_0']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    Using `savez` with \\**kwds, the arrays are saved with the keyword names.

    >>> outfile = TemporaryFile()
    >>> savez(outfile, x=x, y=y)
    >>> outfile.seek(0)
    >>> npzfile = np.load(outfile)
    >>> npzfile.files
    ['y', 'x']
    >>> npzfile['x']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    Using `savez` with append keyword

    >>> outfile, noutfile = mkstemp()
    >>> noutfile = noutfile + '.npz'
    >>> savez(noutfile, x=x, y=y, compress=True)
    >>> npzfile = np.load(noutfile)
    >>> npzfile.files
    ['y', 'x']
    >>> npzfile['x']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> npzfile.close()
    >>> x2 = 2*x
    >>> y2 = 2*y
    >>> savez(noutfile, x2=x2, y2=y2, append=True)
    >>> npzfile = np.load(noutfile)
    >>> npzfile.files
    ['y', 'x', 'x2', 'y2']
    >>> npzfile['x2']
    array([ 0,  2,  4,  6,  8, 10, 12, 14, 16, 18])
    >>> npzfile.close()
    >>> import os
    >>> os.remove(noutfile)
    
    """
    _savez(file, args, kwds, False)


def savez_compressed(file, *args, **kwds):
    """
    Save several arrays into a single file in compressed ``.npz`` format.

    If keyword arguments are given, then filenames are taken from the keywords.
    If arguments are passed in with no keywords, then stored file names are
    arr_0, arr_1, etc.

    This is a copy of numpy.savez_compressed but with a possible keyword argument ``append``,
    which appends arrays in an already existing npz-file.

    Parameters
    ----------
    file : str
        File name of ``.npz`` file.
    args : Arguments
        Function arguments.
    kwds : Keyword arguments
        Keywords.
    append : True = append to existing ``.npz`` file; False = overwrite possible existing ``.npz`` file.

    See Also
    --------
    savez : Save several arrays into an uncompressed ``.npz`` file format

    """
    _savez(file, args, kwds, True)


def _savez(file, args, kwds, compress):
    # Import is postponed to here since zipfile depends on gzip, an optional
    # component of the so-called standard library.
    import zipfile
    # Import deferred for startup time improvement
    import tempfile

    if isinstance(file, basestring):
        if not file.endswith('.npz'):
            file = file + '.npz'

    namedict = kwds
    for i, val in enumerate(args):
        key = 'arr_%d' % i
        if key in namedict.keys():
            raise ValueError(
                "Cannot use un-named variables and keyword %s" % key)
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

    # append if keyword is True
    mode = "w"
    if 'append' in namedict.keys():
	if namedict['append']:
	    mode = "a"
	del namedict['append']
    zipf = zipfile_factory(file, mode=mode, compression=compression)

    # check if new arrays already exist in file in append mode
    inzipf = zipf.namelist()
    if len(inzipf) != 0:
	for key in namedict:
	    fname = key + '.npy'
	    if fname in inzipf:
		raise ValueError("array name already in npz-file: %s" % key)

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


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
