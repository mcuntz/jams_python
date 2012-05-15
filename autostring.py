#!/usr/bin/env python
import numpy as np

def autostring(num, prec=0, zero=False, set_printoptions=False, pp=False, join=False, joinall=False, sep=' '):
    """
        Format number (array) with given decimal precision.

        Definition
        ----------
        def autostring(num, prec=0, zero=False, set_printoptions=False, pp=False, join=False, joinall=False, sep=' '):
          There is a wrapper function for convenience with the short name 'astr' that calls autostring
        def astr(num, prec=0, zero=False, set_printoptions=False, pp=False, join=False, joinall=False, sep=' '):


        Input
        -----
        num                 number array


        Optional Input
        --------------
	prec                number of decimal places of formatted values
	                    minimum field width for integers (default: 0)
	zero                if True, pad values with zeros rather than blanks (default: False)
	set_printoptions    if True, sets linewidth to the format times size of 1st dimension (default: False)
	pp                  shortcut for set_printoptions (default: False)
                            it will be checked for (pp | set_printoptions)
	join                if True, joins all individual strings of last (fastest) dimension into one string (default: False)
	joinall             if True, joins all individual strings into single string,
                            i.e. first flattens the array and then joins it (default: False, overwrites join)
	sep                 separator used when joining (default: space=' ')


        Output
        ------
	string (array) of formatted numbers


        Restrictions
        ------------
	None


        Examples
        --------
        >>> print autostring(3.5967, 3)
        3.597

        >>> print autostring(3.5967)
        4

        >>> print autostring(3, 3)
          3

        >>> print autostring(np.array([3.5967, 3.5964]), 3)
	['3.597' '3.596']

        >>> print autostring(np.array([3.59, 1.123456e12]), 3)
	['3.590e+00' '1.123e+12']

        >>> print autostring(np.array([3.59, 11.1234]), 3, zero=True)
	['03.590' '11.123']

        >>> print autostring(np.array([3, 11]))
	[' 3' '11']

        >>> print autostring(np.array([3, 11]), 3)
	['  3' ' 11']

        >>> print autostring(np.zeros((2,2), dtype=np.float), 1)
	[['0.0' '0.0']
         ['0.0' '0.0']]

	>>> np.set_printoptions(threshold=10)
        >>> print autostring(np.zeros((2,10), dtype=np.float), 1)
	[['0.0' '0.0' '0.0' ..., '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' ..., '0.0' '0.0' '0.0']]

        >>> print autostring(np.zeros((2,10), dtype=np.float), 1, set_printoptions=True)
	[['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print autostring(np.zeros((2,10), dtype=np.float), 1, pp=True)
	[['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print autostring(np.zeros((2,10), dtype=np.float), 1, set_printoptions=False, pp=True)
	[['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print autostring(np.array([3.5967, 3.5964]), 3, join=True)
        3.597 3.596

        >>> print autostring(np.zeros((2,10), dtype=np.float), 1, join=True, sep=';')
	['0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0'
	 '0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0']

        >>> print autostring(np.reshape(np.arange(20,dtype=np.float),(2,10)), 1, joinall=True, sep=';')
         0.0; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 9.0;10.0;11.0;12.0;13.0;14.0;15.0;16.0;17.0;18.0;19.0

        >>> print astr(np.reshape(np.arange(20,dtype=np.float),(2,10)), 1, joinall=True, sep=';')
         0.0; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 9.0;10.0;11.0;12.0;13.0;14.0;15.0;16.0;17.0;18.0;19.0


        History
        -------
        Written,  MC, Nov 2011 - from autostring.pro
        Modified, MC, May 2012 - pp
    """
    #
    # Check input
    isarr = np.size(np.shape(num))
    if (isarr > 2):
	print "AUTOSTRING WARNING: autostring only works with scalars, 1D- and 2D arrays: return original array."
	return num
    # Only treat int and float
    if (isarr==0):
	try:
	    typ = num.dtype
	except AttributeError:
	    if (type(num) == float):
		typ = np.float64
	    elif (type(num) == int):
		typ = np.int32
    else:
	typ = num.dtype
    if np.__version__ >= "1.6":
        if (typ in [np.float16, np.float32, np.float64, np.float128]):
            isfloat = True
        elif (typ in [np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64]):
            isfloat = False
        else:
            print "AUTOSTRING WARNING: autostring cannot work with input type: return original array."
            return num
    else:
        if (typ in [np.float32, np.float64, np.float128]):
            isfloat = True
        elif (typ in [np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64]):
            isfloat = False
        else:
            print "AUTOSTRING WARNING: autostring cannot work with input type: return original array."
            return num        
    # Scalar to array if necessary; Special treatment of -0.0
    if (isarr==0):
	if (num == 0):
	    num = num - num
    else:
	num = np.where(num == 0, 0, num)
    # Zero padding
    if zero:
	nix = '0'
    else:
	nix = ''
    #
    # If we deal with an array of numbers we take the largest for the format
    abs_num = np.amax(np.abs(np.array(num)))
    # leave room for the decimal point and the negative sign, if any
    if (np.amin(num) < 0.):
	num_sign_chars = 1
    else:
	num_sign_chars = 0
    #
    # Floating point
    if isfloat: # number is a float, more or less
	if abs_num >= 1.e6:
	    num_prefix_chars  = 1
	    num_sci_not_chars = 4
	    format_type       = 'e'
	elif ((abs_num < 1.e6) & (abs_num >= 1.)):
	    nprefix = np.int_(np.log10(np.int32(abs_num)))+1
            # special treatment: the output prefix digits could
	    # be one digit longer as the input prefix digits: e.g. 99.99 => 100.0
	    val               = np.around(abs_num*(10.**prec))/(10.**prec)
	    nprefixval        = np.int_(np.log10(val))+1
	    nprefix           = np.amax(np.array([nprefix,nprefixval], dtype=np.int))
	    num_prefix_chars  = nprefix
	    num_sci_not_chars = 0
	    format_type       = 'f'
	elif ((abs_num < 1.) & (abs_num >= 1.e-3)):
	    num_prefix_chars  = 1
	    num_sci_not_chars = 0
	    format_type       = 'f'
	elif (abs_num == 0):
	    num_prefix_chars  = 1
	    num_sci_not_chars = 0
	    format_type       = 'f'
	else:
	    num_prefix_chars  = 1
	    num_sci_not_chars = 4
	    format_type       = 'e'
	#
	num_postfix_chars = prec
	num_total_chars   = num_sign_chars + num_prefix_chars + 1 + num_postfix_chars + num_sci_not_chars
	if (prec == 0): # no dot if prec=0
	    num_total_chars -= 1
	format_string     = ("{0:s}{1:s}{2:d}{3:s}{4:d}{5:s}{6:s}".format('{0:', nix, num_total_chars, 
									  '.', num_postfix_chars, format_type, '}'))
    else: # number is an integer
	format_type = 'd'
	if abs_num != 0:
	    num_digits = np.int_(np.log10(abs_num))+1
	else:
	    num_digits = 1 
	num_total_chars = np.maximum(num_digits + num_sign_chars, prec)
	format_string = ("{0:s}{1:s}{2:d}{3:s}{4:s}".format('{0:', nix, num_total_chars, format_type, '}'))
    #
    if (isarr == 0):
	out = format_string.format(num)
    else:
	fnum = num.flatten()
	nnum = np.size(fnum)
	styp = 'S{0:d}'.format(num_total_chars)
	out = np.empty(nnum, dtype=styp)
	for i in xrange(nnum):
	    out[i] = format_string.format(fnum[i])
	out = np.reshape(out, np.shape(num))
	if (set_printoptions | pp):
	    # num_total_chars+3 for '' and space, +isarr for []
	    np.set_printoptions(linewidth=np.size(num,-1)*(num_total_chars+3)+isarr, threshold=nnum+1)
	if (join | joinall): # There should be reduction routines in numpy
	    if ((isarr == 1) | ((isarr==2) & joinall)):
		if (isarr == 2):
		    out = out.flatten()
		for i in xrange(np.size(out)):
		    if (i==0):
			outc = out[i]
		    else:
			outc = outc+sep+out[i]
	    else:
		sform = 'S{0:d}'.format((len(out[0,0])+len(sep))*np.size(out,1))
		outc = np.zeros(np.size(out,0), dtype=sform)
		for j in xrange(np.size(out,0)):
		    for i in xrange(np.size(out,1)):
			if (i==0):
			    outc[j] = out[j,i]
			else:
			    outc[j] = outc[j]+sep+out[j,i]
	    out = outc
    # return formatted string
    return out

def astr(num, prec=0, zero=False, set_printoptions=False, pp=False, join=False, joinall=False, sep=' '):
    """
        Wrapper function for autostring
        def astr(num, prec=0, zero=False, set_printoptions=False, pp=False, join=False, joinall=False, sep=' '):
    """
    return autostring(num, prec, zero, set_printoptions, pp, join, joinall, sep)

 
if __name__ == '__main__':
    import doctest
    doctest.testmod()
