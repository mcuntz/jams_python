#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Encrypt and decrypt text using a key system as well as a cipher.


    Constants
    ---------
    file_cipher       Name of file with cipher.
    file_pass         Name of the file with e-mail address list.


    Functions
    ---------
    set_up_cipher     Create or change cipher.
    wordEncrypt       Encrypt word into list of keys using cipher.
    wordDecrypt       Decrypts list of keys into original word using cipher.
    sendfail          Send e-mail to e-mail addresses in file_pass, decrypting with file_cipher.


    Example
    -------
    # Setup the encryption and decryption of password
    if not os.path.exists(file_cipher):
        set_up_cipher(file_cipher)

    # Store encrypted password
    if not os.path.exists(file_pass):
        print('')
        prompt = 'Enter password: '
        ipass = getpass.getpass(prompt)
        cpass = wordEncrypt(ipass)
        f = open(pass_file, 'w')
        f.write(' '.join([ str(i) for i in cpass ])+'\n')
        f.close()

    # Read encrypted password
    f = open(file_pass, 'r')
    cpass = f.readline()[:-1]
    f.close()
    password = wordDecrypt([ int(c) for c in cpass.split()])

    # Send fail e-mail
    sendfail('This test', 'Did not work', sender='me@ufz.de')


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License. The JAMS Python package originates from the former UFZ Python library,
    Department of Computational Hydrosystems, Helmholtz Centre for Environmental
    Research - UFZ, Leipzig, Germany.

    Copyright (c) 2014-2015 Matthias Cuntz - mc (at) macu (dot) de

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
    Written,  MC, Jun-Dec 2014 - modified from
              http://code.activestate.com/recipes/577954-encrypt-and-decrypt-text-and-text-files-beta/
    Modified, MC, Oct 2015     - sendfail
"""
from .encrypt import file_cipher, file_pass, set_up_cipher, wordEncrypt, wordDecrypt, sendfail

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.1'
__revision__ = "Revision: 2349"
__date__     = 'Date: 08.10.2015'
