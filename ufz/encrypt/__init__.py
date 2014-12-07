#!/usr/bin/env python
"""
    Encrypt and decrypt text using a key system as well as a cipher.


    Constants
    ---------
    cipher_file       Name of file with cipher


    Functions
    ---------
    set_up_cipher     Create or change cipher
    wordEncrypt       Encrypt word into list of keys using cipher
    wordDecrypt       Decrypts list of keys into original word using cipher


    Example
    -------
    # Setup the encryption and decryption of password
    if not os.path.exists(cipher_file):
        set_up_cipher(cipher_file)

    # Store encrypted password
    pass_file = os.path.join(os.path.expanduser('~'),'.pass.cipher')
    if not os.path.exists(pass_file):
        print('')
        prompt = 'Enter password: '
        ipass = getpass.getpass(prompt)
        cpass = wordEncrypt(ipass)
        f = open(pass_file, 'w')
        f.write(' '.join([ str(i) for i in cpass ])+'\n')
        f.close()

    # Read encrypted password
    f = open(pass_file, 'rb')
    cpass = f.readline()[:-1]
    f.close()
    password = wordDecrypt([ int(c) for c in cpass.split()])


    License
    -------
    This file is part of the UFZ Python package.

    The UFZ Python package is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The UFZ Python package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
    If not, see <http://www.gnu.org/licenses/>.

    Copyright 2014 Matthias Cuntz


    History
    -------
    Written,  MC, Jun-Dec 2014 - modified from
              http://code.activestate.com/recipes/577954-encrypt-and-decrypt-text-and-text-files-beta/
"""
from .encrypt import *

# Information
__author__   = "Matthias Cuntz"
__version__  = '1.0'
__revision__ = "Revision: 1931"
__date__     = 'Date: 08.12.2014'
