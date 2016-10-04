#!/usr/bin/env python
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
    f = open(file_pass, 'rb')
    cpass = f.readline()[:-1]
    f.close()
    password = wordDecrypt([ int(c) for c in cpass.split()])

    # Send fail e-mail
    sendfail('This test', 'Did not work', sender='me@ufz.de')


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

    Copyright 2014-2015 Matthias Cuntz


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
