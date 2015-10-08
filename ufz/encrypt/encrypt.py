#!/usr/bin/env python
from __future__ import print_function
from random import randint, choice
from math import ceil, log
import os
import getpass
from ufz import sendmail

__all__ = ['file_cipher', 'file_pass', 'set_up_cipher', 'wordEncrypt', 'wordDecrypt', 'sendfail']

# modified from http://code.activestate.com/recipes/577954-encrypt-and-decrypt-text-and-text-files-beta/

file_cipher = os.path.join(os.path.expanduser('~'),'.ufz.cipher')
file_pass   = os.path.join(os.path.expanduser('~'),'.ufz.email.cipher')
getVar = lambda searchList, ind:  [searchList[i] for i in ind]
find   = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem]
mod    = lambda n, m:             n % m

# --------------------------------------------------------------------

def baseExpansion(n,c,b):
    '''baseExpansion - Power of ten transformation of keys'''
    i      = len(n)
    base10 = sum([pow(c,i-k-1)*n[k] for k in range(i)])
    j      = int(ceil(log(base10 + 1,b)))
    baseExpanded = [mod(base10//pow(b,j-p),b) for p in range(1,j+1)]
    return baseExpanded

# --------------------------------------------------------------------

def set_up_cipher():
    """
        Write or rewrite the cipher file: file_cipher

        This has to be done once before encryption and decryption can be used.


        Definition
        ----------
        def set_up_cipher():


        Output
        ------
        file with cipher in home directory


        Examples
        --------
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
        Written,  MC, Dec 2014 - modified from
                  http://code.activestate.com/recipes/577954-encrypt-and-decrypt-text-and-text-files-beta/
    """
    # key
    alpha = '1234567890qwertyuiop[]asdfghjkl;zxcvbnm,.!@#$%^&*()_+-=-{}:<>|QWERTYUIOPASDFGHJKLZXCVBNM ~`?'
    # cipher
    cipher = "".join([list(alpha)[randint(0,len(list(alpha))-1)] for i in range(5000)])

    if os.path.exists(file_cipher):
        os.remove(file_cipher)
    f = open(file_cipher, 'wb')
    f.write(cipher)
    f.close()

# --------------------------------------------------------------------

def wordEncrypt(word):
    """
        Encrypt a word into list of keys using the cipher in file_cipher.


        Definition
        ----------
        def wordEncrypt(word):


        Input
        -----
        word     string


        Output
        ------
        list with numeric keys


        Examples
        --------
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
        Written,  MC, Dec 2014 - modified from
                  http://code.activestate.com/recipes/577954-encrypt-and-decrypt-text-and-text-files-beta/
    """
    cipher = open(file_cipher).read()+'\n'
    cipherWord = find(cipher,list(word))
    
    keys = [randint(5001,7000), randint(2,5000)]

    encryptedWord = baseExpansion(list(map(choice, cipherWord)),keys[0],keys[1])
    encryptedWord.extend(keys)

    return list(map(int,encryptedWord))

# --------------------------------------------------------------------

def wordDecrypt(encryptedList):
    """
        Decrypt list of keys back into original word using the cipher in file_cipher.


        Definition
        ----------
        def wordDecrypt(encryptedList):


        Input
        -----
        encryptedList   list with numeric keys


        Output
        ------
        Original word (string)


        Examples
        --------
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
        Written,  MC, Dec 2014 - modified from
                  http://code.activestate.com/recipes/577954-encrypt-and-decrypt-text-and-text-files-beta/
    """
    
    cipher        = open(file_cipher).read()+'\n'
    encryptedWord = encryptedList[0:len(encryptedList)-2]
    
    keys = encryptedList[len(encryptedWord):len(encryptedList)]
    
    decryptedList = map(int,baseExpansion(encryptedWord, keys[1], keys[0]))
    
    return "".join(getVar(cipher,decryptedList))


def sendfail(subject, message, from='ufz.encrypt.sendfail'):
    """
        Send e-mail to e-mail addresses in file_pass, decrypting with file_cipher.


        Definition
        ----------
        sendfail(subject, message, from='ufz.encrypt.sendfail')


        Input
        -----
        subject   Subject of the e-mail
        message   Message body of e-mail


        Optional Input
        --------------
        from      Sender name of e-mail


        Examples
        --------
        sendfail('This did not work.', from='me@ufz.de')
        sendfail('This did not work.', 'Really.')


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

        Copyright 2014-2015 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
        Modified, MC, Oct 2015 - transfer to ufz library; message body
    """
    f = open(file_pass, 'rb')
    email = f.readline()[:-1]
    cpass = f.readline()[:-1]
    f.close()
    user = getpass.getuser()
    password = wordDecrypt([ int(c) for c in cpass.split()])
    _ = sendmail(from,
                 to       = email,
                 subject  = subject,
                 message  = message,
                 login    = user,
                 password = password)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
