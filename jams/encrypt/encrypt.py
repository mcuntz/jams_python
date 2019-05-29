#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
from random import randint, choice
from math import ceil, log
import os
import getpass
from jams.sendmail import sendmail

__all__ = ['file_cipher', 'file_pass', 'set_up_cipher', 'wordEncrypt', 'wordDecrypt', 'sendfail']

# modified from http://code.activestate.com/recipes/577954-encrypt-and-decrypt-text-and-text-files-beta/

file_cipher = os.path.join(os.path.expanduser('~'),'.jams.cipher')
file_pass   = os.path.join(os.path.expanduser('~'),'.jams.email.cipher')
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
        f = open(pass_file, 'r')
        cpass = f.readline()[:-1]
        f.close()
        password = wordDecrypt([ int(c) for c in cpass.split()])



        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de

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
    f = open(file_cipher, 'w')
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
        f = open(pass_file, 'r')
        cpass = f.readline()[:-1]
        f.close()
        password = wordDecrypt([ int(c) for c in cpass.split()])


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de

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
        f = open(pass_file, 'r')
        cpass = f.readline()[:-1]
        f.close()
        password = wordDecrypt([ int(c) for c in cpass.split()])


        License
        -------
        This file is part of the JAMS Python package.

        Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de

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


def sendfail(subject, message='', sender='jams.encrypt.sendfail'):
    """
        Send e-mail to e-mail addresses in file_pass, decrypting with file_cipher.


        Definition
        ----------
        sendfail(subject, message, sender='jams.encrypt.sendfail')


        Input
        -----
        subject   Subject of the e-mail
        message   Message body of e-mail


        Optional Input
        --------------
        sender      Sender name of e-mail


        Examples
        --------
        sendfail('This did not work.', sender='me@ufz.de')
        sendfail('This did not work.', 'Really.')


        License
        -------
        This file is part of the JAMS Python package.

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

        Copyright 2014-2015 Matthias Cuntz


        History
        -------
        Written,  MC, Dec 2014
        Modified, MC, Oct 2015 - transfer to ufz library; message body
    """
    f = open(file_pass, 'r')
    email = f.readline()[:-1]
    cpass = f.readline()[:-1]
    f.close()
    user = getpass.getuser()
    password = wordDecrypt([ int(c) for c in cpass.split()])
    _ = sendmail(sender,
                 to       = email,
                 subject  = subject,
                 message  = message,
                 login    = user,
                 password = password)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
