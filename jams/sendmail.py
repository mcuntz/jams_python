#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import smtplib

__all__ = ['sendmail']

def sendmail(fromm, to, subject, message, login, password, cc=None, smtpserver='smtp.ufz.de:587'):
    """
        Send an e-mail.

        The default smtp-server is smtp.ufz.de on port 587 (SSL/START-TLS).


        Definition
        ----------
        def sendmail(fromm, to, subject, message, login, password, cc=None, smtpserver='smtp.ufz.de:587'):


        Input
        -----
        fromm        string, Name or e-mail for From field.
        to           string, Comma separated list of e-mail addresses.
        subject      string, Subject line
        message      string, Message text
        login        login at smtp-server
        password     password for smtp-server


        Optional input parameters
        -------------------------
        cc           string, Comma separated list of e-mail addresses.
        smtpserver   string, SMTP-server:port (default: smtp.ufz.de:587)


        Output
        ------
        Return of smtplib sendmail function


        Examples
        --------
        import getpass
        email = input('Enter e-mail address in (double) quotes (\' \' or " "): ')
        user = getpass.getuser()
        prompt = 'Enter password for user '+user+': '
        ipass = getpass.getpass(prompt)

        a = sendmail('sendmail@ufz.de',
                     to       = email,
                     subject  = 'sendmail test',
                     message  = '''
                                This is the actual message.
                                It can have several lines, etc.

                                Regards
                                Matthias
                                ''',
                    login    = user,
                    password = ipass)
        print('Sendmail return: ', a ,'.')


        Notes
        -----
        modified from http://rosettacode.org/wiki/Send_an_email#Python


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

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


        History
        -------
        Written,  MC, Dec 2014
    """

# --------------------------------------------------------------------
# Send E-Mail
# modified from http://rosettacode.org/wiki/Send_an_email#Python
    header  = 'From: {:s}\n'.format(fromm)
    if type(to) == list:
        too     = to                                 # to field in sendmail needs list
        header += 'To: {:s}\n'.format(','.join(to))  # To: in message needs string
    else:
        too     = to.split(',')
        header += 'To: {:s}\n'.format(to)
    if (cc is not None) and (cc != ''):
        if type(cc) == list:
            header += 'Cc: {:s}\n'.format(','.join(cc))
        else:
            header += 'Cc: {:s}\n'.format(cc)
    header += 'Subject: {:s}\n\n'.format(subject)
    message = header + message

    server = smtplib.SMTP(smtpserver)
    server.starttls()
    server.login(login, password)
    problems = server.sendmail(fromm, too, message)
    server.quit()
    return problems


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

#     import getpass
#     email = input('Enter e-mail address in (double) quotes (\' \' or " "): ')
#     user = getpass.getuser()
#     prompt = 'Enter password for user '+user+': '
#     ipass = getpass.getpass(prompt)

#     a = sendmail('sendmail@ufz.de',
#                  to       = email,
#                  subject  = 'sendmail test',
#                  message  = '''
# This is the actual message.
# It can have several lines, etc.

# Regards
# Matthias
#                             ''',
#                 login    = user,
#                 password = ipass)
#     print('Sendmail return: ', a ,'.')
