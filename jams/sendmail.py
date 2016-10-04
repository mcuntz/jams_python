#!/usr/bin/env python
from __future__ import print_function
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

        Copyright 2014 Matthias Cuntz


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
