/***************************************************************************
 *   Copyright (C) 2012 by CADS/UV Software Team,                          *
 *                         Indian Institute of Astrophysics                *
 *                         Bangalore 560034                                *
 *                         cads_AT_iiap.res.in                             *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/*--------------------------------------------------------------------------
    File   : signal_handler.c
    Author : Reks
    Date   : 09/07/2005
    Version: 1.11
    File Description: This file contains a function to identify a crash
                      signal, trap it and exit with status 1.
                      the distribution. The Following functions are defined
                      in this file:
                      1. set_signal ()
                      2. signal_handler ()

    URL: More information on Cads software is available at
         http://cads.iiap.res.in/software
 *-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
    Revision         : Author, dd/mm/yyyy
    Revision history :----

    Reks, 09th Jul 2005: Created this file
    Reks, 08th Dec 2005: Re-ordered functions to stop complaints from ANSI
 *-------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

void signals_handler(int signal)
{

    /* Try to figure out what that signal means... */
    switch (signal)
    {
    case SIGSEGV:
        fprintf(stderr, "\nFATAL  : Segmentation Fault\n");
        break;
    case SIGFPE:
        fprintf(stderr, "\nFATAL  : Integer or Floating Point Exception\n");
        break;
    case SIGILL:
        fprintf(stderr, "\nFATAL  : Illegal Operation Exception\n");
        break;
    default:
        fprintf(stderr, "\nFATAL  : Unknown signal trapped. Argument = %d\n",
                signal);
        break;
    }

    /* No point continuing... print out a suggestion and packup.. */
    printf("\n-----------------------------------------------------------\n");
    printf("You seem to have come across a bug in the software. Please \n");
    printf("report this incident to http://cads.iiap.res.in/bugzilla/\n");
    printf("-----------------------------------------------------------\n\n");

    exit(EXIT_FAILURE);
}

void set_signals(void)
{

    /* Any of these signals and send them to signal_handler */
    signal(SIGSEGV, signals_handler);
    signal(SIGFPE, signals_handler);
    signal(SIGILL, signals_handler);
}

