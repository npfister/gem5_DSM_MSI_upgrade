
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

/*  This file contains the declaration of the GlobalMemory 
structure and the maximum number of molecules allowed 
by the program. */

#define MAXLCKS	4096

struct GlobalMemory {
    parmacs_lock_t (IOLock);
    parmacs_lock_t (IndexLock);
    parmacs_lock_t (IntrafVirLock);
    parmacs_lock_t (InterfVirLock);
    parmacs_lock_t (FXLock);
    parmacs_lock_t (FYLock);
    parmacs_lock_t (FZLock);
    parmacs_lock_t (KinetiSumLock);
    parmacs_lock_t (PotengSumLock);
    parmacs_lock_t (MolLock);
    parmacs_barrier_t (start);
    parmacs_barrier_t (InterfBar);
    parmacs_barrier_t (PotengBar);
    int Index;
    double VIR;
    double SUM[3];
    double POTA, POTR, POTRF;
    unsigned long createstart,createend,computestart,computeend;
    unsigned long trackstart, trackend, tracktime;
    unsigned long intrastart, intraend, intratime;
    unsigned long interstart, interend, intertime;
};

extern struct GlobalMemory *gl;

