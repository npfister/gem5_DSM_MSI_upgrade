
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

#include "defs.h"
#include "memory.h"



#ifndef _PARMACS_COMMON_INCLUDES_
#define _PARMACS_COMMON_INCLUDES_

#include <pthread.h>
#include <parmacs_config.h>
#include <parmacs_types.h>
#include <parmacs_decl.h>
#include <parmacs.h>

#endif /* _PARMACS_COMMON_INCLUDES_ */



g_mem *G_Memory;
local_memory Local[MAX_PROCS];

/*
 *  InitGlobalMemory ()
 *
 *  Args : none.
 *
 *  Returns : nothing.
 *
 *  Side Effects : Allocates all the global storage for G_Memory.
 *
 */
void
InitGlobalMemory ()
{
   int i;

   G_Memory = (g_mem *) _parmacs_memory_alloc((int)(sizeof(g_mem)));;
   G_Memory->i_array = (int *) _parmacs_memory_alloc((int)(Number_Of_Processors * sizeof(int)));;
   G_Memory->d_array = (double *) _parmacs_memory_alloc((int)(Number_Of_Processors
					 * sizeof(double)));;
   if (G_Memory == NULL) {
      printf("Ran out of global memory in InitGlobalMemory\n");
      exit(-1);
   }
   G_Memory->count = 0;
   G_Memory->id = 0;
   {
    _parmacs_lockinit(&(G_Memory->io_lock));
};
   {
    _parmacs_lockinit(&(G_Memory->mal_lock));
};
   {
    _parmacs_lockinit(&(G_Memory->single_lock));
};
   {
    _parmacs_lockinit(&(G_Memory->count_lock));
};
   {
    _parmacs_alockinit(&(G_Memory->lock_array), (MAX_LOCKS));
};
   {
    _parmacs_barinit(&(G_Memory->synch));
};
   G_Memory->max_x = -MAX_REAL;
   G_Memory->min_x = MAX_REAL;
   G_Memory->max_y = -MAX_REAL;
   G_Memory->min_y = MAX_REAL;
}


