
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

/*************************************************************************/
/*                                                                       */
/*  SPLASH Ocean Code                                                    */
/*                                                                       */
/*  This application studies the role of eddy and boundary currents in   */
/*  influencing large-scale ocean movements.  This implementation uses   */
/*  dynamically allocated four-dimensional arrays for grid data storage. */
/*                                                                       */
/*  Command line options:                                                */
/*                                                                       */
/*     -nN : Simulate NxN ocean.  N must be (power of 2)+2.              */
/*     -pP : P = number of processors.  P must be power of 2.            */
/*     -eE : E = error tolerance for iterative relaxation.               */
/*     -rR : R = distance between grid points in meters.                 */
/*     -tT : T = timestep in seconds.                                    */
/*     -s  : Print timing statistics.                                    */
/*     -o  : Print out relaxation residual values.                       */
/*     -h  : Print out command line options.                             */
/*                                                                       */
/*  Default: OCEAN -n130 -p1 -e1e-7 -r20000.0 -t28800.0                  */
/*                                                                       */
/*  NOTE: This code works under both the FORK and SPROC models.          */
/*                                                                       */
/*************************************************************************/



#ifndef _PARMACS_COMMON_INCLUDES_
#define _PARMACS_COMMON_INCLUDES_

#include <pthread.h>
#include <parmacs_config.h>
#include <parmacs_types.h>
#include <parmacs_decl.h>
#include <parmacs.h>

#endif /* _PARMACS_COMMON_INCLUDES_ */



#define DEFAULT_N      258
#define DEFAULT_P        1
#define DEFAULT_E        1e-7
#define DEFAULT_T    28800.0
#define DEFAULT_R    20000.0
#define UP               0
#define DOWN             1
#define LEFT             2
#define RIGHT            3
#define UPLEFT           4
#define UPRIGHT          5
#define DOWNLEFT         6
#define DOWNRIGHT        7
#define PAGE_SIZE     4096

#include <stdio.h>
#include <math.h>
#include <time.h>

struct multi_struct {
   double err_multi;
} *multi;

struct global_struct {
   int id;
   int starttime;
   int trackstart;
   double psiai;
   double psibi;
} *global;

double ****psi;
double ****psim;
double ***psium;
double ***psilm;
double ***psib;
double ***ga;
double ***gb;
double ****work1;
double ***work2;
double ***work3;
double ****work4;
double ****work5;
double ***work6;
double ****work7;
double ****temparray;
double ***tauz;
double ***oldga;
double ***oldgb;
double *f;
double ****q_multi;
double ****rhs_multi;

struct locks_struct {
   parmacs_lock_t (idlock);
   parmacs_lock_t (psiailock);
   parmacs_lock_t (psibilock);
   parmacs_lock_t (donelock);
   parmacs_lock_t (error_lock);
   parmacs_lock_t (bar_lock);
} *locks;

struct bars_struct {
   parmacs_barrier_t (iteration);
   parmacs_barrier_t (gsudn);
   parmacs_barrier_t (p_setup); 
   parmacs_barrier_t (p_redph); 
   parmacs_barrier_t (p_soln); 
   parmacs_barrier_t (p_subph); 
   parmacs_barrier_t (sl_prini);
   parmacs_barrier_t (sl_psini);
   parmacs_barrier_t (sl_onetime);
   parmacs_barrier_t (sl_phase_1);
   parmacs_barrier_t (sl_phase_2);
   parmacs_barrier_t (sl_phase_3);
   parmacs_barrier_t (sl_phase_4);
   parmacs_barrier_t (sl_phase_5);
   parmacs_barrier_t (sl_phase_6);
   parmacs_barrier_t (sl_phase_7);
   parmacs_barrier_t (sl_phase_8);
   parmacs_barrier_t (sl_phase_9);
   parmacs_barrier_t (sl_phase_10);
   parmacs_barrier_t (error_barrier);
} *bars;

void subblock();
void slave();
int log_2(int);
void printerr(char *);

int nprocs = DEFAULT_P;
double h1 = 1000.0;
double h3 = 4000.0;
double h = 5000.0;
double lf = -5.12e11;
double res = DEFAULT_R;
double dtau = DEFAULT_T;
double f0 = 8.3e-5;
double beta = 2.0e-11;
double gpr = 0.02;
int im = DEFAULT_N;
int jm;
double tolerance = DEFAULT_E;
double eig2;
double ysca;
int jmm1;
double pi;
double t0 = 0.5e-4 ;
double outday0 = 1.0;
double outday1 = 2.0;
double outday2 = 2.0;
double outday3 = 2.0;
double factjacob;
double factlap;
int numlev;
int *imx;
int *jmx;
double *lev_res;
double *lev_tol;
double maxwork = 10000.0;

struct Global_Private {
  char pad[PAGE_SIZE];
  int *rel_num_x;
  int *rel_num_y;
  int *eist;     
  int *ejst;     
  int *oist;     
  int *ojst;     
  int *rlist;    
  int *rljst;    
  int *rlien;    
  int *rljen;    
  int rownum;
  int colnum;
  int neighbors[8];
  double multi_time;
  double total_time;
} *gp;

double *i_int_coeff;
double *j_int_coeff;
int xprocs;
int yprocs;
int *xpts_per_proc;
int *ypts_per_proc;
int minlevel;
int do_stats = 0;
int do_output = 0;

void main(argc, argv)

int argc;
char *argv[];

{
   int i;
   int j;
   int k;
   double work_multi;
   int my_num;
   int x_part;
   int y_part;
   int d_size;
   int itemp;
   int jtemp;
   double procsqrt;
   FILE *fileptr;
   int iindex;
   int temp = 0;
   char c;
   double min_total;
   double max_total;
   double avg_total;
   double min_multi;
   double max_multi;
   double avg_multi;
   double min_frac;
   double max_frac;
   double avg_frac;
   int ch;
   extern char *optarg;
   unsigned int computeend;
   unsigned int start;

   {(start) = _parmacs_clock();}

   while ((ch = getopt(argc, argv, "n:p:e:r:t:soh")) != -1) {
     switch(ch) {
     case 'n': im = atoi(optarg);
               if (log_2(im-2) == -1) {
                 printerr("Grid must be ((power of 2)+2) in each dimension\n");
                 exit(-1);
               }
               break;
     case 'p': nprocs = atoi(optarg);
               if (nprocs < 1) {
                 printerr("P must be >= 1\n");
                 exit(-1);
               }
               if (log_2(nprocs) == -1) {
                 printerr("P must be a power of 2\n");
                 exit(-1);
               }
               break;
     case 'e': tolerance = atof(optarg); break;
     case 'r': res = atof(optarg); break;
     case 't': dtau = atof(optarg); break;
     case 's': do_stats = !do_stats; break;
     case 'o': do_output = !do_output; break;
     case 'h': printf("Usage: OCEAN <options>\n\n");
               printf("options:\n");
               printf("  -nN : Simulate NxN ocean.  N must be (power of 2)+2.\n");
               printf("  -pP : P = number of processors.  P must be power of 2.\n");
               printf("  -eE : E = error tolerance for iterative relaxation.\n");
               printf("  -rR : R = distance between grid points in meters.\n");
               printf("  -tT : T = timestep in seconds.\n");
               printf("  -s  : Print timing statistics.\n");
               printf("  -o  : Print out relaxation residual values.\n");
               printf("  -h  : Print out command line options.\n\n");
               printf("Default: OCEAN -n%1d -p%1d -e%1g -r%1g -t%1g\n",
                       DEFAULT_N,DEFAULT_P,DEFAULT_E,DEFAULT_R,DEFAULT_T);
               exit(0);
               break;
     }
   }

   {
    register unsigned long int _parmacs_m4_memory_;
    
    _parmacs_m4_memory_ = (unsigned long int)(
      (60000000)
    );
    _parmacs_main_initenv(_parmacs_m4_memory_);
} 

   jm = im;
   printf("\n");
   printf("Ocean simulation with W-cycle multigrid solver\n");
   printf("    Processors                         : %1d\n",nprocs);
   printf("    Grid size                          : %1d x %1d\n",im,jm);
   printf("    Grid resolution (meters)           : %0.2f\n",res);
   printf("    Time between relaxations (seconds) : %0.0f\n",dtau);
   printf("    Error tolerance                    : %0.7g\n",tolerance);
   printf("\n");

   xprocs = 0;
   yprocs = 0;
   procsqrt = sqrt((double) nprocs);
   j = (int) procsqrt;
   while ((xprocs == 0) && (j > 0)) {
     k = nprocs / j;
     if (k * j == nprocs) {
       if (k > j) {
         xprocs = j;
         yprocs = k;
       } else {
         xprocs = k;
         yprocs = j;
       }
     }
     j--;
   }
   if (xprocs == 0) {
     printerr("Could not find factors for subblocking\n");
     exit(-1);
   }  

   minlevel = 0;
   itemp = 1;
   jtemp = 1;
   numlev = 0;
   minlevel = 0;
   while (itemp < (im-2)) {
     itemp = itemp*2;
     jtemp = jtemp*2;
     if ((itemp/yprocs > 1) && (jtemp/xprocs > 1)) {
       numlev++;
     }
   }  
   
   if (numlev == 0) {
     printerr("Must have at least 2 grid points per processor in each dimension\n");
     exit(-1);
   }

   imx = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
   jmx = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
   lev_res = (double *) _parmacs_memory_alloc((int)(numlev*sizeof(double)));;
   lev_tol = (double *) _parmacs_memory_alloc((int)(numlev*sizeof(double)));;
   i_int_coeff = (double *) _parmacs_memory_alloc((int)(numlev*sizeof(double)));;
   j_int_coeff = (double *) _parmacs_memory_alloc((int)(numlev*sizeof(double)));;
   xpts_per_proc = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
   ypts_per_proc = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;

   imx[numlev-1] = im;
   jmx[numlev-1] = jm;
   lev_res[numlev-1] = res;
   lev_tol[numlev-1] = tolerance;

   for (i=numlev-2;i>=0;i--) {
     imx[i] = ((imx[i+1] - 2) / 2) + 2;
     jmx[i] = ((jmx[i+1] - 2) / 2) + 2;
     lev_res[i] = lev_res[i+1] * 2;
   }

   for (i=0;i<numlev;i++) {
     xpts_per_proc[i] = (jmx[i]-2) / xprocs;
     ypts_per_proc[i] = (imx[i]-2) / yprocs;
   }  
   for (i=numlev-1;i>=0;i--) {
     if ((xpts_per_proc[i] < 2) || (ypts_per_proc[i] < 2)) {
       minlevel = i+1;
       break;
     }
   }    
 
   for (i=0;i<numlev;i++) {
     temp += imx[i];
   }
   temp = 0;
   j = 0;
   for (k=0;k<numlev;k++) {
     for (i=0;i<imx[k];i++) {
       j++;
       temp += jmx[k];
     }
   }

   d_size = nprocs*sizeof(double ***);
   psi = (double ****) _parmacs_memory_alloc((int)(d_size));;
   psim = (double ****) _parmacs_memory_alloc((int)(d_size));;
   work1 = (double ****) _parmacs_memory_alloc((int)(d_size));;
   work4 = (double ****) _parmacs_memory_alloc((int)(d_size));;
   work5 = (double ****) _parmacs_memory_alloc((int)(d_size));;
   work7 = (double ****) _parmacs_memory_alloc((int)(d_size));;
   temparray = (double ****) _parmacs_memory_alloc((int)(d_size));;

   d_size = 2*sizeof(double **);
   for (i=0;i<nprocs;i++) {
     psi[i] = (double ***) _parmacs_memory_alloc((int)(d_size));;
     psim[i] = (double ***) _parmacs_memory_alloc((int)(d_size));;
     work1[i] = (double ***) _parmacs_memory_alloc((int)(d_size));;
     work4[i] = (double ***) _parmacs_memory_alloc((int)(d_size));;
     work5[i] = (double ***) _parmacs_memory_alloc((int)(d_size));;
     work7[i] = (double ***) _parmacs_memory_alloc((int)(d_size));;
     temparray[i] = (double ***) _parmacs_memory_alloc((int)(d_size));;
   }

   d_size = nprocs*sizeof(double **);
   psium = (double ***) _parmacs_memory_alloc((int)(d_size));;
   psilm = (double ***) _parmacs_memory_alloc((int)(d_size));;
   psib = (double ***) _parmacs_memory_alloc((int)(d_size));;
   ga = (double ***) _parmacs_memory_alloc((int)(d_size));;
   gb = (double ***) _parmacs_memory_alloc((int)(d_size));;
   work2 = (double ***) _parmacs_memory_alloc((int)(d_size));;
   work3 = (double ***) _parmacs_memory_alloc((int)(d_size));;
   work6 = (double ***) _parmacs_memory_alloc((int)(d_size));;
   tauz = (double ***) _parmacs_memory_alloc((int)(d_size));;
   oldga = (double ***) _parmacs_memory_alloc((int)(d_size));;
   oldgb = (double ***) _parmacs_memory_alloc((int)(d_size));;

   gp = (struct Global_Private *) _parmacs_memory_alloc((int)((nprocs+1)*sizeof(struct Global_Private)));;
   for (i=0;i<nprocs;i++) {
     gp[i].rel_num_x = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].rel_num_y = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].eist = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].ejst = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].oist = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].ojst = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].rlist = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].rljst = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].rlien = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].rljen = (int *) _parmacs_memory_alloc((int)(numlev*sizeof(int)));;
     gp[i].multi_time = 0;
     gp[i].total_time = 0;
   }

   subblock();

   x_part = (jm - 2)/xprocs + 2;
   y_part = (im - 2)/yprocs + 2;

   d_size = x_part*y_part*sizeof(double) + y_part*sizeof(double *);

   global = (struct global_struct *) _parmacs_memory_alloc((int)(sizeof(struct global_struct)));;  
   for (i=0;i<nprocs;i++) {
     psi[i][0] = (double **) _parmacs_memory_alloc((int)(d_size));;
     psi[i][1] = (double **) _parmacs_memory_alloc((int)(d_size));;
     psim[i][0] = (double **) _parmacs_memory_alloc((int)(d_size));;
     psim[i][1] = (double **) _parmacs_memory_alloc((int)(d_size));;
     psium[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     psilm[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     psib[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     ga[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     gb[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work1[i][0] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work1[i][1] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work2[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work3[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work4[i][0] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work4[i][1] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work5[i][0] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work5[i][1] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work6[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work7[i][0] = (double **) _parmacs_memory_alloc((int)(d_size));;
     work7[i][1] = (double **) _parmacs_memory_alloc((int)(d_size));;
     temparray[i][0] = (double **) _parmacs_memory_alloc((int)(d_size));;
     temparray[i][1] = (double **) _parmacs_memory_alloc((int)(d_size));;
     tauz[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     oldga[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
     oldgb[i] = (double **) _parmacs_memory_alloc((int)(d_size));;
   }
   f = (double *) _parmacs_memory_alloc((int)(im*sizeof(double)));;

   multi = (struct multi_struct *) _parmacs_memory_alloc((int)(sizeof(struct multi_struct)));;

   d_size = numlev*sizeof(double **);
   if (numlev%2 == 1) {         /* To make sure that the actual data 
                                   starts double word aligned, add an extra
                                   pointer */
     d_size += sizeof(double **);
   }
   for (i=0;i<numlev;i++) {
     d_size += ((imx[i]-2)/yprocs+2)*((jmx[i]-2)/xprocs+2)*sizeof(double)+
              ((imx[i]-2)/yprocs+2)*sizeof(double *);
   }

   d_size *= nprocs;

   if (nprocs%2 == 1) {         /* To make sure that the actual data 
                                   starts double word aligned, add an extra
                                   pointer */
     d_size += sizeof(double ***);
   }

   d_size += nprocs*sizeof(double ***);
   q_multi = (double ****) _parmacs_memory_alloc((int)(d_size));;
   rhs_multi = (double ****) _parmacs_memory_alloc((int)(d_size));;

   locks = (struct locks_struct *) _parmacs_memory_alloc((int)(sizeof(struct locks_struct)));;
   bars = (struct bars_struct *) _parmacs_memory_alloc((int)(sizeof(struct bars_struct)));;

   {
    _parmacs_lockinit(&(locks->idlock));
}
   {
    _parmacs_lockinit(&(locks->psiailock));
}
   {
    _parmacs_lockinit(&(locks->psibilock));
}
   {
    _parmacs_lockinit(&(locks->donelock));
}
   {
    _parmacs_lockinit(&(locks->error_lock));
}
   {
    _parmacs_lockinit(&(locks->bar_lock));
}

   {
    _parmacs_barinit(&(bars->iteration));
}
   {
    _parmacs_barinit(&(bars->gsudn));
}
   {
    _parmacs_barinit(&(bars->p_setup));
} 
   {
    _parmacs_barinit(&(bars->p_redph));
} 
   {
    _parmacs_barinit(&(bars->p_soln));
} 
   {
    _parmacs_barinit(&(bars->p_subph));
} 
   {
    _parmacs_barinit(&(bars->sl_prini));
}
   {
    _parmacs_barinit(&(bars->sl_psini));
}
   {
    _parmacs_barinit(&(bars->sl_onetime));
}
   {
    _parmacs_barinit(&(bars->sl_phase_1));
}
   {
    _parmacs_barinit(&(bars->sl_phase_2));
}
   {
    _parmacs_barinit(&(bars->sl_phase_3));
}
   {
    _parmacs_barinit(&(bars->sl_phase_4));
}
   {
    _parmacs_barinit(&(bars->sl_phase_5));
}
   {
    _parmacs_barinit(&(bars->sl_phase_6));
}
   {
    _parmacs_barinit(&(bars->sl_phase_7));
}
   {
    _parmacs_barinit(&(bars->sl_phase_8));
}
   {
    _parmacs_barinit(&(bars->sl_phase_9));
}
   {
    _parmacs_barinit(&(bars->sl_phase_10));
}
   {
    _parmacs_barinit(&(bars->error_barrier));
}

   link_all();

   multi->err_multi = 0.0;
   i_int_coeff[0] = 0.0;
   j_int_coeff[0] = 0.0;
   for (i=0;i<numlev;i++) {
     i_int_coeff[i] = 1.0/(imx[i]-1);
     j_int_coeff[i] = 1.0/(jmx[i]-1);
   }

/* initialize constants and variables

   id is a global shared variable that has fetch-and-add operations
   performed on it by processes to obtain their pids.   */

   global->id = 0;
   global->psibi = 0.0;
   pi = atan(1.0);
   pi = 4.*pi;

   factjacob = -1./(12.*res*res);
   factlap = 1./(res*res);
   eig2 = -h*f0*f0/(h1*h3*gpr);

   jmm1 = jm-1 ;
   ysca = ((double) jmm1)*res ;

   im = (imx[numlev-1]-2)/yprocs + 2;
   jm = (jmx[numlev-1]-2)/xprocs + 2;

   for (i=1;i<nprocs;i++) {
     {
    _parmacs_create(slave);
}  
   }

   if (do_output) {
     printf("                       MULTIGRID OUTPUTS\n");
   }

   slave();
   {
    _parmacs_wait_for_end(nprocs-1);
}
   {(computeend) = _parmacs_clock();}

   printf("\n");
   printf("                       PROCESS STATISTICS\n");
   printf("                  Total          Multigrid         Multigrid\n");
   printf(" Proc             Time             Time            Fraction\n");
   printf("    0   %15.0f    %15.0f        %10.3f\n",
          gp[0].total_time,gp[0].multi_time,
          gp[0].multi_time/gp[0].total_time);

   if (do_stats) {
     min_total = max_total = avg_total = gp[0].total_time;
     min_multi = max_multi = avg_multi = gp[0].multi_time;
     min_frac = max_frac = avg_frac = gp[0].multi_time/gp[0].total_time;
     for (i=1;i<nprocs;i++) {
       if (gp[i].total_time > max_total) {
         max_total = gp[i].total_time;
       }
       if (gp[i].total_time < min_total) {
         min_total = gp[i].total_time;
       }
       if (gp[i].multi_time > max_multi) {
         max_multi = gp[i].multi_time;
       }
       if (gp[i].multi_time < min_multi) {
         min_multi = gp[i].multi_time;
       }
       if (gp[i].multi_time/gp[i].total_time > max_frac) {
         max_frac = gp[i].multi_time/gp[i].total_time;
       }
       if (gp[i].multi_time/gp[i].total_time < min_frac) {
         min_frac = gp[i].multi_time/gp[i].total_time;
       }
       avg_total += gp[i].total_time;
       avg_multi += gp[i].multi_time;
       avg_frac += gp[i].multi_time/gp[i].total_time;
     }
     avg_total = avg_total / nprocs;
     avg_multi = avg_multi / nprocs;
     avg_frac = avg_frac / nprocs;
     for (i=1;i<nprocs;i++) {
       printf("  %3d   %15.0f    %15.0f        %10.3f\n",
	      i,gp[i].total_time,gp[i].multi_time,
	      gp[i].multi_time/gp[i].total_time);
     }
     printf("  Avg   %15.0f    %15.0f        %10.3f\n",
            avg_total,avg_multi,avg_frac);
     printf("  Min   %15.0f    %15.0f        %10.3f\n",
            min_total,min_multi,min_frac);
     printf("  Max   %15.0f    %15.0f        %10.3f\n",
            max_total,max_multi,max_frac);
   }
   printf("\n");

   global->starttime = start;
   printf("                       TIMING INFORMATION\n");
   printf("Start time                        : %16d\n",
           global->starttime);
   printf("Initialization finish time        : %16d\n",
           global->trackstart);
   printf("Overall finish time               : %16d\n",
           computeend);
   printf("Total time with initialization    : %16d\n",
           computeend-global->starttime);
   printf("Total time without initialization : %16d\n",
           computeend-global->trackstart);
   printf("    (excludes first timestep)\n");
   printf("\n");

   {
    _parmacs_main_end();
}
}

int log_2(number)

int number;

{
  int cumulative = 1;
  int out = 0;
  int done = 0;

  while ((cumulative < number) && (!done) && (out < 50)) {
    if (cumulative == number) {
      done = 1;
    } else {
      cumulative = cumulative * 2;
      out ++;
    }
  }

  if (cumulative == number) {
    return(out);
  } else {
    return(-1);
  }
}

void printerr(s)

char *s;

{
  fprintf(stderr,"ERROR: %s\n",s);
}

