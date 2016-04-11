
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

#define MASTER            0
#define RED_ITER          0
#define BLACK_ITER        1
#define UP                0
#define DOWN              1
#define LEFT              2
#define RIGHT             3
#define UPLEFT            4
#define UPRIGHT           5
#define DOWNLEFT          6
#define DOWNRIGHT         7
#define PAGE_SIZE      4096

extern struct multi_struct {
   double err_multi;
} *multi;

extern struct global_struct {
   int id;
   int starttime;
   int trackstart;
   double psiai;
   double psibi;
} *global;

extern double eig2;
extern double ysca;
extern int jmm1;
extern double pi;
extern double t0;

extern double ****psi;
extern double ****psim;
extern double ***psium;
extern double ***psilm;
extern double ***psib;
extern double ***ga;
extern double ***gb;
extern double ****work1;
extern double ***work2;
extern double ***work3;
extern double ****work4;
extern double ****work5;
extern double ***work6;
extern double ****work7;
extern double ****temparray;
extern double ***tauz;
extern double ***oldga;
extern double ***oldgb;
extern double *f;
extern double ****q_multi;
extern double ****rhs_multi;

extern struct locks_struct {
   parmacs_lock_t (idlock);
   parmacs_lock_t (psiailock);
   parmacs_lock_t (psibilock);
   parmacs_lock_t (donelock);
   parmacs_lock_t (error_lock);
   parmacs_lock_t (bar_lock);
} *locks;
 
extern struct bars_struct {
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

extern double factjacob;
extern double factlap;

extern struct Global_Private {
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

extern double *i_int_coeff;
extern double *j_int_coeff;
extern int xprocs;
extern int yprocs;

extern int numlev;
extern int *imx;
extern int *jmx;
extern double *lev_res;
extern double *lev_tol;
extern double maxwork;
extern int *xpts_per_proc;
extern int *ypts_per_proc;
extern int minlevel;
extern double outday0;
extern double outday1;
extern double outday2;
extern double outday3;

extern int nprocs;
extern double h1;
extern double h3;
extern double h;
extern double lf;
extern double res;
extern double dtau;
extern double f0;
extern double beta;
extern double gpr;
extern int im;
extern int jm;
extern int do_stats;
extern int do_output;
extern int *multi_times;
extern int *total_times;
