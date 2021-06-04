
#include "paul.h"

static struct domain theDomain = {0};

void start_clock( struct domain * );
void read_par_file( struct domain * );
void setupGrid( struct domain * );
void setupDomain( struct domain * );
void setupCells( struct domain * );

void start() {
   start_clock( &theDomain ); 
   read_par_file( &theDomain ); 
   setupGrid( &theDomain );   
   setupDomain( &theDomain );
   setupCells( &theDomain );
   FILE * rFile = fopen("LOGS/report.dat","w");
   fclose(rFile);
}

void set_wcell( struct domain * );
double getmindt( struct domain * );
void check_dt( struct domain * , double * );
void possiblyOutput( struct domain * , int );
void timestep( struct domain * , double );

void steps(int max_numsteps) {
   int numsteps = 0;
   while( !(theDomain.final_step) ){
      set_wcell( &theDomain );
      double dt = getmindt( &theDomain );
      check_dt( &theDomain , &dt );
      possiblyOutput( &theDomain , 0 );
      timestep( &theDomain , dt );
      numsteps += 1;
      if (numsteps >= max_numsteps) return;
   }
}

void generate_log( struct domain * );
void freeDomain( struct domain * );

void finish() {
   possiblyOutput( &theDomain , 1 );
   generate_log( &theDomain );
   freeDomain( &theDomain );
}

int copy_data(      
   int selector, int lrdat, double rdat[]) 
{
   int nz = theDomain.nz;
   if (nz > lrdat) return(-1);
   struct cell * theCells = theDomain.theCells;
   int i;
   for( i=0 ; i<nz ; ++i ){
      struct cell * c = theCells+i;
      double x;
      if (selector == 0) { x = c->riph; } // r at outer face
      else if (selector == 1) { x = c->miph; } // m at outer face
      else if (selector == 2) { x = c->wiph; } // v at outer face
      else if (selector == 3) { x = c->pot; } // grav potential at outer face
      else if (selector == 4) { x = c->dm; } // dm from outer to inner face
      else if (selector == 5) { x = c->dr; } // dr from outer to inner face
      else if (selector == 6) { x = c->prim[RHO]; } // cell density
      else if (selector == 7) { x = c->prim[PPP]; } // cell pressure
      else if (selector == 8) { x = c->prim[VRR]; } // cell velocity
      else if (selector == 9) { x = c->prim[XXX]; } // cell mass fraction of species X
      else if (selector == 10) { x = c->prim[AAA]; } // cell RTI alpha
      else if (selector == 11) { x = c->cons[DDD]; } // cell mass
      else if (selector == 11) { x = c->cons[SRR]; } // cell momentum
      else if (selector == 11) { x = c->cons[TAU]; } // cell total energy (KE + IE + PE)
      else if (selector == 11) { x = c->cons[XXX]; } // cell mass of species X
      else { return(-1); }
      int k = nz - 1 - i;
      rdat[k] = x;
   }
   return (0);
}

int dispatch(      
   int *iop, 
   int *lrpar, double rpar[], 
   int *lipar, int ipar[],
   int *lrdat, double rdat[], 
   int *lidat, int idat[]) {
   int op = *iop;
   if (op == 0) {
      start();
   } else if (op == 1) {
      steps(ipar[0]);
      rdat[0] = theDomain.t;
      idat[0] = theDomain.final_step;
      idat[1] = theDomain.count_steps;
   } else if (op == 2) {
      finish();
   } else if (op == 3) {
      return copy_data(ipar[0], *lrdat, rdat);
   } else {
      return -1;
   }
   return 0;
}


