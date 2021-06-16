#include "paul.h"

void set_wcell( struct domain * );
double getmindt( struct domain * );
void check_dt( struct domain * , double * );
void possiblyOutput( struct domain * , int );

void radial_flux( struct domain * , double );
void add_source( struct domain * , double );
void move_cells( struct domain * , double , double );
void calc_dr( struct domain * );
void calc_prim( struct domain * );
void boundary( struct domain * );
void adjust_RK_cons( struct domain * , double );
void AMR( struct domain * ); 

void timestep( struct domain * theDomain ){







   set_wcell( theDomain );
   
   
   
   
   
   
   
   double dt = getmindt( theDomain );
   check_dt( theDomain , &dt );
   possiblyOutput( theDomain , 0 );
   
   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   int i;
   for( i=0 ; i<nz ; ++i ){
      struct cell * c = theCells+i;
      memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
   }

   //onestep( theDomain , 0.0 ,     dt , 1 , 0 );
   
      radial_flux( theDomain , dt );
      add_source( theDomain , dt );
      move_cells( theDomain , 0.0 , dt );      
      calc_dr( theDomain );
      calc_prim( theDomain );
      boundary( theDomain );
   
   adjust_RK_cons( theDomain , 0.5 );
   
   //onestep( theDomain , 0.5 , 0.5*dt , 0 , 1 );

      radial_flux( theDomain , 0.5*dt );
      add_source( theDomain , 0.5*dt );
      calc_dr( theDomain );
      calc_prim( theDomain );
      AMR( theDomain );
      boundary( theDomain );

   theDomain->t += dt;   
   theDomain->dt = dt;   
   theDomain->count_steps += 1;

}
