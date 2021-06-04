
#include "paul.h"

double get_moment_arm( double , double );
void initial( double * , double ); 
double get_g( struct cell * );
double get_GMr( struct cell * );
double get_dV( double , double );
void prim2cons( double * , double * , double , double );

void boundary( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   int gE = theDomain->theParList.grav_e_mode;
   
   struct cell * cB = theCells+nz-1;
   double rp = cB->riph;
   double rm = rp-cB->dr;
   double r = get_moment_arm(rp,rm);
   initial( cB->prim , r );
   double dV = get_dV( rp , rm );

   double GMr = 0.0;
   if( gE == 3 ) GMr = get_GMr( cB );
   prim2cons( cB->prim , cB->cons , GMr , dV ); 

}
