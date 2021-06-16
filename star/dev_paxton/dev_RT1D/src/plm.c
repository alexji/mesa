
#include "paul.h"

double minmod( double a , double b , double c ){
   double m = a;
   if( a*b < 0.0 ) m = 0.0;
   if( fabs(b) < fabs(m) ) m = b;
   if( b*c < 0.0 ) m = 0.0;
   if( fabs(c) < fabs(m) ) m = c;
   return(m);
}

double get_g( struct cell * );
void calculate_pot( struct domain * );
void calculate_mass( struct domain * );

void plm( struct domain * theDomain , int PPP_flag ){ // piece-wise linear minmod

   if ( PPP_flag ) {
      int gE = theDomain->theParList.grav_e_mode;
      if( gE ) calculate_mass(theDomain);
      if( gE == 2 ){
         calculate_pot( theDomain );
      }
   }

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   double PLM = theDomain->theParList.PLM;
   int Bal = theDomain->theParList.grav_bal;

   int i,q;
   #pragma omp for private(i,q)
   for( i=0 ; i<nz ; ++i ){
      int im = i-1;
      int ip = i+1;
      if( i==0 ) im = 0;
      if( i==nz-1 ) ip = nz-1;
      struct cell * c  = theCells+i;
      struct cell * cL = theCells+im;
      struct cell * cR = theCells+ip;
      double drL = cL->dr;
      double drC = c->dr;
      double drR = cR->dr;
      for( q=0 ; q<NUM_Q ; ++q ){
         if ( (PPP_flag && q == PPP) || (!PPP_flag && q != PPP) ){
            double pL = cL->prim[q];
            double pC = c->prim[q];
            double pR = cR->prim[q];
         
            if( q==PPP && Bal ){//&& i>0 && i<nz-1){
               // factor out HSE P differences from reconstruction of face P.
               // Remove HSE solution from pL and pR.  
               // restore it in riemann_flux.
               double gHSE_C = c->prim[RHO]*get_g(c);
               double gHSE_L = cL->prim[RHO]*get_g(cL);
               double gHSE_R = cR->prim[RHO]*get_g(cR);

               double pHSE_L = pC - .5*gHSE_C*drC - .5*gHSE_L*drL;
               double pHSE_R = pC + .5*gHSE_C*drC + .5*gHSE_R*drR;

               pL -= pHSE_L;
               pR -= pHSE_R;
               pC = 0.0;
            }
            double sL = pC - pL;
            sL /= .5*( drC + drL );
            double sR = pR - pC;
            sR /= .5*( drR + drC );
            double sC = pR - pL;
            sC /= .5*( drL + drR ) + drC;
            c->grad[q] = minmod( PLM*sL , sC , PLM*sR ); // slope d/dr across cell
         }
      }
   }
}

