
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "paul.h"

void setupGrid( struct domain * theDomain ){

   int Num_R = theDomain->theParList.Num_R;
   int LogZoning = theDomain->theParList.LogZoning;

   double Rmin = theDomain->theParList.rmin;
   double Rmax = theDomain->theParList.rmax;

   int nz = Num_R;

   theDomain->nz = nz;
   theDomain->theCells = (struct cell *) malloc( nz*sizeof(struct cell));

   int i;

   double dx = 1./(double)Num_R;
   double x0 = 0.;
   double R0 = theDomain->theParList.LogRadius;
   for( i=0 ; i<nz ; ++i ){
      double xm = x0 + ((double)i   )*dx;
      double xp = x0 + ((double)i+1.)*dx;
      double rp,rm;
      if( LogZoning == 0 ){
         rp = Rmin + xp*(Rmax-Rmin);
         rm = Rmin + xm*(Rmax-Rmin);
      }else if( LogZoning == 1 ){
         rp = Rmin*pow(Rmax/Rmin,xp);
         rm = Rmin*pow(Rmax/Rmin,xm);
      }else{
         rp = Rmax*(pow(Rmax/R0,xp)-1.)/(Rmax/R0-1.) + Rmin;
         rm = Rmax*(pow(Rmax/R0,xm)-1.)/(Rmax/R0-1.) + Rmin;
      }
      theDomain->theCells[i].riph = rp;
      theDomain->theCells[i].dr   = rp - rm;
   }

}


