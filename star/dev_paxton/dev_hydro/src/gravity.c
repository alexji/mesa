
#include "paul.h"

static int grav_E_mode = 0;
static double grav_G = 0.0;

void setGravityParams( struct domain * theDomain ){
   grav_E_mode = theDomain->theParList.grav_e_mode;
   grav_G = theDomain->theParList.grav_G;
}

double get_dV( double , double );

double get_g( struct cell * c ){

   double G = grav_G;
   double rp = c->riph;
   double rm = c->riph - c->dr;


   double Mm = c->miph - c->dm;
   double rho = c->dm/(4./3.*M_PI*(rp*rp*rp-rm*rm*rm));

   double r2_3 = (rp*rp + rm*rm + rp*rm)/3.;
   double r3_4 = (rp*rp*rp + rp*rp*rm + rp*rm*rm + rm*rm*rm)/4.;

   double g = -G/r2_3*( Mm + 4./3.*M_PI*rho*(r3_4 - rm*rm*rm) );

   return( g );
 
}

double get_GMr( struct cell * c ){
   double G = grav_G;
   double rp = c->riph;
   double rm = rp - c->dr;
   double rc = .5*(rp+rm);
   double M = c->miph - c->dm;
   double dV = 4./3.*M_PI*(rp*rp*rp-rm*rm*rm);
   double rho = c->dm/dV;
   double r2_3 = (rp*rp + rm*rm + rp*rm)/3.;
   double r4_5 = (rp*rp*rp*rp + rp*rp*rp*rm + rp*rp*rm*rm + rp*rm*rm*rm + rm*rm*rm*rm )/5.;

   double mt = M - 4./3.*M_PI*rm*rm*rm*rho;

   double eps = G/r2_3*( mt*rc + 4./3.*M_PI*rho*r4_5 );

   return( eps );
}

void calculate_mass( struct domain * theDomain ){
   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   double M = theDomain->point_mass;
   int i;
   for( i=0 ; i<nz ; ++i ){
      struct cell * c = theCells+i;
      c->dm = c->cons[DDD];
      M += c->dm;
      c->miph = M;
   }
}

void calculate_pot( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   int i;
   double rmax = theCells[nz-1].riph;
   double M    = theCells[nz-1].miph;
   double pot  = -grav_G*M/rmax;

   for( i=nz-1 ; i>=0 ; --i ){
      struct cell * c = theCells+i;
      double dr = c->dr;
      double g = get_g( c );
      c->pot = pot + .5*g*dr; // at center of cell
      pot += g*dr; // at next face
   }

}

void grav_src( struct cell * c , double dt ){

   double rho = c->prim[RHO];
   double v   = c->prim[VRR];
   double f = get_g( c );
   double rp, rm;
   rp = c->riph;
   rm = rp - c->dr;
   double dV = get_dV( rp , rm );
   double dVdt = dV*dt;
   c->cons[SRR] += rho*f*dVdt;
   if( grav_E_mode == 0 ) c->cons[TAU] += rho*v*f*dVdt;

}


