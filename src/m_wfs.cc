/* Copyright (C) 2000 Patrick Eriksson <patrick@rss.chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */


/////////////////////////////////////////////////////////////////////////////
//
// This file contains functions to calculate weightings functions (WFs)
// for a number of variables.
//
/////////////////////////////////////////////////////////////////////////////


#include "arts.h"
#include "md.h"
#include "vecmat.h"
#include "messages.h"          
#include "wsv.h"          
#include "math_funcs.h"          
#include "atm_funcs.h"          
#include "absorption.h"          

extern const Numeric PLANCK_CONST;
extern const Numeric BOLTZMAN_CONST;



//==========================================================================
//=== Function to join WF matrices and related variables
//==========================================================================

// k_join
//
// Patrick Eriksson 14.06.00

void k_join (
                    MATRIX&          ktot,
              const MATRIX&          knew,
                    ARRAYofstring&   ktot_names,
              const ARRAYofstring&   knew_names,
                    MATRIX&          ktot_index,
              const MATRIX&          knew_index,
                    MATRIX&          ktot_aux, 
              const MATRIX&          knew_aux )
{
  // Size of KTOT
  const size_t  ny1 = ktot.dim(1);  // length of measurement vector (y)
  const size_t  nx1 = ktot.dim(2);  // length of state vector (x)

  // If KTOT is empty, just copy the new data
  if ( (ny1==0) && (nx1==0) )
  {
    ktot       = knew;
    ktot_names = knew_names;
    ktot_index = knew_index;
    ktot_aux   = knew_aux;
  }

  // If KTOT not empty, append data
  else
  {
    // More sizes
    const size_t  ny2  = knew.dim(1);
    const size_t  nx2  = knew.dim(2);
    const size_t  nri1 = ktot_names.dim();  // number of retrieval identities
    const size_t  nri2 = knew_names.dim();

    // Indices used   
    size_t  iy, ix, ix2, iri, iri2;

    // Check that sizes match
    if ( ny1 != ny2 )
      throw runtime_error("The two WF matrices have different number of rows."); 

    // Make temporary storage for old data and reallocate the KTOT vars.
    const MATRIX           ktemp       = ktot;
    const MATRIX           ktemp_index = ktot_index;
    const MATRIX           ktemp_aux   = ktot_aux;
    const ARRAYofstring    ktemp_names = ktot_names;;
    ktot.newsize(ny1,nx1+nx2);
    ktot_index.newsize(nri1+nri2,2);
    ktot_aux.newsize(nx1+nx2,3);
    ktot_names.newsize(nri1+nri2);

    // Put in old KTOT "to the left" and KNEW "to the right"
    for ( ix=1; ix<=nx1; ix++ )
    {
      for ( iy=1; iy<=ny1; iy++ )
        ktot(iy,ix)  = ktemp(iy,ix);
      ktot_aux(ix,1) = ktemp_aux(ix,1);
      ktot_aux(ix,2) = ktemp_aux(ix,2);
      ktot_aux(ix,3) = ktemp_aux(ix,3);
    }
    for ( ix=1; ix<=nx2; ix++ )
    {
      ix2 = nx1 + ix;
      for ( iy=1; iy<=ny1; iy++ )
        ktot(iy,ix2)  = knew(iy,ix);
      ktot_aux(ix2,1) = knew_aux(ix,1);
      ktot_aux(ix2,2) = knew_aux(ix,2);
      ktot_aux(ix2,3) = knew_aux(ix,3);
    }  
    for ( iri=1; iri<=nri1; iri++ )
    {
      ktot_names(iri)   = ktemp_names(iri);
      ktot_index(iri,1) = ktemp_index(iri,1);
      ktot_index(iri,2) = ktemp_index(iri,2);
    } 
    for ( iri=1; iri<=nri2; iri++ )
    {
      iri2 = nri1 + iri;
      ktot_names(iri2)   = knew_names(iri);
      ktot_index(iri2,1) = knew_index(iri,1) + ktemp_index(nri1,2);
      ktot_index(iri2,2) = knew_index(iri,2) + ktemp_index(nri1,2);
    } 
  }
}



//==========================================================================
//=== Help functions for grid conversions
//==========================================================================

// GRID2GRID_INDEX
//
// This function gives the indeces of the grid points affected by a change at 
// a set of retrieval points. That is, a first step to calculate dkappa/dkp
// (see Equation 4 of the AUG section on "Basic atmospheric WFs). The second
// stap is performed by the function GRID2GRID_WEIGHTS found below.
//
// This function can be used to convert LOS WFs to WFs for vertical altitudes.
// 
// The profile/function for the retrieved quantity is treated to be piecewise
// linear. With other words, the base functions are tenth functions. The 
// functions areis treated to be 0 outside the end points. 
// The function returns a 2 column matrix, where column 1 is the lowest
// index affected and column 2 the highest index. If no points are affected
// 0's are returned. 
//
// An example:
// x0 = [1 2 3 4 5 6 7] and xp = [0 1 2 3 4.5 5 5.5 6 8] gives
// Is = [ 0     0
//        1     1
//        2     2
//        3     4
//        4     4
//        5     5
//        0     0
//        6     7
//        7     7 ]
// Note that the XP point 5.5 gives 0 as this point not embraces any point 
// of X0 as there are retrieval points at the neighbouring X0 points.
// 
// X0 is the original grid. For example, the LOS grid. 
// XP is the retrieval grid. For example, vertical altitudes
//
// Both vectors most be sorted with increasing values. 
// The length of XP must be > 1.
//
// The algorithm used here is not the most effective. However, more elaborated
// algorithms were tested but they failed for some special case, and  instead 
// a simplistic approach was used.
// 
//
// Patrick Eriksson 14.06.00

void grid2grid_index (
                    MATRIX&   is,
              const VECTOR&   x0,
              const VECTOR&   xp )
{
  const size_t n0 = x0.dim();        // length of original grid
  const size_t np = xp.dim();        // length if retrieval grid
        size_t i0, ip;               // counter for each grid

  // Resize is and set all values to 0
  is.newsize(np,2);
  is = 0.0;
 
  // Some checks
  if ( np < 2 )
    throw logic_error("The retrieval grid must have a length > 1.");
  for ( i0=2; i0<=n0; i0++ )
  {
    if ( x0(i0-1) >= x0(i0) )
      throw logic_error("The grids must be sorted with increasing values.");
  }
  for ( ip=2; ip<=np; ip++ )
  {
    if ( xp(ip-1) >= xp(ip) )
      throw logic_error("The grids must be sorted with increasing values.");
  }

  // Do something only if at least one point of X0 is inside XP
  if ( !( (x0(1)>xp(np)) || (x0(n0)<xp(1)) ) )
  {

    // First point of XP
    i0 = 1;
    for ( ; x0(i0) < xp(1); i0++ ) {}
    if ( x0(i0) < xp(2) )
    {
      is(1,1) = (Numeric) i0;
      for ( ; (i0<=n0) && (x0(i0)<xp(2)); i0++ ) {}
      is(1,2) = (Numeric) (i0 - 1);
    }

    // Points inside XP
    for ( ip=2; ip<np; ip++ )
    {
      i0 = 1;
      for ( ; (i0<=n0) && (x0(i0)<=xp(ip-1)); i0++ ) {}
      if ( (i0<=n0) && (x0(i0)<xp(ip+1)) )
      {
        is(ip,1) = (Numeric) i0;
        for ( ; (i0<=n0) && (x0(i0)<xp(ip+1)); i0++ ) {}
        is(ip,2) = (Numeric) (i0 - 1);
      }
    }

    // Last point of XP
    i0 = 1;
    for ( ; (i0<=n0) && (x0(i0)<=xp(np-1)); i0++ ) {}
    if ( (i0<=n0) && (x0(i0)<=xp(np)) )
    {
      is(np,1) = (Numeric) i0;
      for ( ; (i0<=n0) && (x0(i0)<=xp(np)); i0++ ) {}
      is(np,2) = (Numeric) (i0 - 1);
    }
  }
}



// GRID2GRID_INDEX
//
// This function returns the derivative dkappa/dkp (see Equation 4 of 
// the AUG section on "Basic atmospheric WFs"), here dennoted as the weight,
// for the indeces selected.
// This is an accompanying function to GRID2GRID_INDEX. This latter function
// gives the indeces for which the weight > 0.
//
// The weight is
// w = abs( (x0-xp(ip+-1)) / (xp(ip)-xp(ip+-1)));
// where x0 is the point of the original grid (kappa), xp(ip) is the 
// considered retrieval point and xp(ip+-1) is the retrieval point on the
// other side of x0.
//
// X0 is the original grid. For example, the LOS grid. 
// XP is the retrieval grid. For example, vertical altitudes
// I1 is the first index of X0 to consider
// I2 is the last index of X0 to consider
// IP is the index of the retrieval point
//
// Patrick Eriksson 14.06.00

void grid2grid_weights (
                    VECTOR&   w,
              const VECTOR&   x0,
              const size_t&   i1,
              const size_t&   i2,
              const VECTOR&   xp,
              const size_t&   ip )
{
  const size_t   np =xp.dim();        // number of retrieval points
  const size_t   nw = i2-i1+1;        // number of LOS points affected
        size_t   i;

  // Reallocate w
  w.newsize(nw);

  // First point of the retrieval grid
  if ( ip == 1 )
  {
    for ( i=0; i<nw; i++ )
      w(i+1) = ( xp(2)- x0(i1+i) ) / ( xp(2) - xp(1) );
  }
  // Points inside the retrieval grid
  else if ( ip < np )  
  {
    for ( i=0; i<nw; i++ )
    {
      if ( x0(i1+i) <= xp(ip) )
        w(i+1) = ( x0(i1+i) - xp(ip-1) ) / ( xp(ip) - xp(ip-1) );
      else
        w(i+1) = ( xp(ip+1) - x0(i1+i) ) / ( xp(ip+1) - xp(ip) );
    }
  }
  // Last point of the retrieval grid
  else
  {
    for ( i=0; i<nw; i++ )
      w(i+1) = ( x0(i1+i) - xp(np-1) ) / ( xp(np) - xp(np-1) );
  }
}



//==========================================================================
//=== Help functions for absloswfsCalc to calculate absorption LOS WFs
//==========================================================================

// ABSLOSWFS_1PASS
//
// This function covers cases where each point is passed once.
// That is 1D limb sounding and downward observatiuons are not covered.
//
// The expression used are described in sub-section 2.1 of the AUG section
// "Atmospheric WFs".
//
// Patrick Eriksson 09.06.00

void absloswfs_1pass (
                    MATRIX&   k,
                    VECTOR    y,
              const int&      start_index,
	      const int&      stop_index,   // this variable is used by 1D down
              const Numeric&  lstep,
              const MATRIX&   tr,
              const MATRIX&   s,
              const int&      ground,
              const VECTOR&   e_ground,
              const VECTOR&   y_ground )
{
  const size_t   nf = tr.dim(1);      // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t(nf,1.0);           // transmission to the sensor
           int   q;                   // LOS index (same notation as in AUG)

  if ( (ground==1) || (ground==start_index) )
    throw logic_error("The ground cannot be at one of the end points."); 

  // Resize K
  k.newsize( nf, start_index );

  // We start here at the LOS point closest to the sensor, that is,
  // reversed order compared to RTE_ITERATE  

  // The LOS point closest to the sensor 
  q  = stop_index;
  for ( iv=1; iv<=nf; iv++ )    
  {
    t(iv)   *= tr(iv,q);
    y(iv)    = (y(iv)-s(iv,q)*(1.0-tr(iv,q)))/tr(iv,q);
    k(iv,q)  = -lstep*(y(iv)-s(iv,q))*t(iv)/2;
  }

  // Points inside the LOS
  for ( q=stop_index+1; q<start_index; q++ )
  {
    // Not a ground point
    if ( q != ground )
    {
      for ( iv=1; iv<=nf; iv++ )    
      {
        y(iv)   = (y(iv)-s(iv,q)*(1.0-tr(iv,q)))/tr(iv,q);
        k(iv,q) = -lstep*( 2*(y(iv)-s(iv,q))*tr(iv,q) + s(iv,q) - s(iv,q-1) )
                                                                 *t(iv)/2;
        t(iv)  *= tr(iv,q); 
      }
    }
    // A ground point
    else
    {
      out1 <<
        "WARNING: The function absloswfs_1pass not tested for ground reflections\n";
      for ( iv=1; iv<=nf; iv++ )    
      {
        y(iv)   = ( y(iv) - e_ground(iv)*y_ground(iv) - 
                            s(iv,q)*(1.0-tr(iv,q))*(1-e_ground(iv)) ) / 
                  tr(iv,q) / (1-e_ground(iv));
        k(iv,q) = -lstep*( 2*(y(iv)-s(iv,q))*tr(iv,q)*(1-e_ground(iv)) + 
                           s(iv,q)*(1-e_ground(iv)) + 
                           e_ground(iv)*y_ground(iv) - s(iv,q-1) ) * t(iv) / 2;
        t(iv)  *= tr(iv,q)*(1-e_ground(iv)); 
      }
    }
  }

  // The LOS point furthest away from the sensor
  q = start_index;
  for ( iv=1; iv<=nf; iv++ )    
    k(iv,q)  = -lstep*(y(iv)-s(iv,q-1))*t(iv)/2;

  // To check the function: Y shall now be equal to the radiation entering 
  // the atmosphere and T the total transmission
}



// ABSLOSWFS_LIMB
//
// This function covers 1D limb sounding
//
// The expression used are described in sub-section 2.2 of the AUG section
// "Atmospheric WFs".
//
// Patrick Eriksson 14.06.00

void absloswfs_limb (
                    MATRIX&   k,
                    VECTOR    y,              // = y_q^q
                    VECTOR    yn,             // = y_space
              const int&      start_index,
              const Numeric&  lstep,
              const MATRIX&   tr,
              const MATRIX&   s,
              const int&      ground,
              const VECTOR&   e_ground )
{
  const size_t   nf = tr.dim(1);      // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t1q;                 // transmission tangent point - q
        VECTOR   tqn(nf,1);           // transmission q - sensor
        int      q;                   // LOS index (same notation as in AUG)
        Numeric  tv, tv1;             // transmission value for q and q-1

  // Calculate the square root of the total transmission
  bl( t1q, start_index, start_index, tr, ground, e_ground );
  t1q = sqrt(t1q);

  // Resize K
  k.newsize( nf, start_index );

  // We start at the outermost point
  q  = start_index;       
  for ( iv=1; iv<=nf; iv++ )    
  {
    tv1      = tr(iv,q-1);
    t1q(iv) /= tv1;
    tqn(iv) *= tv1;
    y(iv)    = ( y(iv) - s(iv,q-1)*(1-tv1)*(1+t1q(iv)*t1q(iv)*tv1) ) / tv1;
    k(iv,q)  = -lstep*( ( 2*yn(iv) + s(iv,q-1)*(1-2*tv1) )*t1q(iv)*t1q(iv)*tv1 
                        + y(iv) - s(iv,q-1) )*tv1/2;
  }

  // Points inside the LOS
  for ( q=start_index-1; q>1; q-- )
  {
    for ( iv=1; iv<=nf; iv++ )    
    {
      tv1      = tr(iv,q-1);    
      tv       = tr(iv,q);
      t1q(iv) /= tv1;
      y(iv)    = ( y(iv) - s(iv,q-1)*(1-tv1)*(1+t1q(iv)*t1q(iv)*tv1) ) / tv1;
      k(iv,q)  = -lstep*( ( 4*(yn(iv)-s(iv,q)) + 3*(s(iv,q)-s(iv,q-1)) + 
                            2*s(iv,q-1) )*t1q(iv)*t1q(iv)*tv1 + 
                     2*(y(iv)-s(iv,q-1))*tv1 + s(iv,q-1) - s(iv,q) )*tqn(iv)/2;
      tqn(iv) *= tv1;
      yn(iv)   = yn(iv)*tv + s(iv,q)*(1-tv);
    } 
  }

  // The tangent or ground point
  for ( iv=1; iv<=nf; iv++ )    
    k(iv,q)  = -lstep*( (2*yn(iv)*tv+s(iv,1)*(1-2*tv))*t1q(iv)*t1q(iv) +
                       yn(iv) - s(iv,1) ) * tqn(iv) / 2;

  // To check the function
  // Without ground reflection: T1Q=1 and Y=0
  // With ground reflection: T1Q=sqrt(1-e) and Y=eB
}



// ABSLOSWFS_DOWN
//
// This function covers 1D downward looking observations
//
// The expression used are described in sub-section 2.3 of the AUG section
// "Atmospheric WFs".
//
// Patrick Eriksson 15.06.00

void absloswfs_down (
                    MATRIX&   k,
              const VECTOR    y,
              const VECTOR    y_space,
              const int&      start_index,
              const int&      stop_index,
              const Numeric&  lstep,
              const MATRIX&   tr,
              const MATRIX&   s,
              const int&      ground,
              const VECTOR&   e_ground,
              const VECTOR&   y_ground )
{
  const size_t   nf = tr.dim(1); // number of frequencies
        size_t   iv;             // frequency index
           int   q;              // LOS index (same notation as in AUG)
        VECTOR   y0;             // see below
        MATRIX   k2;             // matrix for calling other LOS WFs functions
        VECTOR   tr0;            // see below

  // Resize K
  k.newsize( nf, start_index );

  // Calculate Y0, the intensity reaching the platform altitude at the far
  // end of LOS, that is, from above.
  y0 = y_space;
  rte_iterate( y0, start_index-1, stop_index, tr, s, nf );

  // Calculate TR0, the transmission from the platform altitude down to the
  // tangent point or the ground, and up to the platform again.
  bl( tr0, stop_index, stop_index, tr, ground, e_ground );


  // The indeces below STOP_INDEX are handled by the limb sounding function.
  // The limb function is given Y0 instead of cosmic radiation 
  absloswfs_limb( k2, y, y0, stop_index, lstep, tr, s, ground, e_ground );  
  for ( iv=1; iv<=nf; iv++ )
  {
    for ( q=1; q<stop_index; q++ )
      k(iv,q) = k2(iv,q);
  }

  // The indeces above STOP_INDEX are handled by the 1pass function.
  // The transmission below STOP_INDEX must here be considered.
  absloswfs_1pass( k2, y0, start_index, stop_index, lstep, tr, s, 
                     ground, e_ground, y_ground );
  for ( iv=1; iv<=nf; iv++ )
  {
    for ( q=stop_index+1; q<=start_index; q++ )
      k(iv,q) = k2(iv,q)*tr0(iv);
  }

  // The platform altitude must be treated seperately
  //
  // Calculate the intensity generated below point q-1, YQQ
  VECTOR yqq;
  rte( yqq, stop_index-1, stop_index-1, tr, s, VECTOR(nf,0.0), 
                                            ground, e_ground, y_ground );
  //
  // Y0 is moved one step upwards and TR0 one step downwards
  q = stop_index; 
  for ( iv=1; iv<=nf; iv++ )
  {
    tr0(iv) /= tr(iv,q-1)*tr(iv,q-1);    
    y0(iv)   = ( y0(iv) - s(iv,q)*(1-tr(iv,q)) ) / tr(iv,q);
    k(iv,q)  = -lstep*( (3*(y0(iv)-s(iv,q))*tr(iv,q-1)*tr(iv,q) + 
                      2*(s(iv,q)-s(iv,q-1))*tr(iv,q-1) + s(iv,q-1) )*tr0(iv) + 
                                           yqq(iv) - s(iv,q-1) )*tr(iv,q-1)/2;
  }
}



//==========================================================================
//=== Functions to calculate source function LOS WFs
//==========================================================================

// SOURCELOSWFS_1PASS
//
// This function covers cases where each point is passed once.
// That is 1D limb sounding and downward observatiuons are not covered.
//
// The expression used are described in sub-section 3.1 of the AUG section
// "Atmospheric WFs".
//
// Patrick Eriksson 11.09.00

void sourceloswfs_1pass (
                    MATRIX&   k,
              const int&      start_index,
	      const int&      stop_index,   // this variable is used by 1D down
              const MATRIX&   tr,
              const int&      ground,
              const VECTOR&   e_ground )
{
  const size_t   nf = tr.dim(1);      // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t(nf,1.0);           // transmission to the sensor
           int   q;                   // LOS index (same notation as in AUG) 

  if ( (ground==1) || (ground==start_index) )
    throw logic_error("The ground cannot be at one of the end points."); 

  // Resize K
  k.newsize( nf, start_index );

  // We start here at the LOS point closest to the sensor, that is,
  // reversed order compared to RTE_ITERATE  

  // The LOS point closest to the sensor 
  q  = stop_index;
  for ( iv=1; iv<=nf; iv++ )    
    k(iv,q)  = ( 1 - tr(iv,q) ) / 2;

  // Points inside the LOS
  for ( q=stop_index+1; q<start_index; q++ )
  {
    // Not a ground point
    if ( q != ground )
    {
      for ( iv=1; iv<=nf; iv++ )    
      {
        k(iv,q) = ( 1 - tr(iv,q-1)*tr(iv,q) ) * t(iv) / 2;
        t(iv)  *= tr(iv,q); 
      }
    }
    // A ground point
    else
    {
      out1 <<
        "WARNING: The function sourceloswfs_1pass not tested for ground reflections\n";
      for ( iv=1; iv<=nf; iv++ )    
      {
        k(iv,q) = ( (1-tr(iv,q))*(1-e_ground(iv))*tr(iv,q-1) + 1 - tr(iv,q-1) )
                                                            * t(iv) / 2;
        t(iv)  *= tr(iv,q)*(1-e_ground(iv)); 
      }
    }
  }

  // The LOS point furthest away from the sensor
  q = start_index;
  for ( iv=1; iv<=nf; iv++ )    
    k(iv,q)  = ( 1 - tr(iv,q-1) ) * t(iv) / 2;
}



// SOURCELOSWFS_LIMB
//
// This function covers 1D limb sounding
//
// The expression used are described in sub-section 3.2 of the AUG section
// "Atmospheric WFs".
//
// Patrick Eriksson 11.09.00

void sourceloswfs_limb (
                    MATRIX&   k,
              const int&      start_index,
              const MATRIX&   tr,
              const int&      ground,
              const VECTOR&   e_ground )
{
  const size_t   nf = tr.dim(1);      // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t1q;                 // transmission tangent point - q
        VECTOR   tqn(nf,1);           // transmission q - sensor
        int      q;                   // LOS index (same notation as in AUG)

  // Calculate the square root of the total transmission
  bl( t1q, start_index, start_index, tr, ground, e_ground );
  t1q = sqrt(t1q);

  // Resize K
  k.newsize( nf, start_index );

  // We start at the outermost point
  q  = start_index;       
  for ( iv=1; iv<=nf; iv++ )    
  {
    t1q(iv) /= tr(iv,q-1);
    k(iv,q)  = ( (1-tr(iv,q-1))*t1q(iv)*t1q(iv)*tr(iv,q-1) + 1 -tr(iv,q-1) )/2;
  }

  // Points inside the LOS
  for ( q=start_index-1; q>1; q-- )
  {
    for ( iv=1; iv<=nf; iv++ )    
    {
      t1q(iv) /= tr(iv,q-1);
      k(iv,q)  = ( (1-tr(iv,q-1)*tr(iv,q))*t1q(iv)*t1q(iv)*tr(iv,q-1)*tr(iv,q)+
                   (1-tr(iv,q-1))*tr(iv,q) + 1 - tr(iv,q) ) * tqn(iv) / 2;
      tqn(iv) *= tr(iv,q-1);
    } 
  }

  // The tangent or ground point
  for ( iv=1; iv<=nf; iv++ )    
    k(iv,q)  = ( (1-tr(iv,1))*(1+t1q(iv)*t1q(iv)*tr(iv,q)) ) * tqn(iv) / 2;
}



// SOURCELOSWFS_DOWN
//
// This function covers 1D downward looking observations
//
// The expression used are described in sub-section 3.3 of the AUG section
// "Atmospheric WFs".
//
// Patrick Eriksson 11.09.00

void sourceloswfs_down (
                    MATRIX&   k,
              const int&      start_index,
              const int&      stop_index,
              const MATRIX&   tr,
              const int&      ground,
              const VECTOR&   e_ground )
{
  const size_t   nf = tr.dim(1); // number of frequencies
        size_t   iv;             // frequency index
           int   q;              // LOS index (same notation as in AUG)
        MATRIX   k2;             // matrix for calling other LOS WFs functions
        VECTOR   tr0;            // see below

  // Resize K
  k.newsize( nf, start_index );

  // Calculate TR0, the transmission from the platform altitude down to the
  // tangent point or the ground, and up to the platform again.
  bl( tr0, stop_index, stop_index, tr, ground, e_ground );

  // The indeces below STOP_INDEX are handled by the limb sounding function.
  sourceloswfs_limb( k2, stop_index, tr, ground, e_ground );  
  for ( iv=1; iv<=nf; iv++ )
  {
    for ( q=1; q<stop_index; q++ )
      k(iv,q) = k2(iv,q);
  }

  // The indecies above STOP_INDEX are handled by the 1pass function.
  // The transmission below STOP_INDEX must here be considered.
  sourceloswfs_1pass( k2, start_index, stop_index, tr, ground, e_ground );
  for ( iv=1; iv<=nf; iv++ )
  {
    for ( q=stop_index+1; q<=start_index; q++ )
      k(iv,q) = k2(iv,q)*tr0(iv);
  }

  // The platform altitude must be treated seperately
  //
  // TR0 is moved one step downwards
  q = stop_index; 
  for ( iv=1; iv<=nf; iv++ )
  {
    tr0(iv) /= tr(iv,q-1)*tr(iv,q-1);    
    k(iv,q)  = ( (1-tr(iv,q-1)*tr(iv,q))*tr0(iv)*tr0(iv)*tr(iv,q-1) + 1 - 
                 tr(iv,q-1) ) / 2;
  }
}



// SOURCELOSWFS
//
// Patrick Eriksson 11.09.00

void sourceloswfs (
                    ARRAYofMATRIX&   sourceloswfs,
              const Los&             los,   
              const ARRAYofMATRIX&   trans,
              const VECTOR&          f_mono,
              const VECTOR&          e_ground,
              const VECTOR&          t_ground )
{
  const size_t  nza = los.start.dim();   // number of zenith angles  

  // Resize the LOS WFs array
  sourceloswfs.newsize(nza);

  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for (size_t i=1; i<=nza; i++ ) 
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Do something only if LOS has any points
    if ( los.p(i).dim() > 0 )
    {
      // The calculations are performed in 3 sub-functions
      //
      // Upward looking (=single pass)
      if ( los.stop(i)==1 )
        sourceloswfs_1pass( sourceloswfs(i), los.start(i), 1, trans(i), 
                                                     los.ground(i), e_ground );
      //
      // 1D limb sounding
      else if ( los.start(i) == los.stop(i) )
        sourceloswfs_limb( sourceloswfs(i), los.start(i), trans(i), 
                                                     los.ground(i), e_ground );
      //
      // 1D downward looking
      else 
        sourceloswfs_down( sourceloswfs(i), los.start(i), los.stop(i), 
                                           trans(i), los.ground(i), e_ground );
    }
  }
  out3 << "\n";
}



//==========================================================================
//=== Core functions for analytical WFs.
//===    Analytical expressions are used for species, temperature without 
//===    hydrostatic eq. and continuum absorption.
//===    These functions are very similar and a change in one of the functions
//===    should most likely be included in the other functions.
//==========================================================================

// K_SPECIES
//
// Patrick Eriksson 04.09.00

void k_species (
                    MATRIX&          ktot,
                    ARRAYofstring&   ktot_names,
                    MATRIX&          ktot_index,
                    MATRIX&          ktot_aux,
              const Los&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,             
              const TagGroups&       tags,
              const ARRAYofMATRIX&   abs_per_tg,
              const ARRAYofVECTOR&   vmrs,
              const VECTOR&          k_grid,
              const ARRAYofsizet&    tg_nr,
              const size_t&          unit )
{
  // Main sizes
  const size_t  nza = los.start.dim();     // number of zenith angles  
  const size_t  nv  = abs_per_tg(1).dim(1);// number of frequencies
  const size_t  ntg = tg_nr.dim();         // number of retrieval tags to do
  const size_t  np  = k_grid.dim();        // number of retrieval altitudes

  // -log(p) is used as altitude variable. The minus is included to get
  // increasing values, a demand for the grid functions. 
  const VECTOR  lgrid=-1.0*log(k_grid);    // -log of the retrieval pressures
        VECTOR  lplos;                     // -log of the LOS pressures

  // Indices
  // IP0 and IF0 are the index off-sets for the total K matrix
        size_t  itg;                       // Tag index
        size_t  iza;                       // Zenith angle index
        size_t  ip, ip0=0;                 // Retrieval point indices
        size_t  iv, iv0;                   // Frequency indices
        size_t  i1, iw;                    // weight indices

  // Other variables
        MATRIX  abs;                       // absorption at k_grid
        MATRIX  is;                        // matrix for storing LOS index
        VECTOR  w;                         // weights for LOS WFs
        VECTOR  a(nv);                     // temporary vector
        size_t  tg;                        // present tag nr
        VECTOR  vmr, p, t;                 // for conversion to VMR and ND 
        VECTOR  nd;                        // number density


  // Set-up local K and additional data. Set all values of K to 0
  MATRIX         k(nza*nv,ntg*np,0.0);
  ARRAYofstring  k_names(ntg);
  MATRIX         k_index(ntg,2);
  MATRIX         k_aux(ntg*np,3);

  // The calculations
  // Loop order:
  //   1 tags
  //   2 zenith angle
  //   3 retrieval altitudes 
  //   4 frequencies
  for ( itg=1; itg<=ntg; itg++ ) 
  {
    // Present tag nr
    tg = tg_nr(itg);

    // Check that the selected tag nr exist
    if ( tg <= 0 )
      throw runtime_error("The tag nr must be > 0."); 
    if ( tg > abs_per_tg.dim() )
      throw runtime_error("You have selected a non-existing tag nr."); 

    if ( ntg==1 )
      out2 << "  Doing tag " << tg << "\n";
    else
      out2 << "  Doing tag " << tg << " (" << itg << " of " << ntg << ")\n";

    // Fill K_NAMES and K_INDEX
    k_index(itg,1) = ip0+1;
    k_index(itg,2) = ip0+np;
    k_names(itg)   = tags(tg)(1).Name();

    // Interpolate to get the total absorption and the VMR values at the 
    // retrieval points and scale the absorption to the selected unit:
    //   1 normalised to a mean profile
    //   2 VMR
    //   3 number density
    interpp( abs, p_abs, abs_per_tg(tg), k_grid );
    interpp( vmr, p_abs, vmrs(tg), k_grid ); 
    if ( unit == 1 )
      for ( ip=1; ip<=np; ip++ )
        k_aux(ip0+ip,2) = 1.0;
    else if ( unit == 2 )
    {
      for ( ip=1; ip<=np; ip++ )
      {
        for ( iv=1; iv<=nv; iv++ )
          abs(iv,ip) /= vmr(ip);
        k_aux(ip0+ip,2) = vmr(ip);
      }
    }  
    else if ( unit == 3 )
    {
      interpp(  nd, p_abs, number_density(p_abs,t_abs), k_grid );
      nd = emult( vmr, nd );
      for ( ip=1; ip<=np; ip++ )
      {
        for ( iv=1; iv<=nv; iv++ )
          abs(iv,ip) /= nd(ip);
        k_aux(ip0+ip,2) = nd(ip);
      }
    }
    else
      throw runtime_error("Allowed retrieval units are 1-3."); 

    // Fill column 2 and 3 of K_AUX
    for ( ip=1; ip<=np; ip++ )
    {
       k_aux(ip0+ip,1) = k_grid(ip);
       k_aux(ip0+ip,3) = vmr(ip);
    }

    // Set frequency zenith angle index off-set to 0
    iv0 = 0;                 

    // Loop zenith angles
    out3 << "    Zenith angle nr:\n      ";
    for ( iza=1; iza<=nza; iza++ ) 
    {
      if ( (iza%20)==0 )
        out3 << "\n      ";
      out3 << " " << iza; cout.flush();
      
      // Do something only if LOS has any points
      if ( los.p(iza).dim() > 0 )
      {
        // Get the LOS points affected by each retrieval point
        lplos = -1.0 * log(los.p(iza));
        grid2grid_index( is, lplos, lgrid );

        // Loop retrieval points
        for ( ip=1; ip<=np; ip++ ) 
        {
          // Check if there is anything to do
          if ( is(ip,1) > 0 )
          {
            // Get the weights for the LOS points
	    grid2grid_weights( w, lplos, (size_t)is(ip,1), (size_t) is(ip,2), 
                                                                 lgrid, ip );

            // Calculate the WFs.
            // A is di/dkappa*dkappa/dkp in a compact form.
            // This is possible as the columns of dkappa/dkp are identical.  
            a = 0.0;                     
            i1 = (size_t)is(ip,1);       // first LOS point to consider
            for ( iv=1; iv<=nv; iv++ )
	    {
              for ( iw=i1; iw<=(size_t)is(ip,2); iw++ )
                a(iv) += absloswfs(iza)(iv,iw) * w(iw-i1+1);
              k(iv0+iv,ip0+ip) = a(iv) * abs(iv,ip);                    
	    }
          }            
        }
      }

      // Update the frequency index offset
      iv0 += nv;
    }  
    out3 << "\n";

    // Increase retrieval altitude index off-set 
    ip0 += np;
  }  

  // Combine total and local K matrices
  k_join( ktot, k, ktot_names, k_names, ktot_index, k_index, ktot_aux, k_aux );
}



// K_CONTABS
//
// Patrick Eriksson 06.09.00

void k_contabs (
                    MATRIX&          ktot,
                    ARRAYofstring&   ktot_names,
                    MATRIX&          ktot_index,
                    MATRIX&          ktot_aux,
              const Los&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          f_mono,
              const VECTOR&          k_grid,
              const size_t&          order )
{
  // Main sizes
  const size_t  nza = los.start.dim();     // number of zenith angles  
  const size_t  nv  = f_mono.dim();        // number of frequencies
  const size_t  np  = k_grid.dim();        // number of retrieval altitudes
  const size_t  npoints = order+1;         // number of off-set points

  // -log(p) is used as altitude variable. The minus is included to get
  // increasing values, a demand for the grid functions. 
  const VECTOR  lgrid=-1.0*log(k_grid);    // -log of the retrieval pressures
        VECTOR  lplos;                     // -log of the LOS pressures

  // Indices
  // IP0 and IF0 are the index off-sets for the total K matrix
        size_t  ipoint;                    // Off-set point index
        size_t  iza;                       // Zenith angle index
        size_t  ip, ip0=0;                 // Retrieval point indices
        size_t  iv, iv0;                   // Frequency indices
        size_t  i1, iw;                    // weight indices

  // Other variables
        VECTOR  fpoints;                   // frequencies of the off-set points
        VECTOR  b(nv);                     // fit base function
        MATRIX  is;                        // matrix for storing LOS index
        VECTOR  w;                         // weights for LOS WFs
        VECTOR  a(nv);                     // temporary vector


  // Check that the selected polynomial order
  if ( order < 0 )
    throw runtime_error("The polynomial order must be >= 0."); 

  // Set-up local K and additional data. Set all values of K to 0
  MATRIX         k(nza*nv,npoints*np,0.0);
  ARRAYofstring  k_names(npoints);
  MATRIX         k_index(npoints,2);
  MATRIX         k_aux(npoints*np,3);

  // Calculate the frequencies of the off-set points
  nlinspace( fpoints, f_mono(1), f_mono(nv), npoints );

  // The calculations
  // Loop order:
  //   1 polynomial order
  //   2 zenith angle
  //   3 retrieval altitudes 
  //   4 frequencies
  //
  // (Another loop order should be somewhat more efficient, but to keep this
  // function consistent with k_species, this loop order was selected.)
  //
  out2 << "  You have selected " << npoints << " off-set fit points.\n";
  for ( ipoint=1; ipoint<=npoints; ipoint++ ) 
  {
    out2 << "  Doing point " << ipoint << "\n";

    // Fill K_NAMES, K_INDEX and K_AUX
    k_index(ipoint,1) = ip0+1;
    k_index(ipoint,2) = ip0+np;
    k_names(ipoint)   = "Continuum absorption, point ?";
    for ( ip=1; ip<=np; ip++ )
    {
       k_aux(ip0+ip,1) = k_grid(ip);
       k_aux(ip0+ip,2) = 0.0;
       k_aux(ip0+ip,3) = 0.0;
    }

    // Set-up base vector for the present fit point 
    b = 1.0;
    if ( npoints > 1 )
    {
      for ( ip=1; ip<=npoints; ip++ )
      {
        if ( ip != ipoint )
	{ 
          for ( iv=1; iv<=nv; iv++ )
            b(iv) *= (f_mono(iv)-fpoints(ip)) / (fpoints(ipoint)-fpoints(ip));
	}
      }
    }

    // Set frequency zenith angle index off-set to 0
    iv0 = 0;                 

    // Loop zenith angles
    out3 << "    Zenith angle nr:\n      ";
    for ( iza=1; iza<=nza; iza++ ) 
    {
      if ( (iza%20)==0 )
        out3 << "\n      ";
      out3 << " " << iza; cout.flush();
      
      // Do something only if LOS has any points
      if ( los.p(iza).dim() > 0 )
      {
        // Get the LOS points affected by each retrieval point
        lplos = -1.0 * log(los.p(iza));
        grid2grid_index( is, lplos, lgrid );

        // Loop retrieval points
        for ( ip=1; ip<=np; ip++ ) 
        {
          // Check if there is anything to do
          if ( is(ip,1) > 0 )
          {
            // Get the weights for the LOS points
	    grid2grid_weights( w, lplos, (size_t)is(ip,1), (size_t) is(ip,2), 
                                                                 lgrid, ip );

            // Calculate the WFs.
            // A is di/dkappa*dkappa/dkp in a compact form.
            // This is possible as the columns of dkappa/dkp are identical.  
            a = 0.0;                     
            i1 = (size_t)is(ip,1);       // first LOS point to consider
            for ( iv=1; iv<=nv; iv++ )
	    {
              for ( iw=i1; iw<=(size_t)is(ip,2); iw++ )
                a(iv) += absloswfs(iza)(iv,iw) * w(iw-i1+1);
              k(iv0+iv,ip0+ip) = a(iv) * b(iv);                    
	    }
          }            
        }
      }

      // Update the frequency index offset
      iv0 += nv;
    }  
    out3 << "\n";

    // Increase retrieval altitude index off-set 
    ip0 += np;
  }  

  // Combine total and local K matrices
  k_join( ktot, k, ktot_names, k_names, ktot_index, k_index, ktot_aux, k_aux );
}



// K_TEMP_NOHYDRO
//
// Patrick Eriksson 06.09.00

void k_temp_nohydro (
                    MATRIX&          ktot,
                    ARRAYofstring&   ktot_names,
                    MATRIX&          ktot_index,
                    MATRIX&          ktot_aux,
              const Los&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          f_mono,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,            
              const ARRAYofVECTOR&   vmrs,
              const ARRAYofARRAYofLineRecord& lines_per_tg,
              const MATRIX&          abs,            
              const ARRAYofMATRIX&   trans,
              const VECTOR&          e_ground,
              const Numeric&         t_ground,
              const VECTOR&          k_grid )
{
  // Main sizes
  const size_t  nza = los.start.dim();     // number of zenith angles  
  const size_t  nv  = f_mono.dim();        // number of frequencies
  const size_t  np  = k_grid.dim();        // number of retrieval altitudes

  // -log(p) is used as altitude variable. The minus is included to get
  // increasing values, a demand for the grid functions. 
  const VECTOR  lgrid=-1.0*log(k_grid);    // -log of the retrieval pressures
        VECTOR  lplos;                     // -log of the LOS pressures

  // Indices
  // IP0 and IF0 are the index off-sets for the total K matrix
        size_t  iza;                       // zenith angle index
        size_t  ip;                        // retrieval point index
        size_t  iv, iv0=0;                 // frequency indices
        size_t  i1, iw;                    // weight indices

  // Other variables
        VECTOR  t;                         // temperature at retrieval points
        MATRIX  abs1k;                     // absorption for t_abs+1K
        MATRIX  dabs_dt;                   // see below
 ARRAYofMATRIX  abs_dummy;                 // dummy absorption array
 ARRAYofMATRIX  sloswfs;                   // source LOS WFs
	MATRIX  is;                        // matrix for storing LOS index
        VECTOR  w;                         // weights for LOS WFs
        VECTOR  a(nv), b(nv), pl;          // temporary vectors
       Numeric  c,d;                       // temporary values


  // Set-up local K and additional data. Set all values of K to 0
  MATRIX         k(nza*nv,np,0.0);
  ARRAYofstring  k_names(1,"Temperature");
  MATRIX         k_index(1,2);
  MATRIX         k_aux(np,3);
  k_index(1,1)   = 1;
  k_index(1,2)   = np;
  interpp( t, p_abs, t_abs, k_grid ); 
  for ( ip=1; ip<=np; ip++ )
  {
     k_aux(ip,1) = k_grid(ip);
     k_aux(ip,2) = t(ip);
     k_aux(ip,3) = 0.0;
  }

  // Calculate absorption for t_abs + 1K to estimate the temperature derivative
  // dabs_dt is the temperature derivative of the absorption at k_grid
  out2 << "  Calculating absorption for t_abs + 1K\n";
  out2 << "  ----- Messages from absCalc: -----\n";
  //
  absCalc( abs1k, abs_dummy, f_mono, p_abs, t_abs+1.0, vmrs, lines_per_tg );
  abs_dummy.clear();
  //
  out2 << "  ----- Back from absCalc ----------\n";
  //
  abs1k = abs1k - abs;
  interpp( dabs_dt, p_abs, abs1k, k_grid );
  abs1k.newsize(0,0);

  // Calculate source LOS WFs
  out2 << "  Calculating source LOS WFs\n";
  sourceloswfs( sloswfs, los, trans, f_mono, e_ground, t_ground );

  // Determine the temperatures at the retrieval points
  out2 << "  Calculating temperature at retrieval points\n";
  interpp( t, p_abs, t_abs, k_grid );

  // The calculations
  // Loop order:
  //   1 zenith angle
  //   2 retrieval altitudes 
  //   3 frequencies
  //
  out2 << "  Calculating the weighting functions\n";
  out3 << "    Zenith angle nr:\n      ";
  for ( iza=1; iza<=nza; iza++ ) 
  {
    if ( (iza%20)==0 )
      out3 << "\n      ";
    out3 << " " << iza; cout.flush();
    
    // Do something only if LOS has any points
    if ( los.p(iza).dim() > 0 )
    {
      // Get the LOS points affected by each retrieval point
      lplos = -1.0 * log(los.p(iza));
      grid2grid_index( is, lplos, lgrid );

      // Loop retrieval points
      for ( ip=1; ip<=np; ip++ ) 
      {
        // Check if there is anything to do
        if ( is(ip,1) > 0 )
        {
          // Get the weights for the LOS points
	  grid2grid_weights( w, lplos, (size_t)is(ip,1), (size_t) is(ip,2), 
                                                                 lgrid, ip );

          // Calculate the WFs.
          // A is di/dsigma*dsigma/dSp in a compact form.
          // B is di/dkappa*dkappa/dkp in a compact form.
          // This is possible as the columns of dkappa/dkp are identical and  
	  // that dkappa/dkp = dsigma/dSp
          // C is just a temporary value
	  // PL is the Planck function for the present temperature value
	  //
          a  = b = 0.0;                    
          c  = PLANCK_CONST / BOLTZMAN_CONST / t(ip);
          planck( pl, f_mono, t(ip) );
          i1 = (size_t)is(ip,1);       // first LOS point to consider
          //
          for ( iv=1; iv<=nv; iv++ )
	  {
            for ( iw=i1; iw<=(size_t)is(ip,2); iw++ )
	    {
              a(iv) += sloswfs(iza)(iv,iw) * w(iw-i1+1);
              b(iv) += absloswfs(iza)(iv,iw) * w(iw-i1+1);
	    }
            d = c * f_mono(iv);
            k(iv0+iv,ip) = a(iv) * d/t(ip) / (1-exp(-d)) * pl(iv) +
                                                      b(iv) * dabs_dt(iv,ip);
	  }
        }            
      }
    }

    // Update the frequency index offset
    iv0 += nv;
  }  
   out3 << "\n";

  // Combine total and local K matrices
  k_join( ktot, k, ktot_names, k_names, ktot_index, k_index, ktot_aux, k_aux );
}



//==========================================================================
//=== Workspace methods
//==========================================================================

// absloswfsCalc
//
// Patrick Eriksson 14.06.00

void absloswfsCalc (
                    ARRAYofMATRIX&   absloswfs,
              const Los&             los,   
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans,
              const VECTOR&          y,
              const VECTOR&          y_space,
              const VECTOR&          f_mono,
              const VECTOR&          e_ground,
              const Numeric&         t_ground )
{
  const size_t  nza = los.start.dim();   // number of zenith angles  
  const size_t  nf  = f_mono.dim();      // number of frequencies  
        VECTOR  yp(nf);                  // part of Y
        size_t  iy, iy0=0;               // y indices

  // Set up vector for ground blackbody radiation
  VECTOR   y_ground; 
  if ( any(los.ground) )
    planck( y_ground, f_mono, t_ground );

  // Resize the LOS WFs array
  absloswfs.newsize(nza);

  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for (size_t i=1; i<=nza; i++ ) 
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Do something only if LOS has any points
    if ( los.p(i).dim() > 0 )
    {

      // Pick out the part of Y corresponding to the present zenith angle
      for ( iy=1; iy<=nf; iy++ )    
        yp(iy) = y(iy0+iy);
      iy0 += nf;        

      // The calculations are performed in 3 sub-functions
      //
      // Upward looking (=single pass)
      if ( los.stop(i)==1 )
        absloswfs_1pass( absloswfs(i), yp, los.start(i), 1, los.l_step(i), 
                     trans(i), source(i), los.ground(i), e_ground, y_ground );
      //
      // 1D limb sounding
      else if ( los.start(i) == los.stop(i) )
        absloswfs_limb( absloswfs(i), yp, y_space, los.start(i), los.l_step(i),
                     trans(i), source(i), los.ground(i), e_ground );
      //
      // 1D downward looking
      else 
        absloswfs_down( absloswfs(i), yp, y_space, los.start(i), los.stop(i), 
                        los.l_step(i), trans(i), source(i), los.ground(i), 
                        e_ground, y_ground );
    }
  }
  out3 << "\n";
}



// absloswfsNoGround
//
// Patrick Eriksson 04.09.00

void absloswfsNoGround (
                    ARRAYofMATRIX&   absloswfs,
              const Los&             los,   
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans,
              const VECTOR&          y,
              const VECTOR&          y_space,
              const VECTOR&          f_mono )
{
  if ( any(los.ground) )  
    throw runtime_error("There is at least one intersection with the ground and this function cannot be used.");

  absloswfsCalc( absloswfs, los,source, trans, y, y_space, f_mono, VECTOR(0), 
                                                                        -1.0 );
}



// kInit
//
// Patrick Eriksson 05.09.00

void kInit (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_index,
                    MATRIX&          k_aux )
{
  k.newsize(0,0);
  k_names.newsize(0);
  k_index.newsize(0,0);
  k_aux.newsize(0,0);
}



// kSpecies
//
// Patrick Eriksson 04.09.00

void kSpecies (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_index,
                    MATRIX&          k_aux,
              const Los&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,             
              const TagGroups&       tags,
              const ARRAYofMATRIX&   abs_per_tg,
              const ARRAYofVECTOR&   vmrs,
              const VECTOR&          k_grid,
              const int&             nr,
              const int&             unit )
{
  ARRAYofsizet  tg_nr(1,nr);
  
  k_species( k, k_names, k_index, k_aux, los, absloswfs, p_abs, t_abs, 
                                 tags, abs_per_tg, vmrs, k_grid, tg_nr, unit );
}



// kSpeciesAll
//
// Patrick Eriksson 04.09.00

void kSpeciesAll (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_index,
                    MATRIX&          k_aux,
              const Los&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,             
              const TagGroups&       tags,
              const ARRAYofMATRIX&   abs_per_tg,
              const ARRAYofVECTOR&   vmrs,
              const VECTOR&          k_grid,
              const int&             unit )
{
  const size_t  ntg = abs_per_tg.dim();     // number of retrieval tags
  ARRAYofsizet  tg_nr;
  
  tg_nr.newsize(ntg);
  for ( size_t i=1; i<=ntg; i++ )
    tg_nr(i) = i;

  k_species( k, k_names, k_index, k_aux, los, absloswfs, p_abs, t_abs, 
                                 tags, abs_per_tg, vmrs, k_grid, tg_nr, unit );
}



// kContAbs
//
// Patrick Eriksson 06.09.00

void kContAbs (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_index,
                    MATRIX&          k_aux,
              const Los&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          f_mono,
              const VECTOR&          k_grid,
              const int&             order )
{
  k_contabs( k, k_names, k_index, k_aux, los, absloswfs, f_mono, k_grid, order );
}



// kTempNoHydro
//
// Patrick Eriksson 06.09.00

void kTempNoHydro (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_index,
                    MATRIX&          k_aux,
              const Los&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          f_mono,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,
              const ARRAYofVECTOR&   vmrs,
              const ARRAYofARRAYofLineRecord& lines_per_tg,
              const MATRIX&          abs,            
              const ARRAYofMATRIX&   trans,
              const VECTOR&          e_ground,
              const Numeric&         t_ground,
              const VECTOR&          k_grid )
{
  k_temp_nohydro( k, k_names, k_index, k_aux, los, absloswfs, f_mono, p_abs, 
          t_abs, vmrs, lines_per_tg, abs, trans, e_ground, t_ground, k_grid );
}



// kTempNoHydroNoGround
//
// Patrick Eriksson 06.09.00

void kTempNoHydroNoGround (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_index,
                    MATRIX&          k_aux,
              const Los&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          f_mono,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,
              const ARRAYofVECTOR&   vmrs,
              const ARRAYofARRAYofLineRecord& lines_per_tg,
              const MATRIX&          abs,            
              const ARRAYofMATRIX&   trans,
              const VECTOR&          k_grid )
{
  k_temp_nohydro( k, k_names, k_index, k_aux, los, absloswfs, f_mono, p_abs, 
          t_abs, vmrs, lines_per_tg, abs, trans, VECTOR(0), 0, k_grid );
}



void hTest(// WS Output:
           SPARSEMATRIX& h)
{
  SPARSEMATRIX A(10000,10000,3,"1.0 1.0 1.0","1 2 3","1 2 3");
  
  h = A;

  out1 << "  h = \n"
       << h << '\n';

}
