/* Copyright (C) 2002-2008 Patrick Eriksson <patrick.eriksson@chalmers.se>

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



/*===========================================================================
  === File description
  ===========================================================================*/

/*!
  \file   ppath.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2002-05-02
  
  \brief  Functions releated to calculation of propagation paths.
  
  Functions to determine propagation paths for different atmospheric
  dimensionalities, with and without refraction.

  The term propagation path is here shortened to ppath.
*/


// Values to apply if some calculation does not provide a solution
const double   R_NOT_FOUND   = -1;       // A value below zero
const double   L_NOT_FOUND   = 99e99;    // Some very large value for l/lat/lon
const double   LAT_NOT_FOUND = 99e99;
const double   LON_NOT_FOUND = 99e99;



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "agenda_class.h"
#include "array.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "math_funcs.h"
#include "messages.h"
#include "mystring.h"
#include "logic.h"
#include "poly_roots.h"
#include "ppath.h"
#include "refraction.h"
#include "special_interp.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



/*===========================================================================
  === Functions related to geometrical propagation paths
  ===========================================================================*/

//! geometrical_ppc
/*! 
   Calculates the propagation path constant for pure geometrical calculations.

   Both positive and negative zenith angles are handled.

   \return         Path constant.
   \param   r      Radius of the sensor position.
   \param   za     Zenith angle of the sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-05-17
*/
double geometrical_ppc( const double& r, const double& za )
{
  assert( r > 0 );
  assert( abs(za) <= 180 );

  return r * sin( DEG2RAD * abs(za) );
}



//! geompath_za_at_r
/*! 
   Calculates the zenith angle for a given radius along a geometrical 
   propagation path.

   For downlooking cases, the two points must be on the same side of 
   the tangent point.

   Both positive and negative zenith angles are handled.

   \return         Zenith angle at the point of interest.
   \param   ppc    Propagation path constant.
   \param   a_za   A zenith angle along the path on the same side of the 
                   tangent point as the point of interest.  
   \param   r      Radius of the point of interest.

   \author Patrick Eriksson
   \date   2002-05-17
*/
double geompath_za_at_r(
       const double&   ppc,
       const double&   a_za,
       const double&   r )
{
  assert( ppc >= 0 );
  assert( abs(a_za) <= 180 );
  assert( r >= ppc - RTOL );

  if( r > ppc )
    {
      double za = RAD2DEG * asin( ppc / r );
      if( abs(a_za) > 90 )
        { za = 180 - za; }
      if( a_za < 0 )
        { za = -za; }
      return za;
    }
  else
    {
      if( a_za > 0 )
        { return 90; }
      else
        { return -90; }
    }
}



//! geompath_r_at_za
/*! 
   Calculates the zenith angle for a given radius along a geometrical 
   propagation path.

   Both positive and negative zenith angles are handled.

   \return         Radius at the point of interest.
   \param   ppc    Propagation path constant.
   \param   za     Zenith angle at the point of interest.

   \author Patrick Eriksson
   \date   2002-06-05
*/
double geompath_r_at_za(
       const double&   ppc,
       const double&   za )
{
  assert( ppc >= 0 );
  assert( abs(za) <= 180 );

  return ppc / sin( DEG2RAD * abs(za) );
}



//! geompath_lat_at_za
/*!
   Calculates the latitude for a given zenith angle along a geometrical 
   propagation path.

   Positive and negative zenith angles are handled. A positive zenith angle
   means a movement towards higher latitudes.

   \return         The latitude of the second point.
   \param   za0    The zenith angle of the starting point.
   \param   lat0   The latitude of the starting point.
   \param   za     The zenith angle of the second point.

   \author Patrick Eriksson
   \date   2002-05-17
*/
double geompath_lat_at_za(
       const double&   za0,
       const double&   lat0,
       const double&   za )
{
  assert( abs(za0) <= 180 );
  assert( abs(za) <= 180 );
  assert( ( za0 >= 0 && za >= 0 )  ||  ( za0 < 0 && za < 0 ) );

  return lat0 + za0 - za;
}



//! geompath_l_at_r
/*!
   Calculates the length from the tangent point for the given radius.

   The tangent point is either real or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   \return         Length along the path from the tangent point.
   \param   ppc    Propagation path constant.
   \param   r      Radius of the point of concern.

   \author Patrick Eriksson
   \date   2002-05-20
*/
double geompath_l_at_r(
       const double&   ppc,
       const double&   r )
{
  assert( ppc >= 0 );
  assert( r >= ppc - RTOL );
  
  if( r > ppc )
    { return   sqrt( r*r - ppc*ppc ); }
  else
    { return   0; }
}



//! geompath_r_at_l
/*!
   Calculates the radius for a distance from the tangent point.

   The tangent point is either rwal or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   \return         Radius. 
   \param   ppc    Propagation path constant.
   \param   l      Length from the tangent point.

   \author Patrick Eriksson
   \date   2002-05-20
*/
double geompath_r_at_l(
       const double&   ppc,
       const double&   l )
{
  assert( ppc >= 0 );
  assert( l >= 0 );
  
  return sqrt( l*l + ppc*ppc );
}



//! geompath_r_at_lat
/*!
   Calculates the radius for a given latitude.

   \return         Radius at the point of interest.
   \param   ppc    Propagation path constant.
   \param   lat0   Latitude at some other point of the path.
   \param   za0    Zenith angle for the point with latitude lat0.
   \param   lat    Latitude of the point of interest.

   \author Patrick Eriksson
   \date   2002-06-05
*/
double geompath_r_at_lat(
       const double&   ppc,
       const double&   lat0,
       const double&   za0,
       const double&   lat )
{
  assert( ppc >= 0 );
  assert( abs(za0) <= 180 );
  assert( ( za0 >= 0 && lat >= lat0 )  ||  ( za0 <= 0 && lat <= lat0 ) );

  // Zenith angle at the new latitude
  const double za = za0 + lat0 -lat;

  return geompath_r_at_za( ppc, za );
}



//! geompath_from_r1_to_r2
/*!
   Determines radii, latitudes and zenith angles between two points of a 
   propagation path.

   Both start and end point are included in the returned vectors.

   \param   r      Output: Radius of propagation path points.
   \param   lat    Output: Latitude of propagation path points.
   \param   za     Output: Zenith angle of propagation path points.
   \param   lstep  Output: Distance along the path between the points. 
   \param   ppc    Propagation path constant.
   \param   r1     Radius for first point.
   \param   lat1   Latitude for first point.
   \param   za1    Zenith angle for first point.
   \param   r2     Radius for second point.
   \param   lmax   Length criterion for distance between path points.
                   A negative value means no length criterion.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void geompath_from_r1_to_r2( 
             Vector&   r,
             Vector&   lat,
             Vector&   za,
             double&   lstep,
       const double&   ppc,
       const double&   r1,
       const double&   lat1,
       const double&   za1,
       const double&   r2,
       const double&   lmax )
{
  // Calculate length along the path for point 1 and 2.
  const double l1 =  geompath_l_at_r( ppc, r1 );
  const double l2 =  geompath_l_at_r( ppc, r2 );
  
  // Calculate needed number of steps, considering a possible length criterion
  Index n;
  if( lmax > 0 )
    {
      // The absolute value of the length distance is needed here
      n = Index( ceil( abs( l2 - l1 ) / lmax ) );
 
      // We can't accept n=0, which is the case if l1 = l2
      if( n < 1 )
        { n = 1; }
    }
  else
    { n = 1; }

  // Length of path steps (note that lstep here can be negative)
  lstep = ( l2 - l1 ) / (double)n;

  // Allocate vectors and put in point 1
  r.resize(n+1);
  lat.resize(n+1);
  za.resize(n+1);
  r[0]   = r1;
  lat[0] = lat1;
  za[0]  = za1;

  // Loop steps (beside last) and calculate radius and zenith angle
  for( Index i=1; i<n; i++ )
    {
      r[i]   = geompath_r_at_l( ppc, l1 + lstep * (double)i );
      za[i]  = geompath_za_at_r( ppc, za1, r[i] );
    }

  // For maximum accuracy, set last radius to be exactly r2.
  r[n]   = r2;
  za[n]  = geompath_za_at_r( ppc, za1, r2 );  // Don't use r[n] here

  // Ensure that zenith and nadir observations keep their zenith angle
  if( abs(za1) < ANGTOL  ||  abs(za1) > 180-ANGTOL )
    { za = za1; }

  // Calculate latitudes
  for( Index i=1; i<=n; i++ )
    { lat[i] = geompath_lat_at_za( za1, lat1, za[i] ); }

  // Take absolute value of lstep
  lstep = abs( lstep );
}



//! za_geom2other_point
/*!
   Calculates the zenith angle for the geometrical propagation path between
   two specified points.

   The returned zenith angle is valid for point 1. That is, the propagation
   path goes from point 1 to point 2.

   \return         Zenith angle.
   \param   r1     Radius for point 1.
   \param   lat1   Latiytude for point 1.
   \param   r2     Radius for point 2.
   \param   lat2   Latitude for point 2.

   \author Patrick Eriksson
   \date   2002-07-03
*/
double za_geom2other_point(
       const double&    r1,
       const double&   lat1,
       const double&    r2,
       const double&   lat2 )
{
  if( lat2 == lat1 )
    {
      if( r1 <= r2 )
        { return 0; }
      else
        { return 180; }
    }
  else
    {
      // Absolute value of the latitude difference
      const double dlat = abs( lat2 - lat1 );

      // The zenith angle is obtained by a combination of the lawes of sine
      // and cosine.
      double za = dlat + RAD2DEG * asin( r1 * sin( DEG2RAD * dlat ) / 
                 sqrt( r1*r1 + r2*r2 - 2 * r1 * r2 * cos( DEG2RAD * dlat ) ) );

      // Consider viewing direction
      if( lat2 < lat1 )
        { za = -za; }

      return za;
    }
}





/*===========================================================================
  === Various functions
  ===========================================================================*/

//! resolve_lon
/*! 
   Resolves which longitude angle that shall be used.

   Longitudes are allowed to vary between -360 and 360 degress, while the
   inverse trigonomtric functions returns values between -180 and 180.
   This function determines if the longitude shall be shifted -360 or
   +360 to fit the longitudes set by the user.
   
   The argument *lon* as input is a value calculated by some inverse
   trigonometric function. The arguments *lon5* and *lon6* are the
   lower and upper limit for the probable range for *lon*. The longitude
   *lon* will be shifted with -360 or +360 degrees if lon is significantly
   outside *lon5* and *lon6*. No error is given if it is not possible to
   obtain a value between *lon5* and *lon6*. 

   The function exists both in a float and double version to avoid
   unnecessary copying of data for Numeric=float.

   \param   lon    In/Out: Longitude, possible shifted when returned.
   \param   lon5   Lower limit of probable range for lon.
   \param   lon6   Upper limit of probable range for lon

   \author Patrick Eriksson
   \date   2003-01-05
*/
void resolve_lon(
              double&   lon,
        const double&   lon5,
        const double&   lon6 )
{
  assert( lon6 >= lon5 ); 

  if( lon < lon5  &&  lon+180 <= lon6 )
    { lon += 360; }
  else if( lon > lon6  &&  lon-180 >= lon5 )
    { lon -= 360; }
}
#ifndef USE_DOUBLE
void resolve_lon(
              float&    lon,
        const double&   lon5,
        const double&   lon6 )
{
  assert( lon6 >= lon5 );

  if( lon < lon5  &&  lon+180 <= lon6 )
    { lon += 360; }
  else if( lon > lon6  &&  lon-180 >= lon5 )
    { lon -= 360; }
}
#endif



//! rotationmat3D
/*! 
   Creates a 3D rotation matrix

   Creates a rotation matrix such that R * x, operates on x by rotating 
   x around the origin a radians around line connecting the origin to the 
   point vrot.

   The function is based on rotationmat3D.m, by Belechi (the function 
   is found in atmlab).

   \param   R     Out: Rotation matrix
   \param   vrot  Rotation axis
   \param   a     Rotation angle

   \author Bileschi and Patrick Eriksson
   \date   2009-10-02
*/
void rotationmat3D( 
           Matrix&   R, 
   ConstVectorView   vrot, 
    const Numeric&   a )
{
  assert( R.ncols() == 3 );
  assert( R.nrows() == 3 );
  assert( vrot.nelem() == 3 );

  const double u = vrot[0];
  const double v = vrot[1];
  const double w = vrot[2];

  const double u2 = u * u;
  const double v2 = v * v;
  const double w2 = w * w;

  assert( sqrt( u2 + v2 + w2 ) );

  const double c = cos( DEG2RAD * a );
  const double s = sin( DEG2RAD * a );
  
  // Fill R
  R(0,0) = u2 + (v2 + w2)*c;
  R(0,1) = u*v*(1-c) - w*s;
  R(0,2) = u*w*(1-c) + v*s;
  R(1,0) = u*v*(1-c) + w*s;
  R(1,1) = v2 + (u2+w2)*c;
  R(1,2) = v*w*(1-c) - u*s;
  R(2,0) = u*w*(1-c) - v*s;
  R(2,1) = v*w*(1-c)+u*s;
  R(2,2) = w2 + (u2+v2)*c;
}



//! cart2zaaa
/*! 
   Converts a cartesian directional vector to zenith and azimuth

   This function and the sister function cart2zaaa handles
   transformation of line-of-sights. This in contrast to the sph/poslos
   functions that handles positions, or combinations of positions and
   line-of-sight.

   The cartesian coordinate system used for these two functions can 
   be defined as
    z : za = 0
    x : za=90, aa=0
    y : za=90, aa=90

   \param   za    Out: LOS zenith angle at observation position.
   \param   aa    Out: LOS azimuth angle at observation position.
   \param   dx    x-part of LOS unit vector.
   \param   dy    y-part of LOS unit vector.
   \param   dz    z-part of LOS unit vector.

   \author Patrick Eriksson
   \date   2009-10-02
*/
void cart2zaaa(
             double&   za,
             double&   aa,
       const double&   dx,
       const double&   dy,
       const double&   dz )
{
  const double r = sqrt( dx*dx + dy*dy + dz*dz );

  assert( r > 0 );

  za = RAD2DEG * acos( dz / r );
  aa = RAD2DEG * atan2( dy, dx );
}



//! zaaa2cart
/*! 
   Converts zenith and azimuth angles to a cartesian unit vector.

   This function and the sister function cart2zaaa handles
   transformation of line-of-sights. This in contrast to the sph/poslos
   functions that handles positions, or combinations of positions and
   line-of-sight.

   The cartesian coordinate system used for these two functions can 
   be defined as
    z : za = 0
    x : za=90, aa=0
    y : za=90, aa=90

   \param   dx    Out: x-part of LOS unit vector.
   \param   dy    Out: y-part of LOS unit vector.
   \param   dz    Out: z-part of LOS unit vector.
   \param   za    LOS zenith angle at observation position.
   \param   aa    LOS azimuth angle at observation position.

   \author Patrick Eriksson
   \date   2009-10-02
*/
void zaaa2cart(
             double&   dx,
             double&   dy,
             double&   dz,
       const double&   za,
       const double&   aa )
{
  const double   zarad  = DEG2RAD * za;
  const double   aarad  = DEG2RAD * aa;

  dz = cos( zarad );
  dx = sin( zarad );
  dy = sin( aarad ) * dx;
  dx = cos( aarad ) * dx;
}



/*! Maps MBLOCK_AA_GRID values to correct ZA and AA

   Sensor LOS azimuth angles and mblock_aa_grid values can not be added in a
   straightforward way due to properties of the polar coordinate system used to
   define line-of-sights. This function performs a "mapping" ensuring that the
   pencil beam directions specified by mblock_za_grid and mblock_aa_grid form
   a rectangular grid (on the unit sphere) for any za.

   za0 and aa0 match the angles of the ARTS WSV sensor_los.
   aa_grid shall hold values "close" to 0. The limit is here set to 5 degrees.

   \param   za         Out: Zenith angle matching aa0+aa_grid
   \param   aa         Out: Azimuth angles matching aa0+aa_grid
   \param   za0        Zenith angle
   \param   aa0        Centre azimuth angle
   \param   aa_grid    MBLOCK_AA_GRID values

   \author Patrick Eriksson
   \date   2009-10-02
*/
void map_daa(
             double&   za,
             double&   aa,
       const double&   za0,
       const double&   aa0,
       const double&   aa_grid )
{
  assert( abs( aa_grid ) <= 5 );

  Vector  xyz(3);
  Vector  vrot(3);
  Vector  u(3);

  // Unit vector towards aa0 at za=90
  //
  zaaa2cart( xyz[0], xyz[1], xyz[2], 90, aa0 );
    
  // Find vector around which rotation shall be performed
  // 
  // We can write this as cross([0 0 1],xyz). It turns out that the result 
  // of this operation is just [-y,x,0].
  //
  vrot[0] = -xyz[1];
  vrot[1] = xyz[0];
  vrot[2] = 0;

  // Unit vector towards aa0+aa at za=90
  //
  zaaa2cart( xyz[0], xyz[1], xyz[2], 90, aa0+aa_grid );

  // Apply rotation
  //
  Matrix R(3,3);
  rotationmat3D( R, vrot, za0-90 );
  mult( u, R, xyz );

  // Calculate za and aa for rotated u
  //
  cart2zaaa( za, aa, u[0], u[1], u[2] );
}






/*===========================================================================
  === Functions related to slope and tilt of the surface and pressure levels
  ===========================================================================*/

//! plevel_slope_2d
/*!
   Calculates the radial slope of the surface or a pressure level for 2D.

   The radial slope is here the derivative of the radius with respect to the
   latitude. The unit is accordingly m/degree.

   Note that the radius is defined to change linearly between grid points,
   and the slope is constant between to points of the latitude grid.

   Note also that the slope is always calculated with respect to increasing
   latitudes, independently of the zenith angle. The zenith angle is
   only used to determine which grid range that is of interest when the
   position is exactly on top of a grid point. 

   \return                 The radial slope [m/degree]
   \param   lat_grid       The latitude grid.
   \param   refellipsoid   As the WSV with the same name.
   \param   z_surf         Geometrical altitude of the surface, or the pressure
                           level of interest, for the latitide dimension
   \param   gp             Latitude grid position for the position of interest
   \param   za             LOS zenith angle.

   \author Patrick Eriksson
   \date   2002-06-03
*/
double plevel_slope_2d(
        ConstVectorView   lat_grid,           
        ConstVectorView   refellipsoid,
        ConstVectorView   z_surf,
        const GridPos&    gp,
        const double&     za )
{
  Index i1 = gridpos2gridrange( gp, abs( za ) >= 0 );
  const double r1 = refell2r( refellipsoid, lat_grid[i1] ) + z_surf[i1];
  const double r2 = refell2r( refellipsoid, lat_grid[i1+1] ) + z_surf[i1+1];
  return ( r2 - r1 ) / ( lat_grid[i1+1] - lat_grid[i1] );
}



//! plevel_slope_2d
/*!
   Calculates the radial slope of the surface or a pressure level for 2D.

   This function returns the same quantity as the function above, but takes
   the radius and latitude at two points of the pressure level, instead
   of vector input. That is, for this function the interesting latitude range
   is known when calling the function.

   \return         The radial slope [m/degree]
   \param   lat1   A latitude.
   \param   lat2   Another latitude.
   \param   r1     Radius at *lat1*.
   \param   r2     Radius at *lat2*.

   \author Patrick Eriksson
   \date   2002-12-21
*/
double plevel_slope_2d(
        const double&   lat1,
        const double&   lat2,
        const double&   r1,
        const double&   r2 )
{
  return   ( r2 - r1 ) / ( lat2 -lat1 );
}



//! rsurf_at_lat
/*!
   Determines the radius of a pressure level or the surface given the
   radius at the corners of a 2D grid cell.

   \return         Radius at the given latitude and longitude.
   \param   lat1   Lower latitude of grid cell.
   \param   lat3   Upper latitude of grid cell.
   \param   r1     Radius at *lat1*
   \param   r3     Radius at *lat3*
   \param   lat    Latitude for which radius shall be determined.

   \author Patrick Eriksson
   \date   2010-03-12
*/
double rsurf_at_lat(
       const double&   lat1,
       const double&   lat3,
       const double&   r1,
       const double&   r3,
       const double&   lat )
{
  return   r1 + ( lat - lat1 ) * ( r3 - r1 ) / ( lat3 - lat1 );
}


//! rsurf_at_latlon
/*!
   Determines the radius of a pressure level or the surface given the
   radius at the corners of a 3D grid cell.

   \return         Radius at the given latitude and longitude.
   \param   lat1   Lower latitude of grid cell.
   \param   lat3   Upper latitude of grid cell.
   \param   lon5   Lower longitude of grid cell.
   \param   lon6   Upper longitude of grid cell.
   \param   r15    Radius at crossing of *lat1* and *lon5*.
   \param   r35    Radius at crossing of *lat3* and *lon5*.
   \param   r36    Radius at crossing of *lat3* and *lon6*.
   \param   r16    Radius at crossing of *lat1* and *lon6*.
   \param   lat    Latitude for which radius shall be determined.
   \param   lon    Longitude for which radius shall be determined.

   \author Patrick Eriksson
   \date   2002-12-30
*/
double rsurf_at_latlon(
       const double&   lat1,
       const double&   lat3,
       const double&   lon5,
       const double&   lon6,
       const double&   r15,
       const double&   r35,
       const double&   r36,
       const double&   r16,
       const double&   lat,
       const double&   lon )
{
  // We can't have any assert of *lat* and *lon* here as we can go outside
  // the ranges when called from *plevel_slope_3d*.

  if( lat == lat1 )
    { return   r15 + ( lon - lon5 ) * ( r16 - r15 ) / ( lon6 -lon5 ); }
  else if( lat == lat3 )
    { return   r35 + ( lon - lon5 ) * ( r36 - r35 ) / ( lon6 -lon5 ); }
  else if( lon == lon5 )
    { return   r15 + ( lat - lat1 ) * ( r35 - r15 ) / ( lat3 -lat1 ); }
  else if( lon == lon6 )
    { return   r16 + ( lat - lat1 ) * ( r36 - r16 ) / ( lat3 -lat1 ); }
  else
    {
      const double   fdlat = ( lat - lat1 ) / ( lat3 - lat1 );
      const double   fdlon = ( lon - lon5 ) / ( lon6 - lon5 );
      return   (1-fdlat)*(1-fdlon)*r15 + fdlat*(1-fdlon)*r35 + 
                                         (1-fdlat)*fdlon*r16 + fdlat*fdlon*r36;
    }
}



//! plevel_slope_3d
/*!
   Calculates the local radial slope of the surface or a pressure level 
   for 3D.

   The function works basically as the non-vector version of
   *plevel_slope_2d*, but the position and viewing direction must
   here be specicified as the slope varies inside the cell grid, in
   constrast to a 2D latitude grid range.

   See further the other version of the function below.

   \return         The radial slope [m/degree]
   \param   lat1   Lower latitude of grid cell.
   \param   lat3   Upper latitude of grid cell.
   \param   lon5   Lower longitude of grid cell.
   \param   lon6   Upper longitude of grid cell.
   \param   r15    Radius at crossing of *lat1* and *lon5*.
   \param   r35    Radius at crossing of *lat3* and *lon5*.
   \param   r36    Radius at crossing of *lat3* and *lon6*.
   \param   r16    Radius at crossing of *lat1* and *lon6*.
   \param   lat    Latitude for which slope shall be determined.
   \param   lon    Longitude for which slope shall be determined.
   \param   aa     Azimuth angle for which slope shall be determined.

   \author Patrick Eriksson
   \date   2002-12-30
*/
double plevel_slope_3d(
        const double&   lat1,
        const double&   lat3,
        const double&   lon5,
        const double&   lon6,
        const double&   r15,
        const double&   r35,
        const double&   r36,
        const double&   r16,
        const double&   lat,
        const double&   lon,
        const double&   aa )
{
  // Size of test angular distance. Unit is degrees.
  const double   dang = 1e-5;

  // Radius at point of interest
  const double   r0 = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                                r15, r35, r36, r16, lat, lon );

  // Convert position and an imaginary LOS to cartesian coordinates
  // double is hard-coded to match declaration of poslos2cart and for
  // maximum precision  
  double   x, y, z, dx, dy, dz;
  poslos2cart( x, y, z, dx, dy, dz, r0, lat, lon, 90, aa );

  // Calculate the distance corresponding to *dang*
  const double   l = r0 * DEG2RAD * dang;

  // Calculate radius of the pressure level at the lat/lon, the distance
  // *l* away.
  double   r2, lat2, lon2, za2, aa2;
  cart2poslos( r2, lat2, lon2, za2, aa2, x+dx*l, y+dy*l, z+dz*l, dx, dy, dz );
  resolve_lon( lon2, lon5, lon6 );              
  r2 = rsurf_at_latlon( lat1, lat3, lon5, lon6, r15, r35, r36, r16, lat2,lon2);

  // Return slope
  return   ( r2 - r0 ) / dang;
}



//! plevel_slope_3d
/*!
   Calculates the radial slope of the surface or a pressure level for 3D.

   The radial slope is here the derivative of the radius with respect
   to an angular change (in degrees) along the great circle along the
   given azimuth angle. That is, how much the radius would change for a
   movement of r*pi/180 in the given azimuth angle (if the
   slope where constant along the distance). The unit is m/degree.

   For a point exactly on a grid value it is not clear if it is the
   range below or above that is of interest. The azimuth angle is used
   to resolve such cases.

   This function is in practice another way to call the non-vector version
   of the function above.

   \return                 The radial slope [m/degree]
   \param   lat_grid       The latitude grid.
   \param   lon_grid       The longitude grid.
   \param   refellipsoid   As the WSV with the same name.
   \param   z_surf         Geometrical altitude of the surface, or the pressure
                           level of interest.
   \param   gp_lat         Latitude grid position for the position of interest.
   \param   gp_lon         Longitude grid position for the position of interest.
   \param   aa             Azimuth angle.

   \author Patrick Eriksson
   \date   2002-06-03
*/
double plevel_slope_3d(
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,  
        ConstVectorView   refellipsoid,
        ConstMatrixView   z_surf,
        const GridPos&    gp_lat,
        const GridPos&    gp_lon,
        const double&     aa )
{
  Index ilat = gridpos2gridrange( gp_lat, abs( aa ) >= 0 );
  Index ilon = gridpos2gridrange( gp_lon, aa >= 0 );

  // Restore latitude and longitude values
  Vector  itw(2);
  double  lat, lon;
  interpweights( itw, gp_lat );
  lat = interp( itw, lat_grid, gp_lat );
  interpweights( itw, gp_lon );
  lon = interp( itw, lon_grid, gp_lon );

  // Extract values that defines the grid cell
  const double   lat1 = lat_grid[ilat];
  const double   lat3 = lat_grid[ilat+1];
  const double   lon5 = lon_grid[ilon];
  const double   lon6 = lon_grid[ilon+1];
  const double   re1  = refell2r( refellipsoid, lat1 );
  const double   re3  = refell2r( refellipsoid, lat3 );
  const double   r15  = re1 + z_surf(ilat,ilon);
  const double   r35  = re3 + z_surf(ilat+1,ilon);
  const double   r36  = re3 + z_surf(ilat+1,ilon+1);
  const double   r16  = re1 + z_surf(ilat,ilon+1);

  return   plevel_slope_3d( lat1, lat3, lon5, lon6, r15, r35, r36, r16, 
                                                                lat, lon, aa );
}



//! plevel_angletilt
/*!
   Calculates the angular tilt of the surface or a pressure level.

   Note that the tilt value is a local value. The tilt for a constant
   slope value, is different for different radii.

   \return        The angular tilt.
   \param    r    The radius for the level at the point of interest.
   \param    c    The radial slope, as returned by e.g. plevel_slope_2d.

   \author Patrick Eriksson
   \date   2002-06-03
*/
double plevel_angletilt(
        const double&   r,
        const double&   c )
{
  // The tilt (in radians) is c/r if c is converted to m/radian. So we get
  // conversion RAD2DEG twice
  return   RAD2DEG * RAD2DEG * c / r;
}



//! is_los_downwards
/*!
   Determines if a line-of-sight is downwards compared to the angular tilt
   of the surface or a pressure level.

   For example, this function can be used to determine if the line-of-sight
   goes into the surface for a starting point exactly on the surface radius.
  
   As the radius of the surface and pressure levels varies as a function of
   latitude, it is not clear if a zenith angle of 90 is above or below e.g.
   the surface.
 
   \return         Boolean that is true if LOS is downwards.
   \param   za     Zenith angle of line-of-sight.
   \param   tilt   Angular tilt of the surface or the pressure level (as
                   returned by plevel_angletilt)

   \author Patrick Eriksson
   \date   2002-06-03
*/
bool is_los_downwards( 
        const double&   za,
        const double&   tilt )
{
  assert( abs(za) <= 180 );

  // Yes, it shall be -tilt in both cases, if you wonder.
  if( za > (90-tilt)  ||  za < (-90-tilt) )
    { return true; }
  else
    { return false; }
}



//! r_crossing_2d
/*!
   Calculates where a 2D LOS crosses the specied radius.

   The function only looks for crossings in the forward direction of
   the given zenith angle (neglecting all solutions giving *l* <= 0).

   For cases with r_start <= r_hit and abs(za) > 90, the tangent point is
   passed and the returned crossing is on the other side of the tangent point.
 
   Both *lat* and *l* are set to a very high value if no crossing is found.

   \param   r         Out: Radius of found crossing.
   \param   lat       Out: Latitude of found crossing.
   \param   l         Out: Length along the path to the crossing.
   \param   r_hit     Radius of the level
   \param   r_start   Radius of start point.
   \param   lat_start Latitude of start point.
   \param   za_start  Zenith angle at start point.
   \param   ppc       Propagation path constant

   \author Patrick Eriksson
   \date   2012-02-18
*/
void r_crossing_2d(
             double&    lat,
             double&    l,
       const double&    r_hit,
       const double&    r_start,      
       const double&    lat_start,    
       const double&    za_start,     
       const double&    ppc )
{
  assert( abs( za_start ) <= 180 );
  assert( r_start >= ppc ); 

  const double absza = abs( za_start );

  // If above and looking upwards or r_hit below tangent point,
  // we have no crossing:
  if( ( r_start >= r_hit  &&  absza <= 90 )  ||  ppc > r_hit )
    { lat = LAT_NOT_FOUND;   l   = L_NOT_FOUND; }

  // Handle the cases of za=0 and za=180. 
  else if( absza < ANGTOL  ||  absza > 180-ANGTOL )
    {
      lat = lat_start;
      l   = abs( r_hit - r_start );
    }

  else
    {
      // Passages of tangent point
      if( absza > 90  &&  r_start <= r_hit )
        {
          double za  = geompath_za_at_r( ppc, sign(za_start)*89.99999, r_hit );
          lat = geompath_lat_at_za( za_start, lat_start, za );
          l = geompath_l_at_r( ppc, r_start ) + geompath_l_at_r( ppc, r_hit );
        }
      
      else
        {
          double za  = geompath_za_at_r( ppc, za_start, r_hit );
          lat = geompath_lat_at_za( za_start, lat_start, za );
          l = abs( geompath_l_at_r( ppc, r_start ) -
                   geompath_l_at_r( ppc, r_hit   ) );
        }
    }
}



//! r_crossing_3d
/*!
   Calculates where a 3D LOS crosses the specied radius

   The solution algorith is described in AUG. See the
   chapter on propagation paths.

   The function only looks for crossings in the forward direction of
   the given zenith angle (neglecting all solutions giving *l* <= 0).
   Note that the tangent point can be passed.
 
   The radius *r* is set to -1, lat and to 999 and *l* to avery high value if
   no crossing is found.

   \param   lat       Out: Latitude of found crossing.
   \param   lon       Out: Longitude of found crossing.
   \param   l         Out: Length along the path to the crossing.
   \param   r_hit     Traget radius.
   \param   r_start   Radius of start point.
   \param   lat_start Latitude of start point.
   \param   lon_start Longitude of start point.
   \param   za_start  Zenith angle at start point.
   \param   ppc       Propagation path constant
   \param   x         x-coordinate of start position.
   \param   y         y-coordinate of start position.
   \param   z         z-coordinate of start position.
   \param   dx        x-part of LOS unit vector.
   \param   dy        y-part of LOS unit vector.
   \param   dz        z-part of LOS unit vector.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void r_crossing_3d(
             double&    lat,
             double&    lon,
             double&    l,
       const double&    r_hit,
       const double&    r_start,    
       const double&    lat_start,  
       const double&    lon_start,  
       const double&    za_start,
       const double&    ppc,
       const double&    x,
       const double&    y,
       const double&    z,
       const double&    dx,
       const double&    dy,
       const double&    dz )
{
  assert( za_start >=   0 );
  assert( za_start <= 180 );

  // If above and looking upwards or r_hit below tangent point,
  // we have no crossing:
  if( ( r_start >= r_hit  &&  za_start <= 90 )  ||  ppc > r_hit )
    { lat = LAT_NOT_FOUND;   lon = LON_NOT_FOUND;   l   = L_NOT_FOUND; }

  // Handle the cases of za=0 and za=180. 
  else if( za_start < ANGTOL  ||  za_start > 180-ANGTOL )
    {
      lat = lat_start;
      lon = lon_start;
      l   = abs( r_hit - r_start );
    }

  else
    {
      const double   p  = x*dx + y*dy + z*dz;
      const double   pp = p * p;
      const double   q  = x*x + y*y + z*z - r_hit*r_hit;

      const double   l1 = -p + sqrt( pp - q );
      const double   l2 = -p - sqrt( pp - q );

      if( l1 < 0  &&  l2 > 0 )
        { l = l2; }
      else if( l1 > 0  &&  l2 < 0 )
        { l = l1; }
      else if( l1 < l2 )
        { l = l1; }
      else
        { l = l2; }

      if( l > 0 )
        {
          lat = RAD2DEG * asin( ( z+dz*l ) / r_hit );
          lon = RAD2DEG * atan2( y+dy*l, x+dx*l );
        }
    }
}



//! plevel_crossing_2d
/*!
   Calculates the angular distance to a crossing of a 2D pressure level
   or the surface.

   The function solves the problem mentioned above for a pressure
   level, or the surface, where the radius changes linearly as a
   function of latitude. No analytical solution to the original
   problem has been found. The problem involves sine and cosine of the
   latitude difference and these functions are replaced with their
   Taylor expansions where the two first terms are kept. This should
   be OK as in practical situations, the latitude difference inside a
   grid cell should not exceed 2 degrees, and the accuracy should be
   sufficient for values up to 3 degrees.

   The problem and its solution is further described in AUG. See the
   chapter on propagation paths.

   Both positive and negative zenith angles are handled. 

   The function only looks for crossings in the forward direction of
   the given zenith angle. This means that if rp>r0 and the absolute
   value of the zenith angle is < 90, no crossing will be found (if
   not the slope of the pressure level happen to be very strong).

   For downlooking cases, only the part down to the tangent point is
   considered.
   
   If the given path point is on the pressure level (rp=r0), the
   solution 0 is rejected.
 
   The latitude difference is set to 999 if no crossing exists.

   The variable names below are the same as in AUG.

   \return         The angular distance to the crossing.
   \param   rp     Radius of a point of the path inside the grid cell
   \param   za     Zenith angle of the path at rp.
   \param   r0     Radius of the pressure level or the surface at the
                   latitude of rp.
   \param   c      Linear slope term, as returned by plevel_slope_2d.

   \author Patrick Eriksson
   \date   2002-06-07
*/
double plevel_crossing_2d(
        const double&   rp,
        const double&   za,
        const double&   r0,
              double    c )
{
  assert( abs(za) <= 180 );

  const double no_crossing = 999;

  // Handle the cases of za=0 and za=180. 
  if( abs(za) < ANGTOL )
    {
      if( rp < r0 )
        { return 0; }
      else
        { return no_crossing; }
    }
  if( abs(za) > 180-ANGTOL )
    {
      if( rp > r0 )
        { return 0; }
      else
        { return no_crossing; }
    }

  // If rp>=r0, check if the given LOS goes in the direction towards
  // the pressure level.  If not, return 999.
  //
  if( rp >= r0   &&  !is_los_downwards( za,  plevel_angletilt( r0, c ) ) )
    { return no_crossing; }

  // The case with c=0 can be handled analytically
  if( c == 0 )
    {
      double   ppc = geometrical_ppc( rp, za );
      if( ( rp < r0  &&  abs(za) <= 90 )  ||  
                                 ( rp > r0  &&  abs(za) > 90  &&  r0 >= ppc ) )
        { return geompath_lat_at_za( za, 0, geompath_za_at_r( ppc, za, r0 ) );}
      else
        { return no_crossing; }
    }


  // Approximative solution:
  else
    {
      // If r0=rp, numerical inaccuracy can give a false solution, very close
      // to 0, that we must throw away.
      double   dmin = 0;
      if( r0 == rp )
        { dmin = 1e-12; }

      // The nadir angle in radians, and cosine and sine of that angle
      const double   beta = DEG2RAD * ( 180 - abs(za) );
      const double   cv = cos( beta );
      const double   sv = sin( beta );

      // Convert slope to m/radian and consider viewing direction
      c *= RAD2DEG;
      if( za < 0 )
        { c = -c; }

      // The vector of polynomial coefficients
      Vector p(5);
      //
      p[0] = ( r0 - rp ) * sv;
      p[1] = r0 * cv + c * sv;
      p[2] = -r0 * sv / 2 + c * cv;
      p[3] = -r0 * cv / 6 - c * sv / 2;
      p[4] = -c * cv / 6;

      // Calculate roots of the polynomial
      Matrix roots(4,2);
      poly_root_solve( roots, p );

      // Find the smallest root with imaginary part = 0, and real part > 0.
      //
      double dlat = no_crossing / RAD2DEG;
      //
      for( Index i=0; i<4; i++ )
        {
          if( roots(i,1) == 0  &&   roots(i,0) > dmin  &&  roots(i,0) < dlat )
            { dlat = roots(i,0); }
        }  

      // Change sign if zenith angle is negative
      if( dlat < no_crossing  &&  za < 0 )
        { dlat = -dlat; }

      return   RAD2DEG * dlat;
    }
}


//! rslope_crossing
/*!
   Calculates the angular distance to a crossing with a level having a 
   radial slope.

   The function solves the problem for a pressure level, or the planet's
   surface, where the radius changes linearly as a function of angle. No
   analytical solution to the original problem has been found. The problem
   involves sine and cosine of the latitude difference and these functions are
   replaced with their Taylor expansions where the two first terms are kept.
   This should be OK for most practical situations (the accuracy should be
   sufficient for values up to 1 degree?). 

   The problem and its solution is further described in AUG. See the chapter on
   propagation paths. The variable names below are the same as in AUG.

   Both positive and negative zenith angles are handled. 

   The function only looks for crossings in the forward direction of
   the given zenith angle. This means that if rp>r0 and the absolute
   value of the zenith angle is < 90, no crossing will be found (if
   not the slope of the pressure level happen to be very strong).

   If the given path point is on the pressure level (rp=r0), the solution 0 is
   rejected. The latitude difference is set to 999 if no crossing exists.

   No special action is taken for special cases such as c=0, za=0 and za=180,
   and for efficieny these cases should be catched by the calling function.

   \return         The angular distance to the crossing.
   \param   rp     Radius of a point of the path inside the grid cell
   \param   za     Zenith angle of the path at rp.
   \param   r0     Radius of the pressure level or the surface at the
                   latitude of rp.
   \param   c      Linear slope term, as returned by plevel_slope_2d.

   \author Patrick Eriksson
   \date   2002-06-07
*/
double rslope_crossing(
        const double&   rp,
        const double&   za,
        const double&   r0,
              double    c )
{
  // If r0=rp, numerical inaccuracy can give a false solution, very close
  // to 0, that we must throw away.
  double   dmin = 0;
  if( r0 == rp )
    { dmin = 1e-12; }

  // The nadir angle in radians, and cosine and sine of that angle
  const double   beta = DEG2RAD * ( 180 - abs(za) );
  const double   cv = cos( beta );
  const double   sv = sin( beta );

  // Convert slope to m/radian and consider viewing direction
  c *= RAD2DEG;
  if( za < 0 )
    { c = -c; }

  // The vector of polynomial coefficients
  Vector p(5);
  //
  p[0] = ( r0 - rp ) * sv;
  p[1] = r0 * cv + c * sv;
  p[2] = -r0 * sv / 2 + c * cv;
  p[3] = -r0 * cv / 6 - c * sv / 2;
  p[4] = -c * cv / 6;

  // Calculate roots of the polynomial
  Matrix roots(4,2);
  poly_root_solve( roots, p );
  
  // Find the smallest root with imaginary part = 0, and real part > 0.
  //
  double dlat = 999 / RAD2DEG;
  //
  for( Index i=0; i<4; i++ )
    {
      if( roots(i,1) == 0  &&  roots(i,0) > dmin  &&  roots(i,0) < dlat )
        { dlat = roots(i,0); }
    }  
  
  // Convert back to degrees
  // Change also sign if zenith angle is negative
  if( dlat < 3.14  &&  za < 0 )
    { dlat = -RAD2DEG * dlat; }
  else
    { dlat = RAD2DEG * dlat; }

  return   dlat;
}



//! plevel_crossing_2d
/*!
   Handles the crossing with a geometric ppaths step and a atmospheric 
   grid box level

   That is, we have a part of a pressure level or the planet's surface,
   extending between two latitudes (lat1 and lat3). The radius at each latitude
   is given (r1 and r3). The function first of determines if the ppath crosses
   the level/surface between the two latitudes. If yes, the radius and the
   latitude of the crossing point are calculated. 

   If the given path point is on the pressure level (rp=r0), the
   solution of zero length is rejected.
 

   \param   r           Out: Radius at crossing.
   \param   lat         Out: Latitude at crossing.
   \param   r_start0    In: Radius of start point.
   \param   lat_start   In: Latitude of start point.
   \param   za_start    In: LOS zenith angle at start point.
   \param   ppc         In: Propagation path constant.
   \param   lat1        In: Latitude of lower end.
   \param   lat3        In: Latitude of upper end.
   \param   r1          In: Radius at lat1.
   \param   r3          In: Radius at lat3.
   \param   above       In: True if ppath start point is above level. 
                        In: Otherwise false.

   \author Patrick Eriksson
   \date   2012-02-19
*/
void plevel_crossing_2d(
              double&   r,
              double&   lat,
        const double&   r_start0,
        const double&   lat_start,
        const double&   za_start,
        const double&   ppc,
        const double&   lat1,
        const double&   lat3,
        const double&   r1,
        const double&   r3,
        const bool&     above )
{
  const double absza = abs( za_start );

  assert( absza <= 180 );
  assert( lat_start >=lat1  &&  lat_start <= lat3 );

  // The case of no slope
  if( r1 == r3 )
    {
      // Set r_start, considering impact of numerical problems
      double r_start = r_start0;
      if( above )
        { if( r_start < r1 ) { r_start = r1; } }
      else
        { if( r_start > r1 ) { r_start = r1; } }

      double l;
      r = r1;
      r_crossing_2d( lat, l, r, r_start, lat_start, za_start, ppc );

      // Check if inside [lat1,lat3]
      if( lat > lat3  ||  lat < lat1 )
        { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND; }  
    }

  // With slope
  else
    {
      const double rmin = min( r1, r3 );
      const double rmax = max( r1, r3 );

      // Set r_start, considering impact of numerical problems
      double r_start = r_start0;
      if( above )
        { if( r_start < rmax ) { r_start = rmax; } }
      else
        { if( r_start > rmin ) { r_start = rmin; } }
      
      double l, za=999;

      // Calculate crossing with closest radius
      if( r_start > rmax )
        {
          r = rmax;
          r_crossing_2d( lat, l, r, r_start, lat_start, za_start, ppc );
        }
      else if( r_start < rmin )
        {
          r = rmin;      
          r_crossing_2d( lat, l, r, r_start, lat_start, za_start, ppc );
        }
      else
        { r = r_start; lat = lat_start; l = 0; za = za_start; }
  
      // latt must be inside [lat1,lat3] if relevant to continue
      if( lat < lat1  ||  lat > lat3 )
        { r = R_NOT_FOUND; }   // lat already set by r_crossing_2d

      // Otherwise continue from found point, considering the level slope 
      else
        {
          // Radius and slope at lat
          const double  cpl = plevel_slope_2d( lat1, lat3, r1, r3 );
          const double  rpl = r1 + cpl*(lat-lat1);

          // Make adjustment if numerical problems
          if( above )
            { if( r < rpl ) { r = rpl; } }
          else
            { if( r > rpl ) { r = rpl; } }

          // za_start = 0
          if( abs(za_start) < ANGTOL )
            {
              if( r >= rpl )
                { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND; }  
              else
                { r = rpl; lat = lat_start; }
            }

          // za_start = 180
          else if( abs(za_start) > 180-ANGTOL )
            {
              if( r <= rpl )
                { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND; }  
              else
                { r = rpl; lat = lat_start; }
            }
          
          // Not nadir or zenith
          else
            {
              // Calculate zenith angle at (r,lat) (if <180, already set above)
              if( za > 180 )  // lat+za preserved (also with negative za)
                { za = lat_start + za_start -lat; };

              // Latitude distance from present point to actual crossing
              const double dlat = rslope_crossing( r, za, rpl, cpl );

              // Update lat and check if still inside [lat1,lat3].
              // If yes, determine r
              lat += dlat;
              if( lat < lat1  ||  lat > lat3 )
                { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND; }  
              else
                { 
                  // Check if the correction ever exceed 0.2 deg?
                  assert( abs(dlat) < 0.2 );
                  r = r1 + cpl*(lat-lat1); 
                }
            }  
        }
    }
}




/*===========================================================================
  === Path through grid cells and 1D pressure range
  ===========================================================================*/

//! do_gridrange_1d
/*!
   Calculates the geometrical path through a 1D grid range.

   This function works as *do_gridcell_2d*, but is valid for 1D cases.

   The coding of variables and end face is as for *do_gridcell_2d*, with
   the excpetion that end faces 2 and 4 do not exist here.

   \param   r_v         Out: Vector with radius of found path points.
   \param   lat_v       Out: Vector with latitude of found path points.
   \param   za_v        Out: Vector with LOS zenith angle at found path points.
   \param   lstep       Out: Vector with length along the path between points.
   \param   endface     Out: Number coding for exit face. See above.
   \param   r_start0    Radius of start point.
   \param   lat_start   Latitude of start point.
   \param   za_start    LOS zenith angle at start point.
   \param   ppc         Propagation path constant.
   \param   lmax        Maximum allowed length along the path. -1 = no limit.
   \param   ra          Radius of lower pressure level.
   \param   rb          Radius of upper pressure level (rb > ra);
   \param   rsurface    Radius for the surface.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void do_gridrange_1d(
              Vector&   r_v,
              Vector&   lat_v,
              Vector&   za_v,
              double&   lstep,
              Index&    endface,
        const double&   r_start0,
        const double&   lat_start,
        const double&   za_start,
        const double&   ppc,
        const double&   lmax,
        const double&   ra,
        const double&   rb,
        const double&   rsurface )
{
  double   r_start = r_start0;

  assert( rb > ra );
  assert( r_start >= ra - RTOL );
  assert( r_start <= rb + RTOL );

  // Shift radius if outside
  if( r_start < ra )
    { r_start = ra; }
  else if( r_start > rb )
    { r_start = rb; }

  endface = -1;

  // If upward, end point radius is always rb
  if( za_start <= 90 )
    { 
      endface = 3;
      // Calculate basic variables from r_start to rb:
      geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, 
                                               lat_start, za_start, rb, lmax );
    }
  else
    {
      // Path reaches ra:
      if( ra > rsurface  &&  ra > ppc ) 
        {
          endface = 1;
          geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, 
                                               lat_start, za_start, ra, lmax );
        }
      // Path reaches the surface:
      else if( rsurface > ppc )
        {
          endface = 7;
          geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, 
                                         lat_start, za_start, rsurface, lmax );
        }
      // The tangent point is passed, on the way back to rb
      else
        {
          endface = 1;
          Vector rx, latx, zax, ry, laty, zay;
          geompath_from_r1_to_r2( rx, latx, zax, lstep, ppc, r_start, 
                                              lat_start, za_start, ppc, lmax );
          const Index lx = rx.nelem()-1;
          geompath_from_r1_to_r2( ry, laty, zay, lstep, ppc, rx[lx], 
                                                 latx[lx], zax[lx], rb, lmax );
          const Index ly = ry.nelem()-1, n  = lx + ly + 1;
          //
          r_v.resize(n); lat_v.resize(n); za_v.resize(n);
          //
          for( Index i=0; i<=lx; i++ )
            { r_v[i] = rx[i];   lat_v[i] = latx[i];   za_v[i] = zax[i]; }
          for( Index i=0; i<=ly; i++ )
            { 
              const Index iv = lx + i;
              r_v[iv] = ry[i];   lat_v[iv] = laty[i];   za_v[iv] = zay[i]; 
            }
        }
    }
  assert( endface > 0 );
}



//! do_gridcell_2d
/*!
   Calculates the geometrical path through a 2D grid cell.

   The function determines the geometrical path from the given start
   point to the boundary of the grid cell. The face where the path
   exits the grid cell is denoted as the end face. The following
   number coding is used for the variable *endface*: <br>
   1: The face at the lower latitude point. <br>
   2: The face at the lower (geometrically) pressure level. <br>
   3: The face at the upper latitude point. <br>
   4: The face at the upper (geometrically) pressure level. <br>
   7: The end point is an intersection with the surface. 

   The corner points are names r[lat][a,b]. For example: r3b.
   The latitudes are numbered to match the end faces. This means that
   the lower latitude has number 1, and the upper number 3. The pressure
   levels are named as a and b: <br>
   a: Lower pressure level (highest pressure). <br>
   b: Upper pressure level (lowest pressure).

   Path points are included if *lmax*>0 and the distance to the end
   point is > than *lmax*.

   The return vectors (*r_v* etc.) can have any length when handed to
   the function.

   \param   r_v         Out: Vector with radius of found path points.
   \param   lat_v       Out: Vector with latitude of found path points.
   \param   za_v        Out: Vector with LOS zenith angle at found path points.
   \param   lstep       Out: Vector with length along the path between points.
   \param   endface     Out: Number coding for exit face. See above.
   \param   r_start0    Radius of start point.
   \param   lat_start   Latitude of start point.
   \param   za_start    LOS zenith angle at start point.
   \param   ppc         Propagation path constant.
   \param   lmax        Maximum allowed length along the path. -1 = no limit.
   \param   lat1        Latitude of left end face (face 1) of the grid cell.
   \param   lat3        Latitude of right end face (face 3) of the grid cell.
   \param   r1a         Radius of lower-left corner of the grid cell.
   \param   r3a         Radius of lower-right corner of the grid cell.
   \param   r3b         Radius of upper-right corner of the grid cell (r3b>r3a)
   \param   r1b         Radius of upper-left corner of the grid cell (r1b>r1a).
   \param   rsurface1   Radius for the surface at *lat1*.
   \param   rsurface3   Radius for the surface at *lat3*.

   \author Patrick Eriksson
   \date   2002-11-28
*/
void do_gridcell_2d(
              Vector&   r_v,
              Vector&   lat_v,
              Vector&   za_v,
              double&   lstep,
              Index&    endface,
        const double&   r_start0,
        const double&   lat_start,
        const double&   za_start,
        const double&   ppc,
        const double&   lmax,
        const double&   lat1,
        const double&   lat3,
        const double&   r1a,
        const double&   r3a,
        const double&   r3b,
        const double&   r1b,
        const double&   rsurface1,
        const double&   rsurface3 )
{
  double   r_start = r_start0;

  // Assert latitudes
  assert( lat_start >= lat1 );
  assert( lat_start <= lat3 );

  // Slopes of pressure levels
  const double  c2 = plevel_slope_2d( lat1, lat3, r1a, r3a );
  const double  c4 = plevel_slope_2d( lat1, lat3, r1b, r3b );
  const double  csurface = plevel_slope_2d( lat1, lat3, rsurface1, rsurface3 );

  // Latitude distance between start point and left grid cell boundary
  const double dlat_left  = lat_start - lat1;

  // Radius of lower and upper pressure level at *lat_start*.
  const double   rlow = r1a + c2*dlat_left;
  const double   rupp = r1b + c4*dlat_left;

  // Latitude distance to latitude end face in the viewing direction
  double dlat_endface;
  if( za_start >= 0 )
    { dlat_endface = lat3 - lat_start; }
  else
    { dlat_endface = -dlat_left; }

  // Assert radius (some extra tolerance is needed for radius)
  assert( r_start >= rlow - RTOL );
  assert( r_start <= rupp + RTOL );

  // Shift radius if outside
  if( r_start < rlow )
    { r_start = rlow; }
  else if( r_start > rupp )
    { r_start = rupp; }

  // The end point is most easily determined by the latitude difference to 
  // the start point. For some cases there exist several crossings and we want 
  // the one closest in latitude to the start point. The latitude distance 
  // for the crossing shall not exceed dlat2end.
  //
  double   dlat2end = 999;
  double   abs_za_start = abs( za_start );
           endface = 0;


  // --- Lower face 
  //
  // We don't need to consider the face if we are standing
  // on the pressure level.
  //
  if( r_start > rlow  &&  abs(za_start) > ANGTOL )
    {
      dlat2end = plevel_crossing_2d( r_start, za_start, rlow, c2 );
      endface = 2;  // This variable will be re-set if there was no crossing
    }


  // --- The surface.
  //
  // Check shall be done only if the surface is, at least partly, inside 
  // the grid cell.
  //
  if( rsurface1 >= r1a  ||  rsurface3 >= r3a )
    {
      double r_surface = rsurface1 + csurface * dlat_left;
     
      assert( r_start >= r_surface );

      double dlat2surface = plevel_crossing_2d( r_start, za_start, r_surface,
                                                                     csurface );
      if( abs(dlat2surface) <= abs(dlat2end) )
        {
          dlat2end = dlat2surface;
          endface  = 7;
        }
    }


  // If dlat2end <= dlat_endface we are ready. Otherwise we have to check
  // remaining cell faces. The same applies after testing upper face.


  // --- Upper face  (pressure level ip+1).
  //
  if( abs(dlat2end) > abs(dlat_endface)  &&  abs_za_start < 180-ANGTOL )
    {
      // For cases when the tangent point is in-between *r_start* and
      // the pressure level, 999 is returned. This case will anyhow
      // be handled correctly.

      double dlattry = plevel_crossing_2d( r_start, za_start, rupp, c4 );

      if( abs(dlattry) < abs(dlat2end) )
        {
          dlat2end = dlattry;
          endface  = 4;
        }
    }

  // Left or right end face
  if( abs(dlat2end) > abs(dlat_endface) )
    { 
      dlat2end = dlat_endface; 
      if( za_start >= 0 )
        { endface  = 3; }
      else
        { endface  = 1; }
    }
  assert( endface );

  // Check if there is a tangent point inside the grid cell. A tangent point
  // requires special treatment
  if( abs_za_start > 90  &&  ( abs_za_start - abs(dlat2end) ) < 90 ) 
    {  
      Vector rx, latx, zax, ry, laty, zay;
      geompath_from_r1_to_r2( rx, latx, zax, lstep, ppc, r_start, 
                                              lat_start, za_start, ppc, lmax );
      const Index lx = rx.nelem()-1;
      do_gridcell_2d( ry, laty, zay, lstep, endface, ppc, latx[lx], zax[lx], 
                      ppc, lmax, lat1, lat3, r1a, r3a, r3b, r1b, 
                      rsurface1, rsurface3 );
      const Index ly = ry.nelem()-1, n  = lx + ly + 1;
      //
      r_v.resize(n); lat_v.resize(n); za_v.resize(n);
      //
      for( Index i=0; i<=lx; i++ )
        { r_v[i] = rx[i];   lat_v[i] = latx[i];   za_v[i] = zax[i]; }
      for( Index i=0; i<=ly; i++ )
        { 
          const Index iv = lx + i;
          r_v[iv] = ry[i];   lat_v[iv] = laty[i];   za_v[iv] = zay[i]; 
        }
    }

  // No tangent point
  else
    {
      // Calculate radius for end point.
      // To obtain best possible accuracy it is calculated to match found,
      // end face, and not based on dlat2end.
      //
      double   r_end = -1;
      //
      if( endface == 1 )
        { r_end = geompath_r_at_lat( ppc, lat_start, za_start, lat1 ); }
      else if( endface == 2 )
        { r_end = r1a + c2 * ( dlat_left + dlat2end ); }
      else if( endface == 3 )
        { r_end = geompath_r_at_lat( ppc, lat_start, za_start, lat3 ); }
      else if( endface == 4 )
        { r_end = r1b + c4 * ( dlat_left + dlat2end ); }
      else if( endface == 7 )
        { r_end = rsurface1 + csurface * ( dlat_left + dlat2end ); }

      // Fill the return vectors
      geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, lat_start, 
                                                       za_start, r_end, lmax );

      // Force end latitude to be exact as possible
      if( endface == 1 )
        { lat_v[lat_v.nelem()-1] = lat1; }
      else if( endface == 3 )
        { lat_v[lat_v.nelem()-1] = lat3; }
      else
        { lat_v[lat_v.nelem()-1] = lat_start + dlat2end; }
    }
}
void do_gridcell_2dv2(
              Vector&   r_v,
              Vector&   lat_v,
              Vector&   za_v,
              double&   lstep,
              Index&    endface,
        const double&   r_start,
        const double&   lat_start,
        const double&   za_start,
        const double&   ppc,
        const double&   lmax,
        const double&   lat1,
        const double&   lat3,
        const double&   r1a,
        const double&   r3a,
        const double&   r3b,
        const double&   r1b,
        const double&   rsurface1,
        const double&   rsurface3 )
{
  // Assert latitudes
  assert( lat_start >= lat1 );
  assert( lat_start <= lat3 );

  // Radius end latitude of end point
  double r=-1, lat;

  endface = 0;

  // Check if crossing with lower pressure level
  plevel_crossing_2d( r, lat, r_start, lat_start, za_start, ppc, lat1, lat3, 
                                                              r1a, r3a, true );
  if( r > 0 )
    { endface = 2; }


  // Check if crossing with surface
  if( rsurface1 >= r1a  ||  rsurface3 >= r3a )
    {
      double rt, latt; 
      plevel_crossing_2d( rt, latt, r_start, lat_start, za_start, ppc, lat1, 
                                            lat3, rsurface1, rsurface3, true );
      if( rt > 0  &&  abs(latt-lat_start) <= abs(lat-lat_start) )
        { endface = 7;   r = rt;   lat = latt; }
    }

  // If crossing found (r>0) we are ready!
  
  // Upper pressure level
  if( r <= 0 )
    {
      double rt, latt ;
      plevel_crossing_2d( rt, latt, r_start, lat_start, za_start, ppc, lat1, 
                                                       lat3, r1b, r3b, false );
      if( latt < LAT_NOT_FOUND )
        { endface = 4;   r = rt;   lat = latt; }
    }
  
  // Latitude endfaces
  if( r <= 0 )
    {
      if( za_start < 0 )
        { endface = 1;  lat = lat1; }
      else
        { endface = 3;  lat = lat3; }
      r = geompath_r_at_lat( ppc, lat_start, za_start, lat ); 
    }

  assert (endface );

  // Check if there is a tangent point inside the grid cell. A tangent point
  // requires special treatment
  const double absza = abs( za_start );
  if( absza > 90  &&  ( absza - abs(lat_start-lat) ) < 90 ) 
    {  
      Vector rx, latx, zax, ry, laty, zay;
      geompath_from_r1_to_r2( rx, latx, zax, lstep, ppc, r_start, lat_start, 
                                                         za_start, ppc, lmax );
      const Index lx = rx.nelem()-1;
      do_gridcell_2dv2( ry, laty, zay, lstep, endface, ppc, latx[lx], zax[lx], 
                        ppc, lmax, lat1, lat3, r1a, r3a, r3b, r1b, 
                        rsurface1, rsurface3 );
      const Index ly = ry.nelem()-1, n  = lx + ly + 1;
      //
      r_v.resize(n); lat_v.resize(n); za_v.resize(n);
      //
      for( Index i=0; i<=lx; i++ )
        { r_v[i] = rx[i];   lat_v[i] = latx[i];   za_v[i] = zax[i]; }
      for( Index i=0; i<=ly; i++ )
        { 
          const Index iv = lx + i;
          r_v[iv] = ry[i];   lat_v[iv] = laty[i];   za_v[iv] = zay[i]; 
        }
    }

  // No tangent point
  else
    {
      // Fill the return vectors
      geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, lat_start, 
                                                           za_start, r, lmax );

      // Re-set last latitude, for higher accuracy (possibly?)
      lat_v[lat_v.nelem()-1] = lat;
    }
}



//! do_gridcell_3d
/*!
   Calculates the geometrical path through a 3D grid cell.

   The function determines the geometrical path from the given start
   point to the boundary of the grid cell. The face where the path
   exits the grid cell is denoted as the end face. The same number
   coding as in *do_gridcell_2d* is used, where the additional longitude
   end faces are numbered as: <br>
   5: The face at the lower longitude point. <br>
   6: The face at the upper longitude point.

   The corner points are numbered as *do_gridcell_2d*, but 5 or 6 is added
   after the latitude number to indicate the longitude. This means that
   r16a, is the corner at lat1, lon6 and pressure level a.

   See further *do_gridcell_2d*.

   \param   r_v         Out: Vector with radius of found path points.
   \param   lat_v       Out: Vector with latitude of found path points.
   \param   lon_v       Out: Vector with longitude of found path points.
   \param   za_v        Out: Vector with LOS zenith angle at found path points.
   \param   aa_v        FIXME: Add documentation.
   \param   lstep       Out: Vector with length along the path between points.
   \param   endface     Out: Number coding for exit face. See above.
   \param   r_start0    Radius of start point.
   \param   lat_start0  Latitude of start point.
   \param   lon_start0  Longitude of start point.
   \param   za_start    LOS zenith angle at start point.
   \param   aa_start    LOS azimuth angle at start point.
   \param   ppc         Propagation path constant.
   \param   lmax        Maximum allowed length along the path. -1 = no limit.
   \param   lat1        Latitude of left end face (face 1) of the grid cell.
   \param   lat3        Latitude of right end face (face 3) of the grid cell.
   \param   lon5        Lower longitude limit of the grid cell.
   \param   lon6        Upper longitude limit of the grid cell.
   \param   r15a        Radius of corner: lower p-level,*lat1* and *lon5*.
   \param   r35a        Radius of corner: lower p-level,*lat3* and *lon5*.
   \param   r36a        Radius of corner: lower p-level,*lat3* and *lon6*.
   \param   r16a        Radius of corner: lower p-level,*lat1* and *lon6*.
   \param   r15b        Radius of corner: upper p-level,*lat1* and *lon5*.
   \param   r35b        Radius of corner: upper p-level,*lat3* and *lon5*.
   \param   r36b        Radius of corner: upper p-level,*lat3* and *lon6*.
   \param   r16b        Radius of corner: upper p-level,*lat1* and *lon6*.
   \param   rsurface15   Radius for the surface at *lat1* and *lon5*.
   \param   rsurface35   Radius for the surface at *lat3* and *lon5*.
   \param   rsurface36   Radius for the surface at *lat3* and *lon6*.
   \param   rsurface16   Radius for the surface at *lat1* and *lon6*.

   \author Patrick Eriksson
   \date   2002-11-28
*/
void do_gridcell_3d(
              Vector&   r_v,
              Vector&   lat_v,
              Vector&   lon_v,
              Vector&   za_v,
              Vector&   aa_v,
              double&   lstep,
              Index&    endface,
        const double&   r_start0, 
        const double&   lat_start0,
        const double&   lon_start0,
        const double&   za_start,
        const double&   aa_start,
        const double&   ppc,
        const double&   lmax,
        const double&   lat1,
        const double&   lat3,
        const double&   lon5,
        const double&   lon6,
        const double&   r15a,
        const double&   r35a,
        const double&   r36a,
        const double&   r16a,
        const double&   r15b,
        const double&   r35b,
        const double&   r36b,
        const double&   r16b,
        const double&   rsurface15,
        const double&   rsurface35,
        const double&   rsurface36,
        const double&   rsurface16 )
{
  double   r_start   = r_start0;
  double   lat_start = lat_start0;
  double   lon_start = lon_start0;

  // Assert latitude and longitude
  assert( lat_start >= lat1 - LATLONTOL );
  assert( lat_start <= lat3 + LATLONTOL );
  assert( !( abs( lat_start) < 90  &&  lon_start < lon5 - LATLONTOL ) );
  assert( !( abs( lat_start) < 90  &&  lon_start > lon6 + LATLONTOL ) );

  // Shift latitude and longitude if outside
  if( lat_start < lat1 )
    { lat_start = lat1; }
  else if( lat_start > lat3 )
    { lat_start = lat3; }
  if( lon_start < lon5 )
    { lon_start = lon5; }
  else if( lon_start > lon6 )
    { lon_start = lon6; }

  // Radius of lower and upper pressure level at the start position
  double   rlow = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                r15a, r35a, r36a, r16a, lat_start, lon_start );
  double   rupp = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                r15b, r35b, r36b, r16b, lat_start, lon_start );

  // Assert radius (some extra tolerance is needed for radius)
  assert( r_start >= rlow - RTOL );
  assert( r_start <= rupp + RTOL );

  // Shift radius if outside
  if( r_start < rlow )
    { r_start = rlow; }
  else if( r_start > rupp )
    { r_start = rupp; }

  // Position and LOS in cartesian coordinates
  double   x, y, z, dx, dy, dz;
  poslos2cart( x, y, z, dx, dy, dz, r_start, lat_start, lon_start, 
                                                          za_start, aa_start );

  // Determine the position of the end point
  //
  endface  = 0;
  bool tanpoint = false;
  //
  double   r_end, lat_end, lon_end, l_end;

  // Zenith and nadir looking are handled as special cases

  // Zenith looking
  if( za_start < ANGTOL )
    {
      r_end   = rupp;
      lat_end = lat_start;
      lon_end = lon_start;
      l_end   = rupp - r_start;
      endface  = 4;
    }

  // Nadir looking
  else if( za_start > 180-ANGTOL )
    {
      const double   rsurface = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                               rsurface15, rsurface35, rsurface36, rsurface16, 
                                                 lat_start, lon_start );
      
      if( rlow > rsurface )
        {
          r_end  = rlow;
          endface = 2;
        }
      else
        {
          r_end  = rsurface;
          endface = 7;
        }
      lat_end = lat_start;
      lon_end = lon_start;
      l_end   = r_start - r_end;
    }

  // Check faces in number order
  else
    {
      // Calculate correction terms for the position to compensate for
      // numerical inaccuracy. 
      //
      double   r_corr, lat_corr, lon_corr;
      //
      cart2sph( r_corr, lat_corr, lon_corr, x, y, z );
      //
      r_corr   -= r_start;
      lat_corr -= lat_start;
      lon_corr -= lon_start;

      // The end point is found by testing different step lengths until the
      // step length has been determined by a precision of *l_acc*.
      //
      // The first step is to found a length that goes outside the grid cell, 
      // to find an upper limit. The lower limit is set to 0. The upper
      // limit is found be testing lengths of 100, 1000 ... m.
      // The search algorith is bisection, the next length to test is the
      // mean value of the minimum and maximum length limits.

               l_end  = 100;

      double   l_acc  = 1e-3;
      double   l_in   = 0, l_out = l_end;
      bool     ready  = false, startup = true;
      double   abs_aa = abs( aa_start );

      double   l_tan = 99e6;
      if( za_start > 90 )
        { l_tan = sqrt( r_start*r_start - ppc*ppc ); }

      bool   do_surface = false;
      if( rsurface15+RTOL >= r15a  ||  rsurface35+RTOL >= r35a  ||  
          rsurface36+RTOL >= r36a  ||  rsurface16+RTOL >= r16a )
        { do_surface = true; }

      while( !ready )
        {
          cart2sph( r_end, lat_end, lon_end, 
                                          x+dx*l_end, y+dy*l_end, z+dz*l_end );
          r_end   -= r_corr;
          lat_end -= lat_corr;
          lon_end -= lon_corr;
          resolve_lon( lon_end, lon5, lon6 );              

          // Special fixes for north-south observations
          if( abs( lat_start ) < 90  &&  ( abs(aa_start) < ANGTOL  ||  
                                           abs(aa_start) > 180-ANGTOL ) )
            { lon_end = lon_start; }
          
          // Special fixes for west-east observations
          if( abs( abs_aa - 90 ) < ANGTOL )
            {
              if( lat_start == 0 )
                { lat_end = 0; }
              else if( lat_start > 0  &&  lat_end > lat_start )
                { lat_end = lat_start; }
              else if( lat_start < 0  &&  lat_end < lat_start )
                { lat_end = lat_start; }
            }

          bool   inside = true;

          rlow = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                r15a, r35a, r36a, r16a, lat_end, lon_end );
          rupp = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                r15b, r35b, r36b, r16b, lat_end, lon_end );
          
          if( do_surface )
            {
              const double   r_surface = 
                           rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                           rsurface15, rsurface35, rsurface36, rsurface16, 
                                                        lat_end, lon_end );

              if( r_surface+RTOL >= rlow  &&  r_end <= r_surface+RTOL )
                { inside = false;   endface = 7; }
            }

          if( inside ) 
            {
              if( lat_end < lat1 )
                { inside = false;   endface = 1; }
              else if( lat_end > lat3 )
                { inside = false;   endface = 3; }
              else if( lon_end < lon5 )
                { inside = false;   endface = 5; }
              else if( lon_end > lon6 )
                { inside = false;   endface = 6; }
              else if( r_end < rlow )
                { inside = false;   endface = 2; }
              else if( r_end > rupp )
                { inside = false;   endface = 4; }
            }              

          if( startup )
            {
              if( inside )
                { 
                  if( l_end < l_tan  &&  l_end*10 > l_tan )
                    { l_end = l_tan; }
                  else
                    { l_end *= 10; }
                }
              else
                { 
                  l_out = l_end;   
                  l_end = ( l_out + l_in ) / 2;
                  startup = false; 
                }
            }
          else
            {
              if( inside )
                { l_in = l_end; }
              else
                { l_out = l_end; }
                            
              if( ( l_out - l_in ) < l_acc )
                { ready = true; }
              else
                { l_end = ( l_out + l_in ) / 2; }
            }
        }

      // Now when we are ready, we remove the correction terms. Otherwise
      // we can end up in an infinite loop if the step length is smaller
      // than the correction.
      r_end   += r_corr;
      lat_end += lat_corr;
      lon_end += lon_corr;
      resolve_lon( lon_end, lon5, lon6 );              

      // Set the relevant coordinate to be consistent with found endface.
      //
      if( endface == 1 )
        { lat_end = lat1; }
      else if( endface == 2 )
        { r_end = rsurf_at_latlon( lat1, lat3, lon5, lon6, r15a, r35a, r36a, 
                                                    r16a, lat_end, lon_end ); }
      else if( endface == 3 )
        { lat_end = lat3; }
      else if( endface == 4 )
        { r_end = rsurf_at_latlon( lat1, lat3, lon5, lon6, r15b, r35b, r36b, 
                                                    r16b, lat_end, lon_end ); }
      else if( endface == 5 )
        { lon_end = lon5; }
      else if( endface == 6 )
        { lon_end = lon6; }
      else if( endface == 7 )
        { r_end = rsurf_at_latlon( lat1, lat3, lon5, lon6, rsurface15,  
                      rsurface35, rsurface36, rsurface16, lat_end, lon_end ); }

      // Check if there is a tangent point inside the grid cell. A tangent 
      // point requires special treatment, which is done last in function
      if( za_start > 90  &&  l_tan < l_end )
        {
          geompath_tanpos_3d( r_end, lat_end, lon_end, l_end, r_start, 
                              lat_start, lon_start, za_start, aa_start, ppc );
          resolve_lon( lon_end, lon5, lon6 );              
          tanpoint = 1;
        }
    }

  //--- Create return vectors
  //
  Index n = 1;
  //
  if( lmax > 0 )
    {
      n = Index( ceil( abs( l_end / lmax ) ) );
      if( n < 1 )
        { n = 1; }
    }
  //
  r_v.resize( n+1 );
  lat_v.resize( n+1 );
  lon_v.resize( n+1 );
  za_v.resize( n+1 );
  aa_v.resize( n+1 );
  //
  r_v[0]   = r_start;
  lat_v[0] = lat_start;
  lon_v[0] = lon_start;
  za_v[0]  = za_start;
  aa_v[0]  = aa_start;
  //
  double  ldouble = l_end / (double)n;
  lstep = ldouble;
  // 
  for( Index j=1; j<=n; j++ )
    {
      const double   l  = ldouble * (double)j;

      // We seperate here float/double to avoid unnecessary copying for
      // Numeric=double. A solution like this is needed to avoid compiler
      // warnings.
#ifdef USE_DOUBLE
      cart2poslos( r_v[j], lat_v[j], lon_v[j], za_v[j], aa_v[j],
                                          x+dx*l, y+dy*l, z+dz*l, dx, dy, dz );
#else
      double   rtmp, lattmp, lontmp, zatmp, aatmp;
      cart2poslos( rtmp, lattmp, lontmp, zatmp, aatmp,
                                          x+dx*l, y+dy*l, z+dz*l, dx, dy, dz );
      r_v[j]   = rtmp;
      lat_v[j] = lattmp;
      lon_v[j] = lontmp;
      za_v[j]  = zatmp;
      aa_v[j]  = aatmp;
#endif
   }

  //--- Set last point especially, which should improve the accuracy
  r_v[n]   = r_end; 
  lat_v[n] = lat_end;
  lon_v[n] = lon_end;

  //--- Set last zenith angle to be as accurate as possible
  if( za_start < ANGTOL  ||  za_start > 180-ANGTOL )
    { za_v[n] = za_start; }
  else if( tanpoint )
    { za_v[n] = 90; }
  else
    { za_v[n] = geompath_za_at_r( ppc, za_start, r_v[n] ); }

  //--- Set last azimuth angle to be as accurate as possible for
  //    zenith and nadir observations
  if( abs( lat_start ) < 90  &&  
          ( abs(aa_start) < ANGTOL  ||  abs( aa_start) > 180-ANGTOL ) )
    {  
      aa_v[n]  = aa_start; 
      lon_v[n] = lon_start;
    }

  // Shall lon values be shifted (value 0 and n+1 are already OK)?
  for( Index j=1; j<n; j++ )
    { resolve_lon( lon_v[j], lon5, lon6 ); }

  // Make second part if above ending at tangent point 
  if( tanpoint )
    {
      Vector rx, latx, lonx, zax, aax, ry, laty, lony, zay, aay;
      rx   = r_v;
      latx = lat_v;
      lonx = lon_v;
      zax  = za_v;
      aax  = aa_v;
      const Index lx = rx.nelem()-1;
      do_gridcell_3d( ry, laty, lony, zay, aay, lstep, endface, ppc, latx[lx], 
                      lonx[lx], zax[lx], aax[lx], ppc, lmax, lat1, lat3, 
                      lon5, lon6, r15a, r35a, r36a, r16a, r15b, r35b, r36b, 
                      r16b, rsurface15, rsurface35, rsurface36, rsurface16 );
      const Index ly = ry.nelem()-1, nxy  = lx + ly + 1;
      //
      r_v.resize(nxy);  lat_v.resize(nxy); lon_v.resize(nxy); 
      za_v.resize(nxy); aa_v.resize(nxy);
      //
      for( Index i=0; i<=lx; i++ )
        { 
          r_v[i]  = rx[i];   lat_v[i] = latx[i];  lon_v[i] = lonx[i];   
          za_v[i] = zax[i];  aa_v[i]  = aax[i]; 
        }
      for( Index i=0; i<=ly; i++ )
        { 
          const Index iv = lx + i;
          r_v[iv]  = ry[i];   lat_v[iv] = laty[i];  lon_v[iv] = lony[i]; 
          za_v[iv] = zay[i];  aa_v[iv] = aay[i]; 
        }
    }
}





/*===========================================================================
  === Functions operating on the Ppath structure
  ===========================================================================*/

//! ppath_init_structure
/*!
   Initiates a Ppath structure to hold the given number of points.

   All fields releated with the surface, symmetry and tangent point are set
   to 0 or empty. The background field is set to background case 0. The
   constant field is set to -1. The refraction field is set to 0.

   The length of the l_step field is set to np-1.

   \param   ppath            Output: A Ppath structure.
   \param   atmosphere_dim   The atmospheric dimensionality.
   \param   np               Number of points of the path.

   \author Patrick Eriksson
   \date   2002-05-17
*/
void ppath_init_structure( 
              Ppath&      ppath,
        const Index&      atmosphere_dim,
        const Index&      np )
{
  assert( atmosphere_dim >= 1 );
  assert( atmosphere_dim <= 3 );

  ppath.dim        = atmosphere_dim;
  ppath.np         = np;
  ppath.constant   = -1;   
  if( atmosphere_dim < 3 )
    {
      ppath.pos.resize( np, 2 );
      ppath.los.resize( np, 1 );
    }
  else
    {
      ppath.pos.resize( np, atmosphere_dim );
      ppath.los.resize( np, 2 );
      ppath.gp_lon.resize( np );
    }
  ppath.gp_p.resize( np );
  if( atmosphere_dim >= 2 )
    { ppath.gp_lat.resize( np ); }
  ppath.r.resize( np );
  if( np > 0 )
    { ppath.l_step.resize( np-1 ); }
  else
    { ppath.l_step.resize( 0 ); }
  ppath_set_background( ppath, 0 );
  ppath.nreal.resize( np );
  ppath.lspace = 0;
}




//! ppath_set_background 
/*!
   Sets the background field of a Ppath structure.

   The different background cases have a number coding to simplify a possible
   change of the strings and checking of the what case that is valid.

   The case numbers are:                    <br>
      0. Not yet set.                       <br>
      1. Space.                             <br>
      2. The surface.                       <br>
      3. The cloud box boundary.            <br>
      4. The interior of the cloud box.       

   \param   ppath            Output: A Ppath structure.
   \param   case_nr          Case number (see above)

   \author Patrick Eriksson
   \date   2002-05-17
*/
void ppath_set_background( 
              Ppath&      ppath,
        const Index&      case_nr )
{
  switch ( case_nr )
    {
    case 0:
      ppath.background = "";
      break;
    case 1:
      ppath.background = "space";
      break;
    case 2:
      ppath.background = "surface";
      break;
    case 3:
      ppath.background = "cloud box level";
      break;
    case 4:
      ppath.background = "cloud box interior";
      break;
    default:
      ostringstream os;
      os << "Case number " << case_nr << " is not defined.";
      throw runtime_error(os.str());
    }
}



//! ppath_what_background
/*!
   Returns the case number for the radiative background.

   See further the function *ppath_set_background*.

   \return                   The case number.
   \param   ppath            A Ppath structure.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Index ppath_what_background( const Ppath&   ppath )
{
  if( ppath.background == "" )
    { return 0; }
  else if( ppath.background == "space" )
    { return 1; }
  else if( ppath.background == "surface" )
    { return 2; }
  else if( ppath.background == "cloud box level" )
    { return 3; }
  else if( ppath.background == "cloud box interior" )
    { return 4; }
  else
    {
      ostringstream os;
      os << "The string " << ppath.background 
         << " is not a valid background case.";
      throw runtime_error(os.str());
    }
}



//! ppath_copy
/*!
   Copy the content in ppath2 to ppath1.

   The ppath1 structure must be allocated before calling the function. The
   structure can be allocated to hold more points than found in ppath2.
   The data of ppath2 is placed in the first positions of ppath1.

   \param   ppath1    Output: Ppath structure.
   \param   ppath2    The ppath structure to be copied.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_copy(
           Ppath&      ppath1,
     const Ppath&      ppath2 )
{
  assert( ppath1.np >= ppath2.np ); 

  // The field np shall not be copied !

  ppath1.dim        = ppath2.dim;
  ppath1.constant   = ppath2.constant;
  ppath1.background = ppath2.background;
  ppath1.lspace     = ppath2.lspace;

  ppath1.pos(Range(0,ppath2.np),joker) = ppath2.pos;
  ppath1.los(Range(0,ppath2.np),joker) = ppath2.los;
  ppath1.r[Range(0,ppath2.np)]         = ppath2.r;
  ppath1.nreal[Range(0,ppath2.np)]     = ppath2.nreal;
  if( ppath2.np > 1 )
    { ppath1.l_step[Range(0,ppath2.np-1)]  = ppath2.l_step; }

  for( Index i=0; i<ppath2.np; i++ )
    {
      gridpos_copy( ppath1.gp_p[i], ppath2.gp_p[i] );
      
      if( ppath1.dim >= 2 )
        { gridpos_copy( ppath1.gp_lat[i], ppath2.gp_lat[i] ); }
      
      if( ppath1.dim == 3 )
        { gridpos_copy( ppath1.gp_lon[i], ppath2.gp_lon[i] ); }
     }
}



//! ppath_append
/*!
   Combines two Ppath structures   

   The function appends a Ppath structure to another structure. 
 
   All the data of ppath1 is kept.

   The first point in ppath2 is assumed to be the same as the last in ppath1.
   Only data in ppath from the fields pos, los, z, l_step, gp_XXX and 
   background are considered.

   \param   ppath1    Output: Ppath structure to be expanded.
   \param   ppath2    The Ppath structure to include in ppath.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_append(
           Ppath&   ppath1,
     const Ppath&   ppath2 )
{
  const Index n1 = ppath1.np;
  const Index n2 = ppath2.np;

  Ppath   ppath;
  ppath_init_structure( ppath, ppath1.dim, n1 );
  ppath_copy( ppath, ppath1 );

  ppath_init_structure( ppath1, ppath1.dim, n1 + n2 - 1 );
  ppath_copy( ppath1, ppath );

  // Append data from ppath2
  Index i1;
  for( Index i=1; i<n2; i++ )
    { 
      i1 = n1 + i - 1;

      ppath1.pos(i1,0)      = ppath2.pos(i,0);
      ppath1.pos(i1,1)      = ppath2.pos(i,1);
      ppath1.los(i1,0)      = ppath2.los(i,0);
      ppath1.r[i1]          = ppath2.r[i];
      ppath1.nreal[i1]      = ppath2.nreal[i];
      gridpos_copy( ppath1.gp_p[i1], ppath2.gp_p[i] );

      if( ppath1.dim >= 2 )
        { gridpos_copy( ppath1.gp_lat[i1], ppath2.gp_lat[i] ); }
      
      if( ppath1.dim == 3 )
        {
          ppath1.pos(i1,2)        = ppath2.pos(i,2);
          ppath1.los(i1,1)        = ppath2.los(i,1);
          gridpos_copy( ppath1.gp_lon[i1], ppath2.gp_lon[i] ); 
        }
      
      ppath1.l_step[i1-1] = ppath2.l_step[i-1];
    }

  if( ppath_what_background( ppath2 ) )
    { ppath1.background = ppath2.background; }

  ppath.lspace += ppath.lspace;
}



//! ppath_fill_1d
/*!
   Fills a 1D Ppath structure with position and LOS values.

   The function fills the fields: pos, los, z, l_step and gp_p.
   The pressure grid positions (gp_p) are filtered through gridpos_check_fd.

   The structure fields must be allocated to correct size before calling the 
   function. The field size must be at least as large as the length of r,
   lat and za vectors.

   The length along the path shall be the same between all points.

   \param   ppath          Output: Ppath structure.
   \param   r              Vector with radius for the path points.
   \param   lat            Vector with latitude for the path points.
   \param   za             Vector with zenith angle for the path points.
   \param   lstep          Length along the path between the points.
   \param   refellipsoid   As the WSV with the same name
   \param   z_field        Geometrical altitudes.
   \param   ip             Pressure grid range.

   \author Patrick Eriksson
   \date   2002-07-18
*/
void ppath_fill_1d(
           Ppath&      ppath,
     ConstVectorView   r,
     ConstVectorView   lat,
     ConstVectorView   za,
     ConstVectorView   lstep,
     ConstVectorView   refellipsoid,
     ConstVectorView   z_field,
     const Index&      ip )
{
  // Help variables that are common for all points.
  const double   r1 = refellipsoid[0] + z_field[ip];
  const double   dr = z_field[ip+1] - z_field[ip];

  for( Index i=0; i<r.nelem(); i++ )
    {
      ppath.r[i]     = r[i];
      ppath.pos(i,0) = r[i] - refellipsoid[0];
      ppath.pos(i,1) = lat[i];
      ppath.los(i,0) = za[i];
      
      ppath.gp_p[i].idx   = ip;
      ppath.gp_p[i].fd[0] = ( r[i] - r1 ) / dr;
      ppath.gp_p[i].fd[1] = 1 - ppath.gp_p[i].fd[0];
      gridpos_check_fd( ppath.gp_p[i] );

      if( i > 0 )
        { ppath.l_step[i-1] = lstep[i-1]; }
    }
}



//! ppath_fill_2d
/*!
   Fills a 2D Ppath structure with position and LOS values.

   The function fills the fields: pos, los, z, l_step, gp_p and gp_lat.

   The structure fields must be allocated to correct size before calling the 
   function. The field size must be at least as large as the length of r,
   lat and za vectors.

   The length along the path shall be the same between all points.

   \param   ppath          Output: Ppath structure.
   \param   r              Vector with radius for the path points.
   \param   lat            Vector with latitude for the path points.
   \param   za             Vector with zenith angle for the path points.
   \param   lstep          Length along the path between the points.
   \param   refellipsoid   As the WSV with the same name.
   \param   z_field        Geometrical altitudes
   \param   lat_grid       Latitude grid.
   \param   ip             Pressure grid range.
   \param   ilat           Latitude grid range.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_fill_2d(
           Ppath&      ppath,
     ConstVectorView   r,
     ConstVectorView   lat,
     ConstVectorView   za,
     const double&     lstep,
     ConstVectorView   refellipsoid,
     ConstMatrixView   z_field,
     ConstVectorView   lat_grid,
     const Index&      ip,
     const Index&      ilat )
{
  // Help variables that are common for all points.
  const double   dlat  = lat_grid[ilat+1] - lat_grid[ilat];
  const double   z1low = z_field(ip,ilat);
  const double   z1upp = z_field(ip+1,ilat);
  const double   dzlow = z_field(ip,ilat+1) -z1low;
  const double   dzupp = z_field(ip+1,ilat+1) - z1upp;
        double   re    = refell2r( refellipsoid, lat_grid[ilat] );
  const double   r1low = re + z1low;
  const double   r1upp = re + z1upp;
                 re    = refell2r( refellipsoid, lat_grid[ilat+1] );
  const double   drlow = re + z_field(ip,ilat+1) - r1low;
  const double   drupp = re + z_field(ip+1,ilat+1) - r1upp;

  for( Index i=0; i<r.nelem(); i++ )
    {
      ppath.r[i]     = r[i];
      ppath.pos(i,1) = lat[i];
      ppath.los(i,0) = za[i];
      
      // Weight in the latitude direction
      double w = ( lat[i] - lat_grid[ilat] ) / dlat;

      // Radius of lower and upper face at present latitude
      const double rlow = r1low + w * drlow;
      const double rupp = r1upp + w * drupp;

      // Geometrical altitude of lower and upper face at present latitude
      const double zlow = z1low + w * dzlow;
      const double zupp = z1upp + w * dzupp;

      ppath.gp_p[i].idx   = ip;
      ppath.gp_p[i].fd[0] = ( r[i] - rlow ) / ( rupp - rlow );
      ppath.gp_p[i].fd[1] = 1 - ppath.gp_p[i].fd[0];
      gridpos_check_fd( ppath.gp_p[i] );

      ppath.pos(i,0) = zlow + ppath.gp_p[i].fd[0] * ( zupp -zlow );

      ppath.gp_lat[i].idx   = ilat;
      ppath.gp_lat[i].fd[0] = ( lat[i] - lat_grid[ilat] ) / dlat;
      ppath.gp_lat[i].fd[1] = 1 - ppath.gp_lat[i].fd[0];
      gridpos_check_fd( ppath.gp_lat[i] );

      if( i > 0 )
        { ppath.l_step[i-1] = lstep; }
    }
}



//! ppath_fill_3d
/*!
   Fills a 3D Ppath structure with position and LOS values.

   The function fills the fields: pos, los, z, l_step, gp_p, gp_lat and gp_lon.

   The structure fields must be allocated to correct size before calling the 
   function. The field size must be at least as large as the length of r,
   lat and za vectors.

   The length along the path shall be the same between all points.

   \param   ppath          Output: Ppath structure.
   \param   r              Vector with radius for the path points.
   \param   lat            Vector with latitude for the path points.
   \param   lon            Vector with longitude for the path points.
   \param   za             Vector with zenith angle for the path points.
   \param   aa             Vector with azimuth angle for the path points.
   \param   lstep          Length along the path between the points.
   \param   refellipsoid   As the WSV with the same name.
   \param   z_field        Geometrical altitudes
   \param   lat_grid       Latitude grid.
   \param   lon_grid       Longitude grid.
   \param   ip             Pressure grid range.
   \param   ilat           Latitude grid range.
   \param   ilon           Longitude grid range.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void ppath_fill_3d(
           Ppath&      ppath,
     ConstVectorView   r,
     ConstVectorView   lat,
     ConstVectorView   lon,
     ConstVectorView   za,
     ConstVectorView   aa,
     const double&     lstep,
     ConstVectorView   refellipsoid,
     ConstTensor3View  z_field,
     ConstVectorView   lat_grid,
     ConstVectorView   lon_grid,
     const Index&      ip,
     const Index&      ilat,
     const Index&      ilon )
{
  // Help variables that are common for all points.
  const double   lat1  = lat_grid[ilat];
  const double   lat3  = lat_grid[ilat+1];
  const double   lon5  = lon_grid[ilon];
  const double   lon6  = lon_grid[ilon+1];
  const double   re1   = refell2r( refellipsoid, lat_grid[ilat] );
  const double   re3   = refell2r( refellipsoid, lat_grid[ilat+1] );
  const double   r15a  = re1 + z_field(ip,ilat,ilon);
  const double   r35a  = re3 + z_field(ip,ilat+1,ilon); 
  const double   r36a  = re3 + z_field(ip,ilat+1,ilon+1); 
  const double   r16a  = re1 + z_field(ip,ilat,ilon+1);
  const double   r15b  = re1 + z_field(ip+1,ilat,ilon);
  const double   r35b  = re3 + z_field(ip+1,ilat+1,ilon); 
  const double   r36b  = re3 + z_field(ip+1,ilat+1,ilon+1);
  const double   r16b  = re1 + z_field(ip+1,ilat,ilon+1);
  const double   dlat  = lat3 - lat1;
  const double   dlon  = lon6 - lon5;

  for( Index i=0; i<r.nelem(); i++ )
    {
      // Radius of pressure levels at present lat and lon
      const double   rlow = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                      r15a, r35a, r36a, r16a, lat[i], lon[i] );
      const double   rupp = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                      r15b, r35b, r36b, r16b, lat[i], lon[i] );

      // Position and LOS
      ppath.r[i]     = r[i];
      ppath.pos(i,1) = lat[i];
      ppath.pos(i,2) = lon[i];
      ppath.los(i,0) = za[i];
      ppath.los(i,1) = aa[i];
      
      // Pressure grid index
      ppath.gp_p[i].idx   = ip;
      ppath.gp_p[i].fd[0] = ( r[i] - rlow ) / ( rupp - rlow );
      ppath.gp_p[i].fd[1] = 1 - ppath.gp_p[i].fd[0];
      gridpos_check_fd( ppath.gp_p[i] );

      // Geometrical altitude
      const double   re = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                           re1, re3, re3, re1, lat[i], lon[i] );
      const double   zlow = rlow - re;
      const double   zupp = rupp - re;
      //
      ppath.pos(i,0) = zlow + ppath.gp_p[i].fd[0] * ( zupp -zlow );

      // Latitude grid index
      ppath.gp_lat[i].idx   = ilat;
      ppath.gp_lat[i].fd[0] = ( lat[i] - lat1 ) / dlat;
      ppath.gp_lat[i].fd[1] = 1 - ppath.gp_lat[i].fd[0];
      gridpos_check_fd( ppath.gp_lat[i] );

      // Longitude grid index
      //
      // The longitude  is undefined at the poles. The grid index is set to
      // the start point.
      //
      if( abs( lat[i] ) < 90 )
        {
          ppath.gp_lon[i].idx   = ilon;
          ppath.gp_lon[i].fd[0] = ( lon[i] - lon5 ) / dlon;
          ppath.gp_lon[i].fd[1] = 1 - ppath.gp_lon[i].fd[0];
          gridpos_check_fd( ppath.gp_lon[i] );
        }
      else
        {
          ppath.gp_lon[i].idx   = 0;
          ppath.gp_lon[i].fd[0] = 0;
          ppath.gp_lon[i].fd[1] = 1;
        }

      if( i > 0 )
        { ppath.l_step[i-1] = lstep; }
    }
}





/*===========================================================================
  === Functions related to propagation paths with refraction
  ===========================================================================*/

//! refraction_ppc
/*! 
   Calculates the propagation path constant for cases with refraction.

   Both positive and negative zenith angles are handled.

   \return               Path constant.
   \param   r            Radius.
   \param   za           LOS Zenith angle.
   \param   refr_index   Refractive index.

   \author Patrick Eriksson
   \date   2002-05-17
*/
double refraction_ppc( 
        const double&   r, 
        const double&   za, 
        const double&   refr_index )
{
  assert( r > 0 );
  assert( abs(za) <= 180 );

  return r * refr_index * sin( DEG2RAD * abs(za) );
}





/*===========================================================================
  === Help functions for the *ppath_step* functions found below
  === These functions are mainly pieces of code that are common for at least
  === two functions (or two places in some function) and for this reason 
  === the headers are not complete. 
  ===========================================================================*/


//! ppath_start_1d
/*! 
   Internal help function for 1D path calculations.

   The function does the asserts and determined some variables that are common
   for geometrical and refraction calculations.

   See the code for details.

   \author Patrick Eriksson
   \date   2002-11-13
*/
void ppath_start_1d(
              double&     r_start,
              double&     lat_start,
              double&     za_start,
              Index&      ip,
        const Ppath&      ppath )
{
  // Number of points in the incoming ppath
  const Index   imax = ppath.np - 1;

  // Extract starting radius, zenith angle and latitude
  r_start   = ppath.r[imax];
  lat_start = ppath.pos(imax,1);
  za_start  = ppath.los(imax,0);

  // Determine index of the pressure level being the lower limit for the
  // grid range of interest.
  //
  ip = gridpos2gridrange( ppath.gp_p[imax], za_start<=90 );
}



//! ppath_end_1d
/*! 
   Internal help function for 1D path calculations.

   The function performs the end part of the calculations, that are common
   for geometrical and refraction calculations.

   See the code for details.

   \author Patrick Eriksson
   \date   2002-11-27
*/
void ppath_end_1d(
              Ppath&      ppath,
        ConstVectorView   r_v,
        ConstVectorView   lat_v,
        ConstVectorView   za_v,
        const double&     lstep,
        ConstVectorView   z_field,
        ConstVectorView   refellipsoid,
        const Index&      ip,
        const Index&      endface,
        const double&     ppc )
{
  // Number of path points
  const Index   np = r_v.nelem();

  // Re-allocate ppath for return results and fill the structure
  //
  ppath_init_structure( ppath, 1, np );
  //
  ppath.constant = ppc;
  //
  ppath_fill_1d( ppath, r_v, lat_v, za_v, Vector(np-1,lstep), refellipsoid, 
                                                                 z_field, ip );

  gridpos_check_fd( ppath.gp_p[np-1] );

  //--- End point is the surface
  if( endface == 7 )
    { ppath_set_background( ppath, 2 ); }

  //--- End point is on top of a pressure level
  else
    { gridpos_force_end_fd( ppath.gp_p[np-1], z_field.nelem() ); }
}



//! ppath_start_2d
/*! 
   Internal help function for 2D path calculations.

   The function does the asserts and determined some variables that are common
   for geometrical and refraction calculations.

   See the code for details.

   \author Patrick Eriksson
   \date   2002-11-18
*/
void ppath_start_2d(
              double&     r_start,
              double&     lat_start,
              double&     za_start,
              Index&      ip,
              Index&      ilat,
              double&     lat1,
              double&     lat3,
              double&     r1a,
              double&     r3a,
              double&     r3b,
              double&     r1b,
              double&     rsurface1,
              double&     rsurface3,
              Ppath&      ppath,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   refellipsoid,
        ConstVectorView   z_surface
        )
{
  // Number of points in the incoming ppath
  const Index imax = ppath.np - 1;

  // Extract starting radius, zenith angle and latitude
  r_start   = ppath.r[imax];
  lat_start = ppath.pos(imax,1);
  za_start  = ppath.los(imax,0);

  // Determine interesting latitude grid range and latitude end points of 
  // the range.
  //
  ilat = gridpos2gridrange( ppath.gp_lat[imax], za_start >= 0 );
  //
  lat1 = lat_grid[ilat];
  lat3 = lat_grid[ilat+1];

  // Determine interesting pressure grid range. Do this first assuming that
  // the pressure levels are not tilted (that is, abs(za_start<=90) always
  // mean upward observation). 
  // Set radius for the corners of the grid cell and the radial slope of 
  // pressure level limits of the grid cell to match the found ip.
  //
  ip = gridpos2gridrange( ppath.gp_p[imax], abs(za_start) <= 90);
  //
  const double re1 = refell2r( refellipsoid, lat_grid[ilat] );
  const double re3 = refell2r( refellipsoid, lat_grid[ilat+1] );
  //
  r1a = re1 + z_field(ip,ilat);        // lower-left
  r3a = re3 + z_field(ip,ilat+1);      // lower-right
  r3b = re3 + z_field(ip+1,ilat+1);    // upper-right
  r1b = re1 + z_field(ip+1,ilat);      // upper-left

  // This part is a fix to catch start postions on top of a pressure level
  // that does not have an end fractional distance for the first step.
  {
    // Radius of lower and upper pressure level at the start position
    const double   rlow = rsurf_at_lat( lat1, lat3, r1a, r3a, lat_start );
    const double   rupp = rsurf_at_lat( lat1, lat3, r1b, r3b, lat_start );
    if( abs(r_start-rlow) < RTOL || abs(r_start-rupp) < RTOL )
      { gridpos_force_end_fd( ppath.gp_p[imax], z_field.nrows() ); }
  }
  
  // Slopes of pressure levels
  double   c2 = plevel_slope_2d( lat1, lat3, r1a, r3a );
  double   c4 = plevel_slope_2d( lat1, lat3, r1b, r3b );

  // Check if the LOS zenith angle happen to be between 90 and the zenith angle
  // of the pressure level (that is, 90 + tilt of pressure level), and in
  // that case if ip must be changed. This check is only needed when the
  // start point is on a pressure level.
  //
  if( is_gridpos_at_index_i( ppath.gp_p[imax], ip )  )
    {
      double tilt = plevel_angletilt( r_start, c2 );

      if( is_los_downwards( za_start, tilt ) )
        {
          ip--;
          r1b = r1a;   r3b = r3a;   c4 = c2;
          r1a = re1 + z_field(ip,ilat);
          r3a = re3 + z_field(ip,ilat+1);
          c2 = plevel_slope_2d( lat1, lat3, r1a, r3a );
        }
    }
  else if( is_gridpos_at_index_i( ppath.gp_p[imax], ip+1 )  )
    {
      double tilt = plevel_angletilt( r_start, c4 );

      if( !is_los_downwards( za_start, tilt ) )
        {
          ip++;
          r1a = r1b;   r3a = r3b;   c2 = c4;
          r3b = re3 + z_field(ip+1,ilat+1);
          r1b = re1 + z_field(ip+1,ilat);    
          c4 = plevel_slope_2d( lat1, lat3, r1b, r3b );
        }
    }

  // Surface radius at latitude end points
  rsurface1 = re1 + z_surface[ilat];
  rsurface3 = re3 + z_surface[ilat+1];
}



//! ppath_end_2d
/*! 
   Internal help function for 2D path calculations.

   The function performs the end part of the calculations, that are common
   for geometrical and refraction calculations.

   See the code for details.

   \author Patrick Eriksson
   \date   2002-11-29
*/
void ppath_end_2d(
              Ppath&      ppath,
        ConstVectorView   r_v,
        ConstVectorView   lat_v,
        ConstVectorView   za_v,
        const double&     lstep,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   refellipsoid,
        const Index&      ip,
        const Index&      ilat,
        const Index&      endface,
        const double&     ppc )
{
  // Number of path points
  const Index   np = r_v.nelem();
  const Index   imax = np-1;

  // Re-allocate ppath for return results and fill the structure
  //
  ppath_init_structure( ppath, 2, np );
  //
  ppath.constant = ppc;
  //
  ppath_fill_2d( ppath, r_v, lat_v, za_v, lstep, refellipsoid, z_field, 
                                                          lat_grid, ip, ilat );

  gridpos_check_fd( ppath.gp_p[imax] );
  gridpos_check_fd( ppath.gp_lat[imax] );

  // Do end-face specific tasks
  if( endface == 7 )
    { ppath_set_background( ppath, 2 ); }

  // Set fractional distance for end point
  //
  if( endface == 1  ||  endface == 3 )
    { gridpos_force_end_fd( ppath.gp_lat[np-1], lat_grid.nelem() ); }
  else if( endface == 2  ||  endface == 4 )
    { gridpos_force_end_fd( ppath.gp_p[np-1], z_field.nrows() ); }

  // Handle cases where exactly a corner is hit, or when slipping outside of
  // the grid box due to numerical inaccuarcy
  if( ppath.gp_p[imax].fd[0] < 0  ||  ppath.gp_p[imax].fd[1] < 0 )
    { 
      gridpos_force_end_fd( ppath.gp_p[imax], z_field.nrows() ); 
    }
  if( ppath.gp_lat[imax].fd[0] < 0  ||  ppath.gp_lat[imax].fd[1] < 0 )
    { 
      gridpos_force_end_fd( ppath.gp_lat[imax], lat_grid.nelem() ); 
    }
}



//! ppath_start_3d
/*! 
   Internal help function for 3D path calculations.

   The function does the asserts and determined some variables that are common
   for geometrical and refraction calculations.

   See the code for details.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void ppath_start_3d(
              double&     r_start,
              double&     lat_start,
              double&     lon_start,
              double&     za_start,
              double&     aa_start,
              Index&      ip,
              Index&      ilat,
              Index&      ilon,
              double&     lat1,
              double&     lat3,
              double&     lon5,
              double&     lon6,
              double&     r15a,
              double&     r35a,
              double&     r36a,
              double&     r16a,
              double&     r15b,
              double&     r35b,
              double&     r36b,
              double&     r16b,
              double&     rsurface15,
              double&     rsurface35,
              double&     rsurface36,
              double&     rsurface16,
              Ppath&      ppath,
              ConstVectorView   lat_grid,
              ConstVectorView   lon_grid,
              ConstTensor3View  z_field,
              ConstVectorView   refellipsoid,
              ConstMatrixView   z_surface )
{
  // Index of last point in the incoming ppath
  const Index imax = ppath.np - 1;

  // Extract starting radius, zenith angle and latitude
  r_start   = ppath.r[imax];
  lat_start = ppath.pos(imax,1);
  lon_start = ppath.pos(imax,2);
  za_start  = ppath.los(imax,0);
  aa_start  = ppath.los(imax,1);

  // Number of lat/lon
  const Index nlat = lat_grid.nelem();
  const Index nlon = lon_grid.nelem();

  // Lower index of lat and lon ranges of interest
  //
  // The longitude is undefined at the poles and as the azimuth angle
  // is defined in other way at the poles.
  //
  if( lat_start == 90 )
    { 
      ilat = nlat - 2;
      GridPos   gp_tmp;
      gridpos( gp_tmp, lon_grid, aa_start );
      if( aa_start < 180 )
        { ilon = gridpos2gridrange( gp_tmp, 1 ); }
      else
        { ilon = gridpos2gridrange( gp_tmp, 0 ); }
    }
  else if( lat_start == -90 )
    { 
      ilat = 0; 
      GridPos   gp_tmp;
      gridpos( gp_tmp, lon_grid, aa_start );
      if( aa_start < 180 )
        { ilon = gridpos2gridrange( gp_tmp, 1 ); }
      else
        { ilon = gridpos2gridrange( gp_tmp, 0 ); }
    }
  else
    { 
      if( lat_start > 0 )
        ilat = gridpos2gridrange( ppath.gp_lat[imax], abs( aa_start ) < 90 );
      else
        ilat = gridpos2gridrange( ppath.gp_lat[imax], abs( aa_start ) <= 90 );
      if( lon_start < lon_grid[nlon-1] )
        { ilon = gridpos2gridrange( ppath.gp_lon[imax], aa_start >= 0 ); }
      else
        { ilon = nlon - 2; }
    }
  //
  lat1 = lat_grid[ilat];
  lat3 = lat_grid[ilat+1];
  lon5 = lon_grid[ilon];
  lon6 = lon_grid[ilon+1];

  // Determine interesting pressure grid range. Do this first assuming that
  // the pressure levels are not tilted (that is, abs(za_start<=90) always
  // mean upward observation). 
  // Set radius for the corners of the grid cell and the radial slope of 
  // pressure level limits of the grid cell to match the found ip.
  //
  ip = gridpos2gridrange( ppath.gp_p[imax], za_start <= 90 );
  //
  const double re1 = refell2r( refellipsoid, lat_grid[ilat] );
  const double re3 = refell2r( refellipsoid, lat_grid[ilat+1] );
  //
  r15a = re1 + z_field(ip,ilat,ilon);
  r35a = re3 + z_field(ip,ilat+1,ilon); 
  r36a = re3 + z_field(ip,ilat+1,ilon+1); 
  r16a = re1 + z_field(ip,ilat,ilon+1);
  r15b = re1 + z_field(ip+1,ilat,ilon);
  r35b = re3 + z_field(ip+1,ilat+1,ilon); 
  r36b = re3 + z_field(ip+1,ilat+1,ilon+1); 
  r16b = re1 + z_field(ip+1,ilat,ilon+1);

  // Check if the LOS zenith angle happen to be between 90 and the zenith angle
  // of the pressure level (that is, 90 + tilt of pressure level), and in
  // that case if ip must be changed. This check is only needed when the
  // start point is on a pressure level.
  //
  if( fabs(za_start-90) <= PTILTMAX )
    {
      if( is_gridpos_at_index_i( ppath.gp_p[imax], ip )  )
        {
          // Slope and angular tilt of lower pressure level
          double c2 = plevel_slope_3d( lat1, lat3, lon5, lon6, 
                       r15a, r35a, r36a, r16a, lat_start, lon_start, aa_start );
          double tilt = plevel_angletilt( r_start, c2 );
	  // when tilt < 1sec assume this is a numeric problem and reset tilt to 0)
	  if( tilt < 2.e-4 )
	    {
	      tilt=0.;
	    }

          if( is_los_downwards( za_start, tilt ) )
            {
              ip--;
              r15b = r15a;   r35b = r35a;   r36b = r36a;   r16b = r16a;
              r15a = re1 + z_field(ip,ilat,ilon);
              r35a = re3 + z_field(ip,ilat+1,ilon); 
              r36a = re3 + z_field(ip,ilat+1,ilon+1); 
              r16a = re1 + z_field(ip,ilat,ilon+1);
            }
        }
      else if( is_gridpos_at_index_i( ppath.gp_p[imax], ip+1 )  )
        {
          // Slope and angular tilt of lower pressure level
          double c4 = plevel_slope_3d( lat1, lat3 ,lon5, lon6, 
                       r15b, r35b, r36b, r16b, lat_start, lon_start, aa_start );
          double tilt = plevel_angletilt( r_start, c4 );

          if( !is_los_downwards( za_start, tilt ) )
            {
              ip++;
              r15a = r15b;   r35a = r35b;   r36a = r36b;   r16a = r16b;
              r15b = re1 + z_field(ip+1,ilat,ilon);
              r35b = re3 + z_field(ip+1,ilat+1,ilon); 
              r36b = re3 + z_field(ip+1,ilat+1,ilon+1); 
              r16b = re1 + z_field(ip+1,ilat,ilon+1);
            }
        }
    }

  // Surface radius at latitude/longitude corner points
  rsurface15 = re1 + z_surface(ilat,ilon);
  rsurface35 = re3 + z_surface(ilat+1,ilon);
  rsurface36 = re3 + z_surface(ilat+1,ilon+1);
  rsurface16 = re1 + z_surface(ilat,ilon+1);
}



//! ppath_end_3d
/*! 
   Internal help function for 3D path calculations.

   The function performs the end part of the calculations, that are common
   for geometrical and refraction calculations.

   See the code for details.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void ppath_end_3d(
              Ppath&      ppath,
        ConstVectorView   r_v,
        ConstVectorView   lat_v,
        ConstVectorView   lon_v,
        ConstVectorView   za_v,
        ConstVectorView   aa_v,
        const double&     lstep,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,
        ConstTensor3View  z_field,
        ConstVectorView   refellipsoid,
        const Index&      ip,
        const Index&      ilat,
        const Index&      ilon,
        const Index&      endface,
        const double&     ppc )
{
  // Number of path points
  const Index   np   = r_v.nelem();
  const Index   imax = np-1;

  // Re-allocate ppath for return results and fill the structure
  //
  ppath_init_structure( ppath, 3, np );
  //
  ppath.constant = ppc;
  //
  ppath_fill_3d( ppath, r_v, lat_v, lon_v, za_v, aa_v, lstep, 
                 refellipsoid, z_field, lat_grid, lon_grid, ip, ilat, ilon );

  // Do end-face specific tasks
  if( endface == 7 )
    { ppath_set_background( ppath, 2 ); }

  // Set fractional distance for end point
  //
  if( endface == 1  ||  endface == 3 )
    { gridpos_force_end_fd( ppath.gp_lat[imax], lat_grid.nelem() ); }
  else if( endface == 2  ||  endface == 4 )
    { gridpos_force_end_fd( ppath.gp_p[imax], z_field.npages() ); }
  else if( endface == 5  ||  endface == 6 )
    { gridpos_force_end_fd( ppath.gp_lon[imax], lon_grid.nelem() ); }

  // Handle cases where exactly a corner is hit, or when slipping outside of
  // the grid box due to numerical inaccuarcy
  if( ppath.gp_p[imax].fd[0] < 0  ||  ppath.gp_p[imax].fd[1] < 0 )
    { 
      gridpos_force_end_fd( ppath.gp_p[imax], z_field.npages() ); 
    }
  if( ppath.gp_lat[imax].fd[0] < 0  ||  ppath.gp_lat[imax].fd[1] < 0 )
    { 
      gridpos_force_end_fd( ppath.gp_lat[imax], lat_grid.nelem() ); 
    }
  if( ppath.gp_lon[imax].fd[0] < 0  ||  ppath.gp_lon[imax].fd[1] < 0 )
    { 
      gridpos_force_end_fd( ppath.gp_lon[imax], lon_grid.nelem() ); 
    }
}



//! interpolate_raytracing_points
/*! 
   Interpolates a set of ray tracing points to a set of points linearly
   spaced along the path.

   All quantities are interpolated linearly.

   Empty vectors can be sent as input for *lon_rt* and *aa_rt* for 1D and 2D.
   The output vectors *lon_v* and *aa_v* are then not filled.

   \author Patrick Eriksson
   \date   2002-11-27
*/
void interpolate_raytracing_points(
             Vector&     r_v,
             Vector&     lat_v,
             Vector&     lon_v,
             Vector&     za_v,
             Vector&     aa_v,
             double&     lstep,
       ConstVectorView   r_rt,
       ConstVectorView   lat_rt,
       ConstVectorView   lon_rt,
       ConstVectorView   za_rt,
       ConstVectorView   aa_rt,
       ConstVectorView   l_rt,
       const double&     lmax )
{
  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. If *lmax* <= 0, then only the end points
  // shall be kept.
  //
  const Index     nrt = r_rt.nelem();
        Index     n = 2;
        double   ltotsum = l_rt.sum();
  // Ensure that there are at least two points
  if( lmax > 0 )
    { n = max (Index( ceil( ltotsum / lmax ) ) + 1, Index(2)); }
  //
  r_v.resize(n);
  lat_v.resize(n);
  za_v.resize(n);
  //
  if( n == 2 )
    {
      r_v[0] = r_rt[0];     r_v[1] = r_rt[nrt-1];
      za_v[0] = za_rt[0];   za_v[1] = za_rt[nrt-1];
      lat_v[0] = lat_rt[0]; lat_v[1] = lat_rt[nrt-1];
      lstep = ltotsum;
      if( lon_rt.nelem() > 0 )
        {
            lon_v.resize(n); lon_v[0] = lon_rt[0]; lon_v[1] = lon_rt[nrt-1];
            aa_v.resize(n);  aa_v[0] = aa_rt[0];   aa_v[1] = aa_rt[nrt-1];
        }
    }
  else
    {
      Vector           lsum(nrt), llin(n);
      ArrayOfGridPos   gp(n);
      Matrix           itw(n,2);
      
      lsum[0] = 0;
      for( Index i=1; i<nrt; i++ )
        { lsum[i] = lsum[i-1] + l_rt[i-1]; }

      nlinspace( llin, 0, ltotsum, n );

      gridpos( gp, lsum, llin );

      interpweights( itw, gp );

      interp( r_v, itw, r_rt, gp );
      interp( za_v, itw, za_rt, gp );
      interp( lat_v, itw, lat_rt, gp );
      lstep = llin[1] - llin[0];

      if( lon_rt.nelem() > 0 )
        {
            lon_v.resize(n); interp( lon_v, itw, lon_rt, gp );
            aa_v.resize(n);  interp( aa_v, itw, aa_rt, gp );
        }
    }
}



//! from_raytracingarrays_to_ppath_vectors_1d_and_2d
/*! 
   A small help function to convert arrays with ray tracing points t
   interpolated values along the path.

   This function is common for 1D and 2D.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void from_raytracingarrays_to_ppath_vectors_1d_and_2d(
             Vector&          r_v,
             Vector&          lat_v,
             Vector&          za_v,
             double&          lstep,
       const Array<double>&   r_array,
       const Array<double>&   lat_array,
       const Array<double>&   za_array,
       const Array<double>&   l_array,
       const Index&           reversed,
       const double&          lmax )
{
  // Copy arrays to vectors for later interpolation
  //
  // Number of ray tracing points
  const Index   nrt=r_array.nelem();
  //
  Vector    r_rt(nrt), lat_rt(nrt), za_rt(nrt), l_rt(nrt-1);    
  //
  if( !reversed )
    {
      for( Index i=0; i<nrt; i++ )
        {
          r_rt[i]   = r_array[i];
          lat_rt[i] = lat_array[i]; 
          za_rt[i]  = za_array[i];
          if( i < (nrt-1) )
            { l_rt[i]   = l_array[i]; }
        }
    }
  else
    {
      for( Index i=0; i<nrt; i++ )
        {
          const Index i1 = nrt - 1 - i;
          r_rt[i]   = r_array[i1];
          lat_rt[i] = lat_array[0]+ lat_array[nrt-1] - lat_array[i1]; 
          za_rt[i]  = 180 - za_array[i1];
          if( i < (nrt-1) )
            { l_rt[i]   = l_array[i1-1]; }
        }
    }

  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. If *lmax* <= 0, then only the end points
  //
  Vector    dummy;   // Dummy vector for output variables not used
  //
  interpolate_raytracing_points( r_v, lat_v, dummy, za_v, dummy, lstep,
      r_rt, lat_rt, Vector(0), za_rt, Vector(0), l_rt, lmax );
}



//! from_raytracingarrays_to_ppath_vectors_3d
/*! 
   A small help function to convert arrays with ray tracing points to
   interpolated values along the path.

   \author Patrick Eriksson
   \date   2003-01-18
*/
void from_raytracingarrays_to_ppath_vectors_3d(
             Vector&          r_v,
             Vector&          lat_v,
             Vector&          lon_v,
             Vector&          za_v,
             Vector&          aa_v,
             double&          lstep,
       const Array<double>&   r_array,
       const Array<double>&   lat_array,
       const Array<double>&   lon_array,
       const Array<double>&   za_array,
       const Array<double>&   aa_array,
       const Array<double>&   l_array,
       const double&          lmax )
{
  // Copy arrays to vectors for later interpolation
  //
  // Number of ray tracing points
  const Index   nrt=r_array.nelem();
  //
  Vector    r_rt(nrt), lat_rt(nrt), lon_rt(nrt);
  Vector    za_rt(nrt), aa_rt(nrt),l_rt(nrt-1);    
  //
  for( Index i=0; i<nrt; i++ )
    {
      r_rt[i]   = r_array[i];
      lat_rt[i] = lat_array[i]; 
      lon_rt[i] = lon_array[i]; 
      za_rt[i]  = za_array[i];
      aa_rt[i]  = aa_array[i];
      if( i < (nrt-1) )
        { l_rt[i]   = l_array[i]; }
    }

  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. If *lmax* <= 0, then only the end points
  //
  interpolate_raytracing_points( r_v, lat_v, lon_v, za_v, aa_v, lstep,
                              r_rt, lat_rt, lon_rt, za_rt, aa_rt, l_rt, lmax );
}



/*===========================================================================
  === Core functions for geometrical *ppath_step* functions
  ===========================================================================*/

//! ppath_step_geom_1d
/*! 
   Calculates 1D geometrical propagation path steps.

   This is the core function to determine 1D propagation path steps by pure
   geometrical calculations. Path points are included for crossings with the
   grids, tangent points and points of intersection with the surface. In
   addition, points are included in the propgation path to ensure that the
   distance along the path between the points does not exceed the selected
   maximum length (lmax). If lmax is <= 0, this means that no length criterion
   shall be applied.

   Note that the input variables are here compressed to only hold data for
   a 1D atmosphere. For example, z_field is z_field(:,0,0).

   For more information read the chapter on propagation paths in AUG.

   \param   ppath             Output: A Ppath structure.
   \param   z_field           Geometrical altitudes corresponding to p_grid.
   \param   refellipsoid      As the WSV with the same name.
   \param   z_surface         Surface altitude.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-05-20
*/
void ppath_step_geom_1d(
              Ppath&      ppath,
        ConstVectorView   z_field,
        ConstVectorView   refellipsoid,
        const double&     z_surface,
        const double&     lmax )
{
  // Starting radius, zenith angle and latitude
  double r_start, lat_start, za_start;

  // Index of the pressure level being the lower limit for the
  // grid range of interest.
  Index ip;

  // Determine the variables defined above, and make asserts of input
  ppath_start_1d( r_start, lat_start, za_start, ip, ppath );

  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  double ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }


  // The path is determined by another function. Determine some variables
  // needed b� that function and call the function.
  //
  // Vars to hold found path points, path step length and coding for end face
  Vector   r_v, lat_v, za_v;
  double   lstep;
  Index    endface;
  //
  do_gridrange_1d( r_v, lat_v, za_v, lstep, endface,
                   r_start, lat_start, za_start, ppc, lmax, 
                   refellipsoid[0]+z_field[ip], refellipsoid[0]+z_field[ip+1], 
                   refellipsoid[0]+z_surface );

  // Fill *ppath*
  //
  ppath_end_1d( ppath, r_v, lat_v, za_v, lstep, z_field, refellipsoid, ip, 
                                                                endface, ppc );
  // nreal is set to 1 in ppath_stepGeometric
}



//! ppath_step_geom_2d
/*! 
   Calculates 2D geometrical propagation path steps.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ppath             Output: A Ppath structure.
   \param   lat_grid          Latitude grid.
   \param   z_field           Geometrical altitudes
   \param   refellipsoid      As the WSV with the same name.
   \param   z_surface         Surface altitudes.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_step_geom_2d(
              Ppath&      ppath,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   refellipsoid,
        ConstVectorView   z_surface,
        const double&     lmax )
{
  // Radius, zenith angle and latitude of start point.
  double   r_start, lat_start, za_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat;

  // Radii and latitudes set by *ppath_start_2d*.
  double   lat1, lat3, r1a, r3a, r3b, r1b, rsurface1, rsurface3;

  // Determine the variables defined above and make all possible asserts
  ppath_start_2d( r_start, lat_start, za_start, ip, ilat, 
                  lat1, lat3, r1a, r3a, r3b, r1b, rsurface1, rsurface3,
                  ppath, lat_grid, z_field, refellipsoid, z_surface );

  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  double ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }

  // Vars to hold found path points, path step length and coding for end face
  Vector   r_v, lat_v, za_v;
  double   lstep;
  Index    endface;

  do_gridcell_2dv2( r_v, lat_v, za_v, lstep, endface,
                  r_start, lat_start, za_start, ppc, lmax, lat1, lat3, 
                                    r1a, r3a, r3b, r1b, rsurface1, rsurface3 );

  // Fill *ppath*
  //
  ppath_end_2d( ppath, r_v, lat_v, za_v, lstep, lat_grid, z_field, refellipsoid,
                                                      ip, ilat, endface, ppc );
  // nreal is set to 1 in ppath_stepGeometric
}



//! ppath_step_geom_3d
/*! 
   Calculates 3D geometrical propagation path steps.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ppath             Output: A Ppath structure.
   \param   lat_grid          Latitude grid.
   \param   lon_grid          Longitude grid.
   \param   z_field           Geometrical altitudes
   \param   refellipsoid      As the WSV with the same name.
   \param   z_surface         Surface altitudes.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void ppath_step_geom_3d(
              Ppath&       ppath,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        ConstTensor3View   z_field,
        ConstVectorView    refellipsoid,
        ConstMatrixView    z_surface,
        const double&      lmax )
{
  // Radius, zenith angle and latitude of start point.
  double   r_start, lat_start, lon_start, za_start, aa_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat, ilon;

  // Radius for corner points, latitude and longitude of the grid cell
  //
  double   lat1, lat3, lon5, lon6;
  double   r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b;
  double   rsurface15, rsurface35, rsurface36, rsurface16;

  // Determine the variables defined above and make all possible asserts
  ppath_start_3d( r_start, lat_start, lon_start, za_start, aa_start, 
                  ip, ilat, ilon, lat1, lat3, lon5, lon6,
                  r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b, 
                  rsurface15, rsurface35, rsurface36, rsurface16,
                  ppath, lat_grid, lon_grid, z_field, refellipsoid, z_surface );


  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  double ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }


  // Vars to hold found path points, path step length and coding for end face
  Vector   r_v, lat_v, lon_v, za_v, aa_v;
  double   lstep;
  Index    endface;

  do_gridcell_3d( r_v, lat_v, lon_v, za_v, aa_v, lstep, endface,
                  r_start, lat_start, lon_start, za_start, aa_start, ppc, lmax,
                  lat1, lat3, lon5, lon6, 
                  r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b,
                  rsurface15, rsurface35, rsurface36, rsurface16 );

  // Fill *ppath*
  //
  ppath_end_3d( ppath, r_v, lat_v, lon_v, za_v, aa_v, lstep, lat_grid, 
                lon_grid, z_field, refellipsoid, ip, ilat, ilon, endface, ppc );
  // nreal is set to 1 in ppath_stepGeometric
}





/*===========================================================================
  === Ray tracing functions
  ===========================================================================*/

//! raytrace_1d_linear_euler
/*! 
   Performs ray tracing for 1D with linear Euler steps.

   A geometrical step with length of *lraytrace* is taken from each
   point. The zenith angle for the end point of that step is
   calculated exactly by the expression c = r*n*sin(theta), and a new
   step is taken. The length of the last ray tracing step to reach the
   end radius is adopted to the distance to the end radius.

   The refractive index is assumed to vary linearly between the pressure
   levels.

   As the ray tracing is performed from the last end point, the found path
   will not be symmetric around the tangent point.

   For more information read the chapter on propagation paths in AUG.
   The algorithm used is described in that part of AUG.

   The array variables *r_array*, *lat_array* and *za_array* shall include
   the start position when calling the function. The length of *l_array*
   will be one smaller than the length of the other arrays.

   \param   ws           Current Workspace
   \param   r_array      Out: Radius of ray tracing points.
   \param   lat_array    Out: Latitude of ray tracing points.
   \param   za_array     Out: LOS zenith angle at ray tracing points.
   \param   l_array      Out: Distance along the path between ray tracing 
                         points.
   \param   endface      FIXME: Add documentation.
   \param   r            Start radius for ray tracing.
   \param   lat          Start latitude for ray tracing.
   \param   za           Start zenith angle for ray tracing.
   \param   rte_pressure     FIXME: Add documentation.
   \param   rte_temperature  FIXME: Add documentation.
   \param   rte_vmr_list     FIXME: Add documentation.
   \param   refr_index   Refractive index.
   \param   refr_index_agenda   FIXME: Add documentation.
   \param   ppc          Propagation path constant.
   \param   lraytrace    Maximum allowed length for ray tracing steps.
   \param   r1           Radius of lower pressure level.
   \param   r3           Radius of upper pressure level (r3 > r1).
   \param   r_surface     Radius of the surface.
   \param   refellipsoid  As the WSV with the same name.
   \param   p_grid       Pressure grid.
   \param   z_field      Geometrical altitudes corresponding to p_grid.
   \param   t_field      Temperatures corresponding to p_grid.
   \param   vmr_field    VMR values corresponding to p_grid.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void raytrace_1d_linear_euler(
              Workspace&       ws,
              Array<double>&   r_array,
              Array<double>&   lat_array,
              Array<double>&   za_array,
              Array<double>&   l_array,
              Index&           endface,
              double           r,
              double           lat,
              double           za,
              Numeric&         rte_pressure,
              Numeric&         rte_temperature,
              Vector&          rte_vmr_list,
              Numeric&         refr_index,
        const Agenda&          refr_index_agenda,
        const double&          ppc,
        const double&          lraytrace,
        const double&          r1,
        const double&          r3,
        const double&          r_surface,
        ConstVectorView        refellipsoid,
        ConstVectorView        p_grid,
        ConstVectorView        z_field,
        ConstVectorView        t_field,
        ConstMatrixView        vmr_field )
{
  // Loop boolean
  bool ready = false;

  // Variables for output from do_gridrange_1d
  Vector    r_v, lat_v, za_v;
  double   lstep, dlat = 9999;

  while( !ready )
    {
      // Constant for the geometrical step to make
      const double   ppc_step = geometrical_ppc( r, za );

      // Where will a geometric path exit the grid cell?
      do_gridrange_1d( r_v, lat_v, za_v, lstep, endface, r, lat, za,
                                              ppc_step, -1, r1, r3, r_surface );

      assert( r_v.nelem() == 2 );

      // If *lstep* is < *lraytrace*, extract the found end point and
      // we are ready. 
      // Otherwise, we make a geometrical step with length *lraytrace*.
      // To avoid that the refrcation calculations end up outside the grid
      // cell when *lstep* is just below *lraytrace*, we allow *lstep*
      // to be somewhat larger than *lraytrace*.

      if( lstep <= 1.01 * lraytrace )
        {
          r     = r_v[1];
          dlat  = lat_v[1] - lat;
          lat   = lat_v[1];
          ready = true;
        }
      else
        {
          if( abs(za) <= 90 )
            { lstep = lraytrace; }
          else
            { lstep = -lraytrace; }

          const double   r_new = geompath_r_at_l( ppc_step, 
                                         geompath_l_at_r(ppc_step,r) + lstep );

          dlat = RAD2DEG * acos( ( r_new*r_new + r*r - 
                                             lstep*lstep ) / ( 2 * r_new*r ) );
          r     = r_new;
          lat   = lat + dlat;
          lstep = abs( lstep );
        }

      // Calculate LOS zenith angle at found point.
      //
      // Refractive index at *r*
      get_refr_index_1d( ws, refr_index, rte_pressure, rte_temperature,
                         rte_vmr_list, refr_index_agenda,
                         p_grid, refellipsoid, z_field, t_field, vmr_field, r );

      const double   ppc_local = ppc / refr_index; 

      if( r >= ppc_local )
        { za = geompath_za_at_r( ppc_local, za, r ); }
      else
        {  
          r = ppc_local;
          za = 90;
        }
  
      // Store found point
      r_array.push_back( r );
      lat_array.push_back( lat );
      za_array.push_back( za );
      l_array.push_back( lstep );
    }  
}



//! raytrace_2d_linear_euler
/*! 
   Performs ray tracing for 2D with linear Euler steps.

   A geometrical step with length of *lraytrace* is taken from each
   point. The zenith angle for the end point of that step is
   calculated considering the gradient of the refractive index. The
   length of the last ray tracing step to reach the end radius is
   adopted to the distance to the end radius.

   The refractive index is assumed to vary linearly along the pressure
   levels and the latitude grid points.

   For more information read the chapter on propagation paths in AUG.
   The algorithm used is described in that part of AUG.

   The array variables *r_array*, *lat_array* and *za_array* shall include
   the start position when calling the function. The length of *l_array*
   will be one smaller than the length of the other arrays.

   \param   ws              Current Workspace
   \param   r_array         Out: Radius of ray tracing points.
   \param   lat_array       Out: Latitude of ray tracing points.
   \param   za_array        Out: LOS zenith angle at ray tracing points.
   \param   l_array         Out: Distance along the path between ray tracing 
                            points.
   \param   endface         Out: Number coding of exit face.
   \param   r               Out: Start radius for ray tracing.
   \param   lat             Out: Start latitude for ray tracing.
   \param   za              Out: Start zenith angle for ray tracing.
   \param   rte_pressure      Out: The WSV with the same name.
   \param   rte_temperature   Out: The WSV with the same name.
   \param   rte_vmr_list      Out: The WSV with the same name.
   \param   refr_index        Out: The WSV with the same name.
   \param   refr_index_agenda    The WSV with the same name.
   \param   lraytrace       Maximum allowed length for ray tracing steps.
   \param   lat1            Latitude of left end face of the grid cell.
   \param   lat3            Latitude of right end face  of the grid cell.
   \param   r1a             Radius of lower-left corner of the grid cell.
   \param   r3a             Radius of lower-right corner of the grid cell.
   \param   r3b             Radius of upper-right corner of the grid cell.
   \param   r1b             Radius of upper-left corner of the grid cell.
   \param   rsurface1        Radius for the surface at *lat1*.
   \param   rsurface3        Radius for the surface at *lat3*.
   \param   p_grid          The WSV with the same name.
   \param   lat_grid        The WSV with the same name.
   \param   refellipsoid    The WSV with the same name.
   \param   z_field         The WSV with the same name.
   \param   t_field         The WSV with the same name.
   \param   vmr_field       The WSV with the same name.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void raytrace_2d_linear_euler(
              Workspace&       ws,
              Array<double>&   r_array,
              Array<double>&   lat_array,
              Array<double>&   za_array,
              Array<double>&   l_array,
              Index&           endface,
              double           r,
              double           lat,
              double           za,
              Numeric&         rte_pressure,
              Numeric&         rte_temperature,
              Vector&          rte_vmr_list,
              Numeric&         refr_index,
        const Agenda&          refr_index_agenda,
        const double&          lraytrace,
        const double&          lat1,
        const double&          lat3,
        const double&          r1a,
        const double&          r3a,
        const double&          r3b,
        const double&          r1b,
        const double&          rsurface1,
        const double&          rsurface3,
        ConstVectorView        p_grid,
        ConstVectorView        lat_grid,
        ConstVectorView        refellipsoid,
        ConstMatrixView        z_field,
        ConstMatrixView        t_field,
        ConstTensor3View       vmr_field )
{
  // Loop boolean
  bool ready = false;

  // Variables for output from do_gridcell_2d
  Vector   r_v, lat_v, za_v;
  double   lstep, dlat = 9999, r_new, lat_new;

  while( !ready )
    {
      // Constant for the geometrical step to make
      const double   ppc_step = geometrical_ppc( r, za );

      // Where will a geometric path exit the grid cell?
      do_gridcell_2d( r_v, lat_v, za_v, lstep, endface,
                      r, lat, za, ppc_step, -1, lat1, lat3, r1a, r3a, r3b, r1b, 
                      rsurface1, rsurface3 );

      assert( r_v.nelem() == 2 );

      // If *lstep* is < *lraytrace*, extract the found end point and
      // we are ready.
      // Otherwise, we make a geometrical step with length *lraytrace*.
      // To avoid that the refrcation calculations end up outside the grid
      // cell when *lstep* is just below *lraytrace*, we allow *lstep*
      // to be somewhat larger than *lraytrace*.

      if( lstep <= 1.01 * lraytrace )
        {
          r_new   = r_v[1];
          dlat    = lat_v[1] - lat;
          lat_new = lat_v[1];
          ready   = true;
        }
      else
        {
          if( abs(za) <= 90 )
            { lstep = lraytrace; }
          else
            { lstep = -lraytrace; }

          r_new = geompath_r_at_l( ppc_step, 
                                         geompath_l_at_r(ppc_step,r) + lstep );
          dlat = RAD2DEG * acos( ( r_new*r_new + r*r - 
                                             lstep*lstep ) / ( 2 * r_new*r ) );
          if( za < 0 )
            { dlat = -dlat; }

          lat_new = lat + dlat;
          lstep   = abs( lstep );

          // For paths along the latitude end faces we can end up outside the
          // grid cell. We simply look for points outisde the grid cell.
          if( lat_new < lat1 )
            { lat_new = lat1; }
          else if( lat_new > lat3 )
            { lat_new = lat3; }
        }

      // Calculate LOS zenith angle at found point.
        {
                Numeric   dndr, dndlat;
          const double   za_rad = DEG2RAD * za;

          refr_gradients_2d( ws, refr_index, dndr, dndlat, rte_pressure, 
                             rte_temperature, rte_vmr_list, refr_index_agenda,
                             p_grid, lat_grid, refellipsoid, z_field, t_field, 
                             vmr_field, r, lat );

          za += -dlat + RAD2DEG * lstep / refr_index * ( -sin(za_rad) * dndr +
                                                        cos(za_rad) * dndlat );

          // Make sure that obtained *za* is inside valid range
          if( za < -180 )
            { za += 360; }
          else if( za > 180 )
            { za -= 360; }
        }

      r   = r_new;
      lat = lat_new;

      // If the path is zenith/nadir along a latitude end face, we must check
      // that the path does not exit with new *za*.
      if( lat == lat1  &&  za < 0 )
        { endface = 1;   ready = 1; }
      else if( lat == lat3  &&  za > 0 )
        { endface = 3;   ready = 1; }
      
      // Store found point
      r_array.push_back( r );
      lat_array.push_back( lat );
      za_array.push_back( za );
      l_array.push_back( lstep );
    }  
}



//! raytrace_3d_linear_euler
/*! 
   Performs ray tracing for 3D with linear Euler steps.

   A geometrical step with length of *lraytrace* is taken from each
   point. The zenith angle for the end point of that step is
   calculated considering the gradient of the refractive index. The
   length of the last ray tracing step to reach the end radius is
   adopted to the distance to the end radius.

   The refractive index is assumed to vary linearly along the pressure
   levels and the latitude grid points.

   For more information read the chapter on propagation paths in AUG.
   The algorithm used is described in that part of AUG.

   The array variables *r_array*, *lat_array*, *lon_array*, *za_array*
   and *aa_array* shall include the start position when calling the
   function. The length of *l_array* will be one smaller than the
   length of the other arrays.

   \param   ws             Current Workspace
   \param   r_array        Out: Radius of ray tracing points.
   \param   lat_array      Out: Latitude of ray tracing points.
   \param   lon_array      Out: Longitude of ray tracing points.
   \param   za_array       Out: LOS zenith angle at ray tracing points.
   \param   aa_array       Out: LOS azimuth angle at ray tracing points.
   \param   l_array        Out: Distance along the path between ray tracing 
                           points.
   \param   endface        Out: Number coding of exit face.
   \param   r              Out: Start radius for ray tracing.
   \param   lat            Out: Start latitude for ray tracing.
   \param   lon            Out: Start longitude for ray tracing.
   \param   za             Out: Start zenith angle for ray tracing.
   \param   aa             Out: Start azimuth angle for ray tracing.
   \param   rte_pressure     Out: The WSV with the same name.
   \param   rte_temperature  Out: The WSV with the same name.
   \param   rte_vmr_list     Out: The WSV with the same name.
   \param   refr_index       Out: The WSV with the same name.
   \param   refr_index_agenda    The WSV with the same name.
   \param   lraytrace      Maximum allowed length for ray tracing steps.
   \param   lat1           Latitude of left end face of the grid cell.
   \param   lat3           Latitude of right end face of the grid cell.
   \param   lon5           Lower longitude of the grid cell.
   \param   lon6           Upper longitude of the grid cell.
   \param   r15a           Radius of corner: lower p-level,*lat1* and *lon5*.
   \param   r35a           Radius of corner: lower p-level,*lat3* and *lon5*.
   \param   r36a           Radius of corner: lower p-level,*lat3* and *lon6*.
   \param   r16a           Radius of corner: lower p-level,*lat1* and *lon6*.
   \param   r15b           Radius of corner: upper p-level,*lat1* and *lon5*.
   \param   r35b           Radius of corner: upper p-level,*lat3* and *lon5*.
   \param   r36b           Radius of corner: upper p-level,*lat3* and *lon6*.
   \param   r16b           Radius of corner: upper p-level,*lat1* and *lon6*.
   \param   rsurface15      Radius for the surface at *lat1* and *lon5*.
   \param   rsurface35      Radius for the surface at *lat3* and *lon5*.
   \param   rsurface36      Radius for the surface at *lat3* and *lon6*.
   \param   rsurface16      Radius for the surface at *lat1* and *lon6*.   
   \param   p_grid         The WSV with the same name.
   \param   lat_grid       The WSV with the same name.
   \param   lon_grid       The WSV with the same name.
   \param   refellipsoid   The WSV with the same name.
   \param   z_field        The WSV with the same name.
   \param   t_field        The WSV with the same name.
   \param   vmr_field      The WSV with the same name.

   \author Patrick Eriksson
   \date   2003-01-18
*/
void raytrace_3d_linear_euler(
              Workspace&       ws,
              Array<double>&   r_array,
              Array<double>&   lat_array,
              Array<double>&   lon_array,
              Array<double>&   za_array,
              Array<double>&   aa_array,
              Array<double>&   l_array,
              Index&           endface,
              double           r,
              double           lat,
              double           lon,
              double           za,
              double           aa,
              Numeric&         rte_pressure,
              Numeric&         rte_temperature,
              Vector&          rte_vmr_list,
              Numeric&         refr_index,
        const Agenda&          refr_index_agenda,
        const double&          lraytrace,
        const double&          lat1,
        const double&          lat3,
        const double&          lon5,
        const double&          lon6,
        const double&          r15a,
        const double&          r35a,
        const double&          r36a,
        const double&          r16a,
        const double&          r15b,
        const double&          r35b,
        const double&          r36b,
        const double&          r16b,
        const double&          rsurface15,
        const double&          rsurface35,
        const double&          rsurface36,
        const double&          rsurface16,
        ConstVectorView        p_grid,
        ConstVectorView        lat_grid,
        ConstVectorView        lon_grid,
        ConstVectorView        refellipsoid,
        ConstTensor3View       z_field,
        ConstTensor3View       t_field,
        ConstTensor4View       vmr_field )
{
  // Loop boolean
  bool ready = false;

  // Variables for output from do_gridcell_2d
  Vector    r_v, lat_v, lon_v, za_v, aa_v;
  double    r_new, lat_new, lon_new, za_new, aa_new;
  double    lstep;

  while( !ready )
    {
      // Constant for the geometrical step to make
      const double   ppc_step = geometrical_ppc( r, za );

      // Where will a geometric path exit the grid cell?
      do_gridcell_3d( r_v, lat_v, lon_v, za_v, aa_v, lstep, endface,
                    r, lat, lon, za, aa, ppc_step, -1, lat1, lat3, lon5, lon6, 
                    r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b,
                      rsurface15, rsurface35, rsurface36, rsurface16 );

      assert( r_v.nelem() == 2 );

      // If *lstep* is < *lraytrace*, extract the found end point and
      // we are ready.
      // Otherwise, we make a geometrical step with length *lraytrace*.
      // To avoid that the refrcation calculations end up outside the grid
      // cell when *lstep* is just below *lraytrace*, we allow *lstep*
      // to be somewhat larger than *lraytrace*.

      if( lstep <= 1.01 * lraytrace )
        {
          r_new   = r_v[1];
          lat_new = lat_v[1];
          lon_new = lon_v[1];
          za_new  = za_v[1];
          aa_new  = aa_v[1];
          ready   = true;
        }
      else
        {
          // Sensor pos and LOS in cartesian coordinates
          double   x, y, z, dx, dy, dz;
          poslos2cart( x, y, z, dx, dy, dz, r, lat, lon, za, aa ); 

          lstep = lraytrace;

          cart2poslos( r_new, lat_new, lon_new, za_new, aa_new, 
                              x+dx*lstep, y+dy*lstep, z+dz*lstep, dx, dy, dz );

          // For paths along some end face we can end up outside the
          // grid cell. We simply look for points outisde the grid cell.
          if( lat_new < lat1 )
            { lat_new = lat1; }
          else if( lat_new > lat3 )
            { lat_new = lat3; }
          if( lon_new < lon5 )
            { lon_new = lon5; }
          else if( lon_new > lon6 )
            { lon_new = lon6; }

          // Checks to improve the accuracy for speciel cases.
          // The checks are only needed for values cacluated locally as
          // the same checks are made inside *do_gridcell_3d*.

          //--- Set zenith angle to be as accurate as possible
          if( za < ANGTOL  ||  za > 180-ANGTOL )
            { za_new = za; }
          else
            { za_new = geompath_za_at_r( ppc_step, za, r_new ); }

          //--- Set azimuth angle and lon. to be as accurate as possible for
          //    north-south observations
          if( abs( lat ) < 90  &&  
                     ( abs(aa) < ANGTOL  ||  abs( aa ) > 180-ANGTOL ) )
            { 
              aa_new  = aa; 
              lon_new = lon;
            }

          // Shall lon values be shifted?
          resolve_lon( lon_new, lon5, lon6 );
        }

      // Calculate LOS zenith angle at found point.
      {
        Numeric   dndr, dndlat, dndlon;

        refr_gradients_3d( ws,
                             refr_index, dndr, dndlat, dndlon, rte_pressure, 
                             rte_temperature, rte_vmr_list, refr_index_agenda, 
                             p_grid, lat_grid, lon_grid, refellipsoid, z_field, 
                             t_field, vmr_field, r, lat, lon );

        const double   aterm = RAD2DEG * lstep / refr_index;
        const double   za_rad = DEG2RAD * za;
        const double   aa_rad = DEG2RAD * aa;
        const double   sinza = sin( za_rad );
        const double   sinaa = sin( aa_rad );
        const double   cosaa = cos( aa_rad );

        if( za < ANGTOL  ||  za > 180-ANGTOL )
          { 
            za_new += aterm * ( cos(za_rad) * 
                                ( cosaa * dndlat + sinaa * dndlon ) );
            aa_new = RAD2DEG * atan2( dndlon, dndlat); 
          }
        else
          { 
            za_new += aterm * ( -sinza * dndr + cos(za_rad) * 
                                       ( cosaa * dndlat + sinaa * dndlon ) );
            aa_new += aterm * sinza * ( cosaa * dndlon - sinaa * dndlat ); 
          }
          
        // Make sure that obtained angles are inside valid ranges
        if( za_new > 180 )
          { za_new = 180 - za_new; }
        else if( za_new < 0 )
          { za_new = -za_new; }
        //
        if( aa_new > 180 )
          { aa_new -= 360; }
        else if( aa_new < -180 )
          { aa_new += 360; }
      }

      r   = r_new;
      lat = lat_new;
      lon = lon_new;
      za  = za_new;
      aa  = aa_new;

      // For some cases where the path goes along an end face, 
      // it could be the case that the refraction bends the path out 
      // of the grid cell.
      if( lon == lon5  &&  aa < 0 )
        { endface = 5;   ready = 1; }
      else if( lon == lon6  &&  aa > 0 )
        { endface = 6;   ready = 1; }
      else if( za > 0  &&  za < 180 )
        {
          if( lat == lat1  &&  lat != -90  &&  abs( aa ) > 90 ) 
            { endface = 1;   ready = 1; }
          else if( lat == lat3  &&  lat != 90  &&  abs( aa ) < 90 ) 
            { endface = 3;   ready = 1; }
        }

      // Store found point
      r_array.push_back( r );
      lat_array.push_back( lat );
      lon_array.push_back( lon );
      za_array.push_back( za );
      aa_array.push_back( aa );
      l_array.push_back( lstep );
    }  
}





/*===========================================================================
  === Core functions for refraction *ppath_step* functions
  ===========================================================================*/

//! ppath_step_refr_1d
/*! 
   Calculates 1D propagation path steps including effects of refraction.

   This function works as the function *ppath_step_geom_1d* but considers
   also refraction. The upper length of the ray tracing steps is set by
   the argument *lraytrace*. This argument controls only the internal
   calculations. The maximum distance between the path points is still
   determined by *lmax*.

   \param   ws                Current Workspace
   \param   ppath             Out: A Ppath structure.
   \param   rte_pressure      Out: The WSV with the same name.
   \param   rte_temperature   Out: The WSV with the same name.
   \param   rte_vmr_list      Out: The WSV with the same name.
   \param   refr_index        Out: The WSV with the same name.
   \param   refr_index_agenda The WSV with the same name.
   \param   p_grid            Pressure grid.
   \param   z_field           Geometrical altitudes corresponding to p_grid.
   \param   t_field           Temperatures corresponding to p_grid.
   \param   vmr_field         VMR values corresponding to p_grid.
   \param   refellipsoid      As the WSV wit the same name.
   \param   z_surface         Surface altitude.
   \param   rtrace_method     String giving which ray tracing method to use.
                              See the function for options.
   \param   lraytrace         Maximum allowed length for ray tracing steps.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-11-26
*/
void ppath_step_refr_1d(
              Workspace&  ws,
              Ppath&      ppath,
              Numeric&    rte_pressure,
              Numeric&    rte_temperature,
              Vector&     rte_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        ConstMatrixView   vmr_field,
        ConstVectorView   refellipsoid,
        const double&     z_surface,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax )
{
  // Starting radius, zenith angle and latitude
  double   r_start, lat_start, za_start;

  // Index of the pressure level being the lower limit for the
  // grid range of interest.
  Index   ip;

  // Determine the variables defined above, and make asserts of input
  ppath_start_1d( r_start, lat_start, za_start, ip, ppath );

  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  // If the sensor is placed outside the atmosphere, the constant is
  // already set.
  double ppc;
  if( ppath.constant < 0 )
    { 
      get_refr_index_1d( ws, refr_index, rte_pressure, rte_temperature,
                         rte_vmr_list, refr_index_agenda, p_grid, refellipsoid,
                         z_field, t_field, vmr_field, r_start );
      ppc = refraction_ppc( r_start, za_start, refr_index ); 
    }
  else
    { ppc = ppath.constant; }


  // Perform the ray tracing
  //
  // Arrays to store found ray tracing points
  // (Vectors don't work here as we don't know how many points there will be)
  Array<double>   r_array, lat_array, za_array, l_array;
  //
  // Store the start point
  r_array.push_back( r_start );
  lat_array.push_back( lat_start );
  za_array.push_back( za_start );
  //
  // Number coding for end face
  Index   endface;
  //
  if( rtrace_method  == "linear_euler" )
    {
      raytrace_1d_linear_euler( ws, r_array, lat_array, za_array, l_array, 
             endface, r_start, lat_start, za_start, rte_pressure, 
             rte_temperature, rte_vmr_list, refr_index, refr_index_agenda, 
             ppc, lraytrace, refellipsoid[0]+z_field[ip], 
             refellipsoid[0]+z_field[ip+1], refellipsoid[0] + z_surface, 
             refellipsoid, p_grid, z_field, t_field, vmr_field );
    }
#ifndef NDEBUG
  else
    {
      bool   known_ray_trace_method = false;
      assert( known_ray_trace_method );
    }
#endif

  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. 
  //
  Vector    r_v, lat_v, za_v;
  double   lstep;
  //
  from_raytracingarrays_to_ppath_vectors_1d_and_2d( r_v, lat_v, za_v, lstep, 
                              r_array, lat_array, za_array, l_array, 0, lmax );

  // Fill *ppath*
  ppath_end_1d( ppath, r_v, lat_v, za_v, lstep, z_field, refellipsoid, ip, 
                                                                endface, ppc );
}



//! ppath_step_refr_2d
/*! 
   Calculates 2D propagation path steps, with refraction, using a simple
   and fast ray tracing scheme.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ws                Current Workspace
   \param   ppath             Out: A Ppath structure.
   \param   rte_pressure      Out: The WSV with the same name.
   \param   rte_temperature   Out: The WSV with the same name.
   \param   rte_vmr_list      Out: The WSV with the same name.
   \param   refr_index        Out: The WSV with the same name.
   \param   refr_index_agenda The WSV with the same name.
   \param   p_grid            Pressure grid.
   \param   lat_grid          Latitude grid.
   \param   z_field           Geometrical altitudes.
   \param   t_field           Atmospheric temperatures.
   \param   vmr_field         VMR values.
   \param   refellipsoid      As the WSV with the same name.
   \param   z_surface         Surface altitudes.
   \param   rtrace_method     String giving which ray tracing method to use.
                              See the function for options.
   \param   lraytrace         Maximum allowed length for ray tracing steps.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void ppath_step_refr_2d(
              Workspace&  ws,
              Ppath&      ppath,
              Numeric&    rte_pressure,
              Numeric&    rte_temperature,
              Vector&     rte_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstMatrixView   t_field,
        ConstTensor3View  vmr_field,
        ConstVectorView   refellipsoid,
        ConstVectorView   z_surface,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax )
{
  // Radius, zenith angle and latitude of start point.
  double   r_start, lat_start, za_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat;

  // Radii and latitudes set by *ppath_start_2d*.
  double   lat1, lat3, r1a, r3a, r3b, r1b, rsurface1, rsurface3;

  // Determine the variables defined above and make all possible asserts
  ppath_start_2d( r_start, lat_start, za_start, ip, ilat, 
                  lat1, lat3, r1a, r3a, r3b, r1b, rsurface1, rsurface3,
                  ppath, lat_grid, z_field, refellipsoid, z_surface );

  // Perform the ray tracing
  //
  // No constant for the path is valid here.
  //
  // Arrays to store found ray tracing points
  // (Vectors don't work here as we don't know how many points there will be)
  Array<double>   r_array, lat_array, za_array, l_array;
  //
  // Store the start point
  r_array.push_back( r_start );
  lat_array.push_back( lat_start );
  za_array.push_back( za_start );
  //
  // Number coding for end face
  Index   endface;
  //
  if( rtrace_method  == "linear_euler" )
    {
      raytrace_2d_linear_euler( ws,
        r_array, lat_array, za_array, l_array, endface,
        r_start, lat_start, za_start, rte_pressure, rte_temperature, 
        rte_vmr_list, refr_index, refr_index_agenda, lraytrace, 
        lat1, lat3, r1a, r3a, r3b, r1b, rsurface1, rsurface3,
        p_grid, lat_grid, refellipsoid, z_field, t_field, vmr_field );
    }
#ifndef NDEBUG
  else
    {
      bool   known_ray_trace_method = false;
      assert( known_ray_trace_method );
    }
#endif

  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. 
  //
  Vector    r_v, lat_v, za_v;
  double   lstep;
  //
  from_raytracingarrays_to_ppath_vectors_1d_and_2d( r_v, lat_v, za_v, lstep, 
                              r_array, lat_array, za_array, l_array, 0, lmax );

  // Fill *ppath*
  ppath_end_2d( ppath, r_v, lat_v, za_v, lstep, lat_grid, z_field, refellipsoid,
                                             ip, ilat, endface, -1 );
}



//! ppath_step_refr_3d
/*! 
   Calculates 3D propagation path steps, with refraction, using a simple
   and fast ray tracing scheme.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ws                Current Workspace
   \param   ppath             Out: A Ppath structure.
   \param   rte_pressure      Out: The WSV with the same name.
   \param   rte_temperature   Out: The WSV with the same name.
   \param   rte_vmr_list      Out: The WSV with the same name.
   \param   refr_index        Out: The WSV with the same name.
   \param   refr_index_agenda The WSV with the same name.
   \param   p_grid            Pressure grid.
   \param   lat_grid          Latitude grid.
   \param   lon_grid          Longitude grid.
   \param   z_field           Geometrical altitudes.
   \param   t_field           Atmospheric temperatures.
   \param   vmr_field         VMR values.
   \param   refellipsoid      As teh WSv with the same name.
   \param   z_surface         Surface altitudes.
   \param   rtrace_method     String giving which ray tracing method to use.
                              See the function for options.
   \param   lraytrace         Maximum allowed length for ray tracing steps.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2003-01-08
*/
void ppath_step_refr_3d(
              Workspace&  ws,
              Ppath&      ppath,
              Numeric&    rte_pressure,
              Numeric&    rte_temperature,
              Vector&     rte_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,
        ConstTensor3View  z_field,
        ConstTensor3View  t_field,
        ConstTensor4View  vmr_field,
        ConstVectorView   refellipsoid,
        ConstMatrixView   z_surface,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax )
{
  // Radius, zenith angle and latitude of start point.
  double   r_start, lat_start, lon_start, za_start, aa_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat, ilon;

  // Radius for corner points, latitude and longitude of the grid cell
  //
  double   lat1, lat3, lon5, lon6;
  double   r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b;
  double   rsurface15, rsurface35, rsurface36, rsurface16;

  // Determine the variables defined above and make all possible asserts
  ppath_start_3d( r_start, lat_start, lon_start, za_start, aa_start, 
                  ip, ilat, ilon, lat1, lat3, lon5, lon6,
                  r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b, 
                  rsurface15, rsurface35, rsurface36, rsurface16,
                  ppath, lat_grid, lon_grid, z_field, refellipsoid, z_surface );

  // Perform the ray tracing
  //
  // No constant for the path is valid here.
  //
  // Arrays to store found ray tracing points
  // (Vectors don't work here as we don't know how many points there will be)
  Array<double>   r_array, lat_array, lon_array, za_array, aa_array, l_array;
  //
  // Store the start point
  r_array.push_back( r_start );
  lat_array.push_back( lat_start );
  lon_array.push_back( lon_start );
  za_array.push_back( za_start );
  aa_array.push_back( aa_start );
  //
  // Number coding for end face
  Index   endface;
  //
  if( rtrace_method  == "linear_euler" )
    {
      raytrace_3d_linear_euler( ws, r_array, lat_array, lon_array, za_array, 
                                aa_array, l_array, endface, r_start,
                                lat_start, lon_start, za_start, aa_start, 
                                rte_pressure, rte_temperature, rte_vmr_list, 
                                refr_index, refr_index_agenda, lraytrace, 
                                lat1, lat3, lon5, lon6, 
                                r15a, r35a, r36a, r15a, r15b, r35b, r36b, r15b,
                                rsurface15, rsurface35, rsurface36, rsurface16,
                                p_grid, lat_grid, lon_grid, 
                                refellipsoid, z_field, t_field, vmr_field );
    }
#ifndef NDEBUG
  else
    {
      bool   known_ray_trace_method = false;
      assert( known_ray_trace_method );
    }
#endif

  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. 
  //
  Vector    r_v, lat_v, lon_v, za_v, aa_v;
  double   lstep;
  //
  from_raytracingarrays_to_ppath_vectors_3d( r_v, lat_v, lon_v, za_v, aa_v, 
                                        lstep, r_array, lat_array, lon_array, 
                                           za_array, aa_array, l_array, lmax );

  // Fill *ppath*
  ppath_end_3d( ppath, r_v, lat_v, lon_v, za_v, aa_v, lstep, lat_grid, lon_grid,
                z_field, refellipsoid, ip, ilat, ilon, endface, -1 );
}




/*===========================================================================
  === Core of ppathCalc and help function(s)
  ===========================================================================*/

//! ppath_start_stepping
/*!
   Initiates a Ppath structure for calculation of a path with *ppath_step*.

   The function performs two main tasks. As mentioned above, it initiates
   a Ppath structure (a), but it also checks that the end point of the path is
   at an allowed location (b).

   (a): The Ppath structure is set to hold the position and LOS of the last
   point of the path inside the atmosphere. This point is either the
   sensor position, or the point where the path leaves the model atmosphere.
   If the path is totally outside the atmosphere, no point is put into the
   structure. If the (practical) end and start points are identical, such
   as when the sensor is inside the cloud box, the background field is set.

   (b): If it is found that the end point of the path is at an illegal position
   a detailed error message is given. Not allowed cases are: <br>  
      1. The sensor is placed below surface level. <br> 
      2. For 2D and 3D, the path leaves the model atmosphere at a latitude or
         longitude end face. <br> 
      3. For 2D and 3D, the path is totally outside the atmosphere and the 
         latitude and longitude of the tangent point is outside the range of
         the corresponding grids. 

   All input variables are identical with the WSV with the same name.
   The output variable is here called ppath for simplicity, but is in
   fact *ppath_step*.

   \param   ppath             Output: A Ppath structure.
   \param   atmosphere_dim    The atmospheric dimensionality.
   \param   p_grid            The pressure grid.
   \param   lat_grid          The latitude grid.
   \param   lon_grid          The longitude grid.
   \param   z_field           The field of geometrical altitudes.
   \param   refellipsoid      As the WSV with the same name.
   \param   z_surface         Surface altitude.
   \param   cloudbox_on       Flag to activate the cloud box.
   \param   cloudbox_limits   Index limits of the cloud box.
   \param   outside_cloudbox  Boolean to flag if the propagation path is 
                              (expected to be) outside the cloudbox. Ordinary
                              clerar sky calculations are selected by the value
                              1. The value 0 means tracking of a propagation 
                              path inside the cloudbox. The path is then
                              tracked to the cloudbox boundary.
   \param   rte_pos           The position of the sensor.
   \param   rte_los           The line-of-sight of the sensor.

   \author Patrick Eriksson
   \date   2002-05-17
*/
void ppath_start_stepping(
          Ppath&          ppath,
    const Index&          atmosphere_dim,
    ConstVectorView       p_grid,
    ConstVectorView       lat_grid,
    ConstVectorView       lon_grid,
    ConstTensor3View      z_field,
    ConstVectorView       refellipsoid,
    ConstMatrixView       z_surface,
    const Index&          cloudbox_on, 
    const ArrayOfIndex&   cloudbox_limits,
    const bool&           outside_cloudbox,
    ConstVectorView       rte_pos,
    ConstVectorView       rte_los,
    const Verbosity&      verbosity)
{
  CREATE_OUT1
  CREATE_OUT2
  
  // This function contains no checks or asserts as it is only a sub-function.

  // Allocate the ppath structure
  ppath_init_structure(  ppath, atmosphere_dim, 1 );

  // Index of last pressure level
  const Index lp = p_grid.nelem() - 1;

  // The different atmospheric dimensionalities are handled seperately

  //-- 1D ---------------------------------------------------------------------
  if( atmosphere_dim == 1 )
    {
      // Sensor is inside the model atmosphere:
      if( rte_pos[0] < z_field(lp,0,0) )
        {
          // Check that the sensor is above the surface
          if( (rte_pos[0] + RTOL) < z_surface(0,0) )
            {
              ostringstream os;
              os << "The ppath starting point is placed " 
                 << (z_surface(0,0) - rte_pos[0])/1e3 
                 << " km below the surface.";
              throw runtime_error(os.str());
            }

          // Set ppath
          ppath.pos(0,0) = rte_pos[0];
          ppath.pos(0,1) = 0; 
          ppath.r[0]     = refellipsoid[0] + rte_pos[0];
          ppath.los(0,0) = rte_los[0];
          //
          gridpos( ppath.gp_p, z_field(joker,0,0), ppath.pos(0,0) );
          gridpos_check_fd( ppath.gp_p[0] );

          // Is the sensor on the surface looking down?
          // If yes and the sensor is inside the cloudbox, the background will
          // be changed below.
          if( ppath.pos(0,0) <= z_surface(0,0)  &&  ppath.los(0,0) > 90 )
            { ppath_set_background( ppath, 2 ); }

          // Check sensor position with respect to cloud box.
          if( cloudbox_on )
            {
              const double fgp = fractional_gp( ppath.gp_p[0] );
              //
              // Is the sensor inside the cloud box?
              if( fgp > cloudbox_limits[0]  && fgp < cloudbox_limits[1] )
                { 
                  if( outside_cloudbox )
                    { ppath_set_background( ppath, 4 ); }
                }
              // Is the sensor at the cloud box boundary?
              else if( ( fgp==cloudbox_limits[0] && rte_los[0] <= 90 ) || 
                       ( fgp==cloudbox_limits[1] && rte_los[0] > 90 ) )
                { ppath_set_background( ppath, 3 ); }
            }
          else
            { assert( outside_cloudbox ); }
        }

      // Sensor is outside the model atmosphere:
      else
        {
          // We can here set ppc as we are outside the atmosphere and n=1
          ppath.constant = geometrical_ppc( refellipsoid[0] + rte_pos[0], 
                                                              rte_los[0] );

          // Totally outside
          if( rte_los[0] <= 90  ||  
              ppath.constant >= refellipsoid[0] + z_field(lp,0,0) )
            {
              ppath.pos(0,0) = rte_pos[0];
              ppath.pos(0,1) = 0; 
              ppath.r[0]     = refellipsoid[0] + rte_pos[0];
              ppath.los(0,0) = rte_los[0];
              //
              ppath_set_background( ppath, 1 );
              out1 << "  --- WARNING ---, path is totally outside of the "
                   << "model atmosphere\n";
            }

          // Path enters the atmosphere. 
          else
            { 
              ppath.r[0]     = refellipsoid[0] + z_field(lp,0,0);
              ppath.pos(0,0) = z_field(lp,0,0);
              ppath.los(0,0) = geompath_za_at_r( ppath.constant, rte_los[0],
                                                                 ppath.r[0] );
              ppath.pos(0,1) = geompath_lat_at_za( rte_los[0], 0, 
                                                             ppath.los(0,0) );
              ppath.lspace = geompath_l_at_r( ppath.constant,
                                              refellipsoid[0] + rte_pos[0] ) -
                             geompath_l_at_r( ppath.constant, ppath.r[0] ); 

              // Here we know the grid position exactly
              ppath.gp_p[0].idx   = lp-1; 
              ppath.gp_p[0].fd[0] = 1;
              ppath.gp_p[0].fd[1] = 0;

              // If cloud box reaching TOA, we have also found the background
              if( cloudbox_on  &&  cloudbox_limits[1] == lp )
                { ppath_set_background( ppath, 3 ); }
           }
        }
    }  // End 1D


  //-- 2D ---------------------------------------------------------------------
  else if( atmosphere_dim == 2 )
    {
      // Index of last latitude
      const Index llat = lat_grid.nelem() -1;

      // Is sensor inside range of lat_grid?
      // If yes, determine TOA altitude at sensor position
      GridPos   gp_lat;
      Vector    itw(2);
      bool      islatin = false;
      double    r_e;  // Ellipsoid radius at sensor position
      double    z_toa   = -99e99;
      if( rte_pos[1] > lat_grid[0]  &&  rte_pos[1] < lat_grid[llat] )
        { 
          islatin = true; 
          gridpos( gp_lat, lat_grid, rte_pos[1] );
          interpweights( itw, gp_lat );
          z_toa = interp( itw, z_field(lp,joker,0), gp_lat );
          r_e   = refell2d( refellipsoid, lat_grid, gp_lat );
        }
      else
        { r_e = refell2r( refellipsoid, rte_pos[1] ); }

      // Sensor is inside the model atmosphere:
      if( islatin  &&  rte_pos[0] < z_toa )
        {
          const double z_s = interp( itw, z_surface(joker,0), gp_lat );

          // Check that the sensor is above the surface
          if( (rte_pos[0] + RTOL) < z_s )
            {
              ostringstream os;
              os << "The ppath starting point is placed " 
                 << (z_s - rte_pos[0])/1e3 << " km below the surface.";
              throw runtime_error(os.str());
            }

          // Set ppath
          ppath.pos(0,0) = rte_pos[0];
          ppath.pos(0,1) = rte_pos[1];
          ppath.r[0]     = r_e + rte_pos[0];
          ppath.los(0,0) = rte_los[0];

          // Create a vector with the geometrical altitude of the pressure 
          // levels for the sensor latitude and use it to get ppath.gp_p.
          Vector z_grid(lp+1);
          z_at_lat_2d( z_grid, p_grid, lat_grid, z_field(joker,joker,0), 
                                                                      gp_lat );
          gridpos( ppath.gp_p, z_grid, rte_pos[0] );
          gridpos_check_fd( ppath.gp_p[0] );

          // Latitude grid position
          gridpos_copy( ppath.gp_lat[0], gp_lat );
          gridpos_check_fd( ppath.gp_lat[0] );

          // Is the sensor on the surface looking down?
          // If yes and the sensor is inside the cloudbox, the background will
          // be changed below.
          if( ppath.pos(0,0) <= z_s )
            { 
              // Calculate radial slope of the surface
              const double rslope = plevel_slope_2d( lat_grid, refellipsoid,
                                 z_surface(joker,0), gp_lat, ppath.los(0,0) );

              // Calculate angular tilt of the surface
              const double atilt = plevel_angletilt( r_e + z_s, rslope );

              // Are we looking down into the surface?
              // If yes and the sensor is inside the cloudbox, the background 
              // will be changed below.
              if( is_los_downwards( ppath.los(0,0), atilt ) )
                { ppath_set_background( ppath, 2 ); }
            }

          // Check sensor position with respect to cloud box.
          if( cloudbox_on )
            {
              const double fgp = fractional_gp( ppath.gp_p[0] );
              const double fgl = fractional_gp( ppath.gp_lat[0] );
              //
              // Is the sensor inside the cloud box?
              if( fgp > cloudbox_limits[0]  &&  fgp < cloudbox_limits[1]  &&
                  fgl > cloudbox_limits[2]  &&  fgl < cloudbox_limits[3]  )
                {
                  if( outside_cloudbox )
                    { ppath_set_background( ppath, 4 ); }
                }
              // Is the sensor at the cloud box boundary?
              else if( ( fgp==cloudbox_limits[0] && abs(rte_los[0]) <= 90 ) ||
                       ( fgp==cloudbox_limits[1] && abs(rte_los[0]) > 90 )  ||
                       ( fgl==cloudbox_limits[2] && rte_los[0] >= 0 ) || 
                       ( fgl==cloudbox_limits[3] && rte_los[0] <= 0 ) )
                { ppath_set_background( ppath, 3 ); }
              else
                { assert( outside_cloudbox ); }
            }
        }      

      // Sensor is outside the model atmosphere:
      else
        {
          // Handle cases when the sensor looks in the wrong way
          if( ( rte_pos[1] <= lat_grid[0]     &&  rte_los[0] <= 0 )  || 
              ( rte_pos[1] >= lat_grid[llat]  &&  rte_los[0] >= 0 ) )
            {
              ostringstream os;
              os << "The sensor is outside (or at the limit) of the model "
                 << "atmosphere but\nlooks in the wrong direction (wrong sign "
                 << "for the zenith angle?).\nThis case includes nadir "
                 << "looking exactly at the latitude end points.";
              throw runtime_error( os.str() );
            }

          // We can here set the path constant, that equals the radius of
          // the geometrical tangent point.
          const Numeric r_p = r_e + rte_pos[0];
          ppath.constant = geometrical_ppc( r_p, rte_los[0] );

          // Determine TOA radii, and min and max value
          Vector   r_toa(llat+1);
          double   r_toa_min = 99e99, r_toa_max = -1;
          for( Index ilat=0; ilat<=llat; ilat++ )
            {
              r_toa[ilat] = refell2r( refellipsoid, lat_grid[ilat] ) +
                                                            z_field(lp,ilat,0);
              if( r_toa[ilat] < r_toa_min )
                { r_toa_min = r_toa[ilat]; } 
              if( r_toa[ilat] > r_toa_max )
                { r_toa_max = r_toa[ilat]; } 
            }
          if( r_p <= r_toa_max )
            {
              ostringstream os;
              os << "The sensor is horizontally outside (or at the limit) of "
                 << "the model\natmosphere, but is at a radius smaller than "
                 << "the maximum value of\nthe top-of-the-atmosphere radii. "
                 << "This is not allowed. Make the\nmodel atmosphere larger "
                 << "to also cover the sensor position?";
              throw runtime_error( os.str() );
            }

          // Upward:
          if( abs(rte_los[0]) <= 90 )
            {
              ppath.pos(0,0) = rte_pos[0];
              ppath.pos(0,1) = rte_pos[1]; 
              ppath.r[0]     = r_e + rte_pos[0];
              ppath.los(0,0) = rte_los[0];
              //
              ppath_set_background( ppath, 1 );
              out1 << "  ------- WARNING -------: path is totally outside of "
                   << "the model atmosphere\n";
            }

          // Downward:
          else
            {
              bool      above=false, ready=false, failed=false;
              double    rt=-1, latt, lt, lt_old = L_NOT_FOUND;
              GridPos   gp_latt;
              Vector    itwt(2);

              // Check if clearly above the model atmosphere
              if( ppath.constant >= r_toa_max )
                { above=true; ready=true; }
              else
                { // Otherwise pick out a suitable first test radius
                  if( islatin  ||  ppath.constant > r_toa_min ) 
                    { rt = r_toa_max; } 
                  else 
                    { rt = r_toa_min; }
                }

              // Iterate until solution found or moving out from model atm.
              //
              while( !ready && !failed )
                {
                  // If rt < ppath.constant, ppath above atmosphere
                  if( rt < ppath.constant )
                    { 
                      above = true; 
                      ready = true; 
                    }

                  else
                    {
                      // Calculate length to entrance point at rt 
                      r_crossing_2d( latt, lt, rt, r_p, rte_pos[1], rte_los[0],
                                                              ppath.constant );
                      assert( lt < 9e9 );

                      // Entrance outside range of lat_grid = fail
                      if( latt < lat_grid[0]  ||  latt > lat_grid[llat] )
                        { failed = true; }

                      // OK iteration
                      else
                        {
                          // Converged? (Is length to entrance point found 
                          // with an accuracy of about 1 mm?)
                          if( abs( lt-lt_old ) < 1e-3 )
                            { ready = true; }

                          // Update rt
                          lt_old = lt;
                          gridpos( gp_latt, lat_grid, latt );
                          interpweights( itwt, gp_latt );
                          rt = interp( itwt, r_toa, gp_latt );
                        }
                    }
                }  // while

              if( failed )
                {
                  ostringstream os;
                  os << "The path does not enter the model atmosphere. It "
                     << "reaches the\ntop of the atmosphere "
                     << "altitude around latitude " << latt << " deg.";
                  throw runtime_error( os.str() );
                }
              else if( above )
                {
                  ppath.pos(0,0) = rte_pos[0];
                  ppath.pos(0,1) = rte_pos[1]; 
                  ppath.r[0]     = r_e + rte_pos[0];
                  ppath.los(0,0) = rte_los[0];
                  //
                  ppath_set_background( ppath, 1 );
                  out1 << "  ------- WARNING -------: path is totally outside "
                       << "of the model atmosphere\n";
                }
              else
                {
                  ppath.r[0]     = rt;
                  ppath.pos(0,0) = interp( itwt, z_field(lp,joker,0), gp_latt );
                  // Calculate za first and use to determine lat
                  ppath.los(0,0) = geompath_za_at_r( ppath.constant, 
                                                     rte_los[0], rt );
                  ppath.pos(0,1) = geompath_lat_at_za( rte_los[0], rte_pos[1], 
                                                              ppath.los(0,0) );
                  ppath.lspace = lt;

                  // Here we know the pressure grid position exactly
                  ppath.gp_p[0].idx   = lp-1; 
                  ppath.gp_p[0].fd[0] = 1;
                  ppath.gp_p[0].fd[1] = 0;

                  // Latitude grid position already calculated
                  gridpos_copy( ppath.gp_lat[0], gp_latt );

                  // Hit with cloudbox reaching TOA?
                  if( cloudbox_on  &&  cloudbox_limits[1] == lp )
                    {
                      Numeric fgp = fractional_gp(gp_latt);
                      if( fgp >= (Numeric)cloudbox_limits[2]  &&
                          fgp <= (Numeric)cloudbox_limits[3] )
                        { ppath_set_background( ppath, 3 ); }
                    }
                }
            }  // Downward 
        }  // Outside
    }  // End 2D


  //-- 3D ---------------------------------------------------------------------
  else
    {
      // Index of last latitude and longitude
      const Index llat = lat_grid.nelem() -1;
      const Index llon = lon_grid.nelem() -1;

      // Is sensor inside range of lat_grid and lon_grid?
      // If yes, determine TOA altitude at sensor position
      GridPos   gp_lat, gp_lon;
      Vector    itw(4);
      bool      islatlonin = false;
      double    r_e;  // Ellipsoid radius at sensor position
      double    z_toa   = -99e99;
      if( rte_pos[1] > lat_grid[0]  &&  rte_pos[1] < lat_grid[llat]  &&
          rte_pos[2] > lon_grid[0]  &&  rte_pos[2] < lon_grid[llon] )
        { 
          islatlonin = true; 
          gridpos( gp_lat, lat_grid, rte_pos[1] );
          gridpos( gp_lon, lon_grid, rte_pos[2] );
          interpweights( itw, gp_lat, gp_lon );
          z_toa = interp( itw, z_field(lp,joker,joker), gp_lat, gp_lon );
          r_e   = refell2d( refellipsoid, lat_grid, gp_lat );
        }
      else
        { r_e = refell2r( refellipsoid, rte_pos[1] ); }

      // Sensor is inside the model atmosphere:
      if( islatlonin  &&  rte_pos[0] < z_toa )
        {
          const double z_s = interp( itw, z_surface, gp_lat, gp_lon );

          // Check that the sensor is above the surface
          if( (rte_pos[0] + RTOL) < z_s )
            {
              ostringstream os;
              os << "The ppath starting point is placed " 
                 << (z_s - rte_pos[0])/1e3 << " km below the surface.";
              throw runtime_error(os.str());
            }

          // Set ppath
          ppath.pos(0,0) = rte_pos[0];
          ppath.pos(0,1) = rte_pos[1];
          ppath.pos(0,2) = rte_pos[2];
          ppath.r[0]     = r_e + rte_pos[0];
          ppath.los(0,0) = rte_los[0];
          ppath.los(0,1) = rte_los[1];

          // Create a vector with the geometrical altitude of the pressure 
          // levels for the sensor latitude and use it to get ppath.gp_p.
          Vector z_grid(lp+1);
          z_at_latlon( z_grid, p_grid, lat_grid, lon_grid, z_field, gp_lat,
                                                                    gp_lon );
          gridpos( ppath.gp_p, z_grid, rte_pos[0] );
          gridpos_check_fd( ppath.gp_p[0] );

          // Latitude and longitude grid positions
          gridpos_copy( ppath.gp_lat[0], gp_lat );
          gridpos_copy( ppath.gp_lon[0], gp_lon );
          gridpos_check_fd( ppath.gp_lat[0] );
          gridpos_check_fd( ppath.gp_lon[0] );

          // Is the sensor on the surface looking down?
          // If yes and the sensor is inside the cloudbox, the background will
          // be changed below.
          if( ppath.pos(0,0) <= z_s )
            { 
              // Calculate radial slope of the surface
              const double rslope = plevel_slope_3d( lat_grid, lon_grid,
                    refellipsoid, z_surface, gp_lat, gp_lon, ppath.los(0,1) );

              // Calculate angular tilt of the surface
              const double atilt = plevel_angletilt( r_e + z_s, rslope );

              // Are we looking down into the surface?
              // If yes and the sensor is inside the cloudbox, the background 
              // will be changed below.
              if( is_los_downwards( ppath.los(0,0), atilt ) )
                { ppath_set_background( ppath, 2 ); }
            }

          // Check sensor position with respect to cloud box.
          if( cloudbox_on )
            {
              const double fgp = fractional_gp( ppath.gp_p[0] );
              const double fgl = fractional_gp( ppath.gp_lat[0] );
              const double fgo = fractional_gp( ppath.gp_lon[0] );
              //
              // Is the sensor inside the cloud box?
              if( fgp > cloudbox_limits[0]  &&  fgp < cloudbox_limits[1]  &&
                  fgl > cloudbox_limits[2]  &&  fgl < cloudbox_limits[3]  &&
                  fgo > cloudbox_limits[4]  &&  fgo < cloudbox_limits[5]  )
                {
                  if( outside_cloudbox )
                    { ppath_set_background( ppath, 4 ); }
                }
              // Is the sensor at the cloud box boundary?
              else if( ( fgp==cloudbox_limits[0] && rte_los[0] <= 90 ) ||
                       ( fgp==cloudbox_limits[1] && rte_los[0] > 90 )  ||
                       ( fgl==cloudbox_limits[2] && abs(rte_los[1]) <= 90 ) || 
                       ( fgl==cloudbox_limits[3] && abs(rte_los[1]) >= 90 ) ||
                       ( fgo==cloudbox_limits[4] && abs(rte_los[1]) >= 0 ) || 
                       ( fgo==cloudbox_limits[5] && abs(rte_los[1]) <= 0 ) )
                { ppath_set_background( ppath, 3 ); }
              else
                { assert( outside_cloudbox ); }
            }
        }      

      // Sensor is outside the model atmosphere:
      else
        {
          // Handle cases when the sensor appears to look the wrong way in
          // the north-south direction
          if( ( rte_pos[1] <= lat_grid[0]     &&  abs( rte_los[1] ) >= 90 )  || 
              ( rte_pos[1] >= lat_grid[llat]  &&  abs( rte_los[1] ) <= 90 ) )
            {
              ostringstream os;
              os << "The sensor is north or south (or at the limit) of the "
                 << "model atmosphere\nbut looks in the wrong direction.";
              throw runtime_error( os.str() );
            }

          // Handle cases when the sensor appears to look the wrong way in
          // the west-east direction. We demand that the sensor is inside the
          // range of lon_grid even if all longitudes are covered.
          if( ( rte_pos[2] <= lon_grid[0]     &&  rte_los[1] < 0 )  || 
              ( rte_pos[2] >= lon_grid[llon]  &&  rte_los[1] > 0 ) )
            {
              ostringstream os;
              os << "The sensor is east or west (or at the limit) of the "
                 << "model atmosphere\nbut looks in the wrong direction.";
              throw runtime_error( os.str() );
            }

          // We can here set the path constant, that equals the radius of
          // the geometrical tangent point.
          const Numeric r_p = r_e + rte_pos[0];
          ppath.constant = geometrical_ppc( r_p, rte_los[0] );

          // Determine TOA radii, and min and max value
          Matrix   r_toa(llat+1,llon+1);
          double   r_toa_min = 99e99, r_toa_max = -1;
          for( Index ilat=0; ilat<=llat; ilat++ )
            {
              const Numeric r_lat = refell2r(refellipsoid,lat_grid[ilat]);
              for( Index ilon=0; ilon<=llon; ilon++ )
                {
                  r_toa(ilat,ilon) = r_lat+ z_field(lp,ilat,ilon);
                  if( r_toa(ilat,ilon) < r_toa_min )
                    { r_toa_min = r_toa(ilat,ilon); } 
                  if( r_toa(ilat,ilon) > r_toa_max )
                    { r_toa_max = r_toa(ilat,ilon); } 
                }
            }

          if( r_p <= r_toa_max )
            {
              ostringstream os;
              os << "The sensor is horizontally outside (or at the limit) of "
                 << "the model\natmosphere, but is at a radius smaller than "
                 << "the maximum value of\nthe top-of-the-atmosphere radii. "
                 << "This is not allowed. Make the\nmodel atmosphere larger "
                 << "to also cover the sensor position?";
              throw runtime_error( os.str() );
            }

          // Upward:
          if( rte_los[0] <= 90 )
            {
              ppath.pos(0,0) = rte_pos[0];
              ppath.pos(0,1) = rte_pos[1]; 
              ppath.pos(0,1) = rte_pos[2]; 
              ppath.r[0]     = r_e + rte_pos[0];
              ppath.los(0,0) = rte_los[0];
              ppath.los(0,1) = rte_los[1];
              //
              ppath_set_background( ppath, 1 );
              out1 << "  ------- WARNING -------: path is totally outside of "
                   << "the model atmosphere\n";
            }

          // Downward:
          else
            {
              bool      above=false, ready=false, failed=false;
              double    rt=-1, latt, lont, lt, lt_old = L_NOT_FOUND;
              GridPos   gp_latt, gp_lont;
              Vector    itwt(4);

              // Check if clearly above the model atmosphere
              if( ppath.constant >= r_toa_max )
                { above=true; ready=true; }
              else
                { // Otherwise pick out a suitable first test radius
                  if( islatlonin  ||  ppath.constant > r_toa_min ) 
                    { rt = r_toa_max; } 
                  else 
                    { rt = r_toa_min; }
                }
              
              // Sensor pos and LOS in cartesian coordinates
              double   x, y, z, dx, dy, dz;
              poslos2cart( x, y, z, dx, dy, dz, r_p, rte_pos[1], rte_pos[2], 
                                                     rte_los[0], rte_los[1] );

              // Iterate until solution found or moving out from model atm.
              //
              while( !ready && !failed )
                {
                  // If rt < ppath.constant, ppath above atmosphere
                  if( rt < ppath.constant )
                    { 
                      above = true; 
                      ready = true; 
                    }

                  else
                    {
                      // Calculate lat and lon for entrance point at rt 
                      r_crossing_3d( latt, lont, lt, rt, r_p, rte_pos[1], 
                                     rte_pos[2], rte_los[0], ppath.constant,
                                     x, y, z, dx, dy, dz );
                      assert( lt < 9e9 );
                      resolve_lon( lont, lon_grid[0], lon_grid[llon] ); 

                      // Entrance outside range of lat/lon_grids = fail
                      if( latt < lat_grid[0]  ||  latt > lat_grid[llat]  ||
                          lont < lon_grid[0]  ||  lont > lon_grid[llon] )
                        { failed = true; }

                      // OK iteration
                      else
                        {
                          // Converged? (Is length to entrance point found 
                          // with an accuracy of about 1 mm?)
                          if( abs( lt-lt_old ) < 1e-3 )
                            { ready = true; }
                          
                          // Update rt
                          lt_old = lt;
                          gridpos( gp_latt, lat_grid, latt );
                          gridpos( gp_lont, lon_grid, lont );
                          interpweights( itwt, gp_latt, gp_lont );
                          rt = interp( itwt, r_toa, gp_latt, gp_lont );
                        }
                    }
                }  // while

              if( failed )
                {
                  ostringstream os;
                  os << "The path does not enter the model atmosphere. It\n"
                     << "reaches the top of the atmosphere altitude around:\n"
                     << "  lat: " << latt << " deg.\n  lon: " << lont 
                     << " deg.";
                  throw runtime_error( os.str() );
                }
              else if( above )
                {
                  ppath.pos(0,0) = rte_pos[0];
                  ppath.pos(0,1) = rte_pos[1]; 
                  ppath.pos(0,1) = rte_pos[2]; 
                  ppath.r[0]     = r_e + rte_pos[0];
                  ppath.los(0,0) = rte_los[0];
                  ppath.los(0,1) = rte_los[1];
                  //
                  ppath_set_background( ppath, 1 );
                  out1 << "  ------- WARNING -------: path is totally outside "
                       << "of the model atmosphere\n";
                }
              else
                {
                  // Calculate lt for last rt, and use it to determine pos/los
                  lt = geompath_l_at_r( ppath.constant, r_p ) -
                       geompath_l_at_r( ppath.constant, rt );
                  cart2poslos( ppath.r[0], ppath.pos(0,1), ppath.pos(0,2),
                                           ppath.los(0,0), ppath.los(0,1),
                               x+dx*lt, y+dy*lt, z+dz*lt, dx, dy, dz );
                  assert( abs( ppath.r[0] -rt ) < RTOL );
                  resolve_lon( ppath.pos(0,2), lon_grid[0], lon_grid[llon] ); 
                  //
                  ppath.pos(0,0) = interp( itwt, z_field(lp,joker,joker), 
                                                            gp_latt, gp_lont );
                  ppath.lspace = lt;
                  
                  // Here we know the pressure grid position exactly
                  ppath.gp_p[0].idx   = lp-1; 
                  ppath.gp_p[0].fd[0] = 1;
                  ppath.gp_p[0].fd[1] = 0;

                  // Lat and lon grid position already calculated
                  gridpos_copy( ppath.gp_lat[0], gp_latt );
                  gridpos_copy( ppath.gp_lon[0], gp_lont );

                  // Hit with cloudbox reaching TOA?
                  if( cloudbox_on  &&  cloudbox_limits[1] == lp )
                    {
                      Numeric fgp1 = fractional_gp(gp_latt);
                      Numeric fgp2 = fractional_gp(gp_lont);
                      if( fgp1 >= (Numeric)cloudbox_limits[2]  &&
                          fgp1 <= (Numeric)cloudbox_limits[3]  && 
                          fgp2 >= (Numeric)cloudbox_limits[4]  &&
                          fgp2 <= (Numeric)cloudbox_limits[5])
                        { ppath_set_background( ppath, 3 ); }
                    }
                }
            }  // Downward 
        }  // Outside
    }  // End 3D
}





//! ppath_calc
/*! 
   This is the core for the WSM ppathCalc.

   This function takes the same input as ppathCalc (that is, those
   input arguments are the WSV with the same name), but there are some
   additional argument(s):

   \param ws                 Current Workspace
   \param ppath              Output: A Ppath structure
   \param ppath_step_agenda  FIXME: Add documentation.
   \param atmosphere_dim     The atmospheric dimensionality.
   \param p_grid             The pressure grid.
   \param lat_grid           The latitude grid.
   \param lon_grid           The longitude grid.
   \param z_field            The field of geometrical altitudes.
   \param refellipsoid       As the WSM with the same name.
   \param z_surface          Surface altitude.
   \param cloudbox_on        Flag to activate the cloud box.
   \param cloudbox_limits    Index limits of the cloud box.
   \param rte_pos            The position of the sensor.
   \param rte_los            The line-of-sight of the sensor.
   \param outside_cloudbox   Boolean to flag if the propagation path is 
                             (expected to be) outside the cloudbox. Ordinary
                             clear sky calculations are selected by the value
                             1. The value 0 means tracking of a propagation 
                             path inside the cloudbox. The path is then
                             tracked to the cloudbox boundary.

   \author Patrick Eriksson
   \date   2003-01-08
*/
void ppath_calc(Workspace&            ws,
                // WS Output:
                Ppath&                ppath,
                // WS Input:
                const Agenda&         ppath_step_agenda,
                const Index&          atmosphere_dim,
                const Vector&         p_grid,
                const Vector&         lat_grid,
                const Vector&         lon_grid,
                const Tensor3&        z_field,
                const Vector&         refellipsoid,
                const Matrix&         z_surface,
                const Index&          cloudbox_on, 
                const ArrayOfIndex&   cloudbox_limits,
                const Vector&         rte_pos,
                const Vector&         rte_los,
                const bool&           outside_cloudbox,
                const Verbosity&      verbosity)
{
  CREATE_OUT2
  
  // This function is a WSM but it is normally only called from RteCalc. 
  // For that reason, this function does not repeat input checks that are
  // performed in yCalc, it only performs checks regarding the sensor 
  // position and LOS.

  //--- Check input -----------------------------------------------------------

  chk_vector_length( "rte_pos", rte_pos, atmosphere_dim );
  if( atmosphere_dim == 1 )
    {
      chk_vector_length( "rte_los", rte_los, 1 );
      chk_if_in_range( "sensor zenith angle", rte_los[0], 0., 180. );
    }
  else if( atmosphere_dim == 2 )
    {
      chk_vector_length( "rte_los", rte_los, 1 );
      chk_if_in_range( "sensor zenith angle", rte_los[0], -180., 180. );
    }
  else
    {
      chk_if_in_range( "sensor latitude", rte_pos[1], -90., 90. );
      chk_if_in_range( "sensor longitude", rte_pos[2], -360., 360. );
      chk_vector_length( "rte_los", rte_los, 2 );
      chk_if_in_range( "sensor zenith angle", rte_los[0], 0., 180. );
      chk_if_in_range( "sensor azimuth angle", rte_los[1], -180., 180. );
    }
  assert( outside_cloudbox  ||  cloudbox_on );
  
  //--- End: Check input ------------------------------------------------------


  // Some messages
  out2 << "  -------------------------------------\n";
  out2 << "  sensor altitude        : " << rte_pos[0]/1e3 << " km\n";
  if( atmosphere_dim >= 2 )    
    out2 << "  sensor latitude        : " << rte_pos[1] << "\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor longitude       : " << rte_pos[2] << "\n";
  out2 << "  sensor zenith angle    : " << rte_los[0] << "\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor azimuth angle   : " << rte_los[1] << "\n";


  // Initiate the partial Ppath structure. 
  // The function doing the work sets ppath_step to the point of the path
  // inside the atmosphere closest to the sensor, if the path is at all inside
  // the atmosphere.
  // If the background field is set by the function this flags that there is no
  // path to follow (for example when the sensor is inside the cloud box).
  // The function checks also that the sensor position and LOS give an
  // allowed path.
  //
  Ppath   ppath_step;
  //
  ppath_start_stepping( ppath_step, atmosphere_dim, p_grid, lat_grid, 
                        lon_grid, z_field, refellipsoid, z_surface,
                        cloudbox_on, cloudbox_limits, outside_cloudbox, 
                        rte_pos, rte_los, verbosity );
  // For debugging:
  //Print( ppath_step, 0, verbosity );

  // The only data we need to store from this initial ppath_step is lspace
  const  Numeric lspace = ppath_step.lspace;  

  out2 << "  -------------------------------------\n";

  // Perform propagation path steps until the starting point is found, which
  // is flagged by ppath_step by setting the background field.
  //
  // The results of each step, returned by ppath_step_agenda as a new 
  // ppath_step, are stored as an array of Ppath structures.
  //
  Array<Ppath> ppath_array(0);
  Index   np     = 1;   // Counter for number of points of the path
  Index   istep  = 0;   // Counter for number of steps
  //
  const Index imax_p   = p_grid.nelem() - 1;
  const Index imax_lat = lat_grid.nelem() - 1;
  const Index imax_lon = lon_grid.nelem() - 1;

  while( !ppath_what_background( ppath_step ) )
    {
      // Call ppath_step agenda. 
      // The new path step is added to *ppath_array* last in the while block
      //
      istep++;
      //
      ppath_step_agendaExecute( ws, ppath_step, atmosphere_dim, 
                                lat_grid, lon_grid, z_field, refellipsoid, 
                                z_surface, ppath_step_agenda );
      // For debugging:
      //Print( ppath_step, 0, verbosity );

      // Number of points in returned path step
      const Index n = ppath_step.np;

      // Increase the total number
      np += n - 1;

      if( istep > 10e3 )
        throw runtime_error(
          "10 000 path points have been reached. Is this an infinite loop?" );
      

      //----------------------------------------------------------------------
      //---  Check if some boundary is reached
      //----------------------------------------------------------------------

      //--- Outside cloud box ------------------------------------------------
      if( outside_cloudbox )
        {

          // Check if the top of the atmosphere is reached
          if( is_gridpos_at_index_i( ppath_step.gp_p[n-1], imax_p ) )
            { 
              ppath_set_background( ppath_step, 1 ); 
            }

          // Check that path does not exit at a latitude or longitude end face
          if( atmosphere_dim == 2 )
            {
              // Latitude 
              if( is_gridpos_at_index_i( ppath_step.gp_lat[n-1], 0 ) )
                {
                  ostringstream os;
                  os << "The path exits the atmosphere through the lower "
                     << "latitude end face.\nThe exit point is at an altitude" 
                     << " of " << ppath_step.pos(n-1,0)/1e3 << " km.";
                  throw runtime_error( os.str() );
                }
              if( is_gridpos_at_index_i( ppath_step.gp_lat[n-1], imax_lat ) )
                {
                  ostringstream os;
                  os << "The path exits the atmosphere through the upper "
                     << "latitude end face.\nThe exit point is at an altitude" 
                     << " of " << ppath_step.pos(n-1,0)/1e3 << " km.";
                  throw runtime_error( os.str() );
                }
            }
          if( atmosphere_dim == 3 )
            {
              // Latitude 
              if( lat_grid[0] > -90  && 
                           is_gridpos_at_index_i( ppath_step.gp_lat[n-1], 0 ) )
                {
                  ostringstream os;
                  os << "The path exits the atmosphere through the lower "
                     << "latitude end face.\nThe exit point is at an altitude" 
                     << " of " << ppath_step.pos(n-1,0)/1e3 << " km.";
                  throw runtime_error( os.str() );
                }
              if( lat_grid[imax_lat] < 90  && 
                    is_gridpos_at_index_i( ppath_step.gp_lat[n-1], imax_lat ) )
                {
                  ostringstream os;
                  os << "The path exits the atmosphere through the upper "
                     << "latitude end face.\nThe exit point is at an altitude" 
                     << " of " << ppath_step.pos(n-1,0)/1e3 << " km.";
                  throw runtime_error( os.str() );
                }

              // Longitude 
              // Note that it must be if and else if here. Otherwise e.g. -180 
              // will be shifted to 180 and then later back to -180.
              if( is_gridpos_at_index_i( ppath_step.gp_lon[n-1], 0 )  &&
                  ppath_step.los(n-1,1) < 0  &&  
                  abs( ppath_step.pos(n-1,1) ) < 90 )
                {
                  // Check if the longitude point can be shifted +360 degrees
                  if( lon_grid[imax_lon] - lon_grid[0] >= 360 )
                    {
                      ppath_step.pos(n-1,2) = ppath_step.pos(n-1,2) + 360;
                      gridpos( ppath_step.gp_lon[n-1], lon_grid, 
                                                       ppath_step.pos(n-1,2) );
                    }
                  else
                    {
                      ostringstream os;
                      os << "The path exits the atmosphere through the lower " 
                         << "longitude end face.\nThe exit point is at an "
                         << "altitude of " << ppath_step.pos(n-1,0)/1e3 
                         << " km.";
                      throw runtime_error( os.str() );
                    }
                }
              else if( 
                   is_gridpos_at_index_i( ppath_step.gp_lon[n-1], imax_lon ) &&
                   ppath_step.los(n-1,1) > 0  &&  
                   abs( ppath_step.pos(n-1,1) ) < 90 )
                {
                  // Check if the longitude point can be shifted -360 degrees
                  if( lon_grid[imax_lon] - lon_grid[0] >= 360 )
                    {
                      ppath_step.pos(n-1,2) = ppath_step.pos(n-1,2) - 360;
                      gridpos( ppath_step.gp_lon[n-1], lon_grid, 
                                                       ppath_step.pos(n-1,2) );
                    }
                  else
                    {
                      ostringstream os;
                      os << "The path exits the atmosphere through the upper "
                         << "longitude end face.\nThe exit point is at an "
                         << "altitude of " << ppath_step.pos(n-1,0)/1e3 
                         << " km.";
                      throw runtime_error( os.str() );
                    }
                }
            }
          
        
          // Check if there is an intersection with an active cloud box
          if( cloudbox_on )
            {
              double ipos = fractional_gp( ppath_step.gp_p[n-1] );

              if( ipos >= double( cloudbox_limits[0] )  && 
                  ipos <= double( cloudbox_limits[1] ) )
                {
                  if( atmosphere_dim == 1 )
                    { ppath_set_background( ppath_step, 3 ); }
                  else
                    {
                      ipos = fractional_gp( ppath_step.gp_lat[n-1] );

                      if( ipos >= double( cloudbox_limits[2] )  && 
                          ipos <= double( cloudbox_limits[3] ) )
                        {
                          if( atmosphere_dim == 2 )
                            { ppath_set_background( ppath_step, 3 ); }
                          else
                            {
                          
                              ipos = fractional_gp( ppath_step.gp_lon[n-1] );

                              if( ipos >= double( cloudbox_limits[4] )  && 
                                  ipos <= double( cloudbox_limits[5] ) )
                                { ppath_set_background( ppath_step, 3 ); } 
                            }
                        }
                    }
                }
            }
        }

      //--- Inside cloud box -------------------------------------------------
      else
        {
          // A first version just checked if point was at or outside any
          // boundary but numerical problems could cause that the start point
          // was taken as the exit point. So check of ppath direction had to be
          // added. Fractional distances used for this. 

          // Pressure dimension
          double ipos1 = fractional_gp( ppath_step.gp_p[n-1] );
          double ipos2 = fractional_gp( ppath_step.gp_p[n-2] );
          if( ipos1 <= double( cloudbox_limits[0] )  &&  ipos1 < ipos2 )
            { ppath_set_background( ppath_step, 3 ); }
              
          else if( ipos1 >= double( cloudbox_limits[1] )  &&  ipos1 > ipos2 )
            { ppath_set_background( ppath_step, 3 ); }

          else if( atmosphere_dim > 1 )
            {
              // Latitude dimension
              ipos1 = fractional_gp( ppath_step.gp_lat[n-1] );
              ipos2 = fractional_gp( ppath_step.gp_lat[n-2] );
              if( ipos1 <= double( cloudbox_limits[2] )  &&  ipos1 < ipos2 )  
                { ppath_set_background( ppath_step, 3 ); }

              else if( ipos1 >= double( cloudbox_limits[3] ) && ipos1 > ipos2 )
                { ppath_set_background( ppath_step, 3 ); }

              else if ( atmosphere_dim > 2 )
                {
                  // Longitude dimension
                  ipos1 = fractional_gp( ppath_step.gp_lon[n-1] );
                  ipos2 = fractional_gp( ppath_step.gp_lon[n-2] );
                  if( ipos1 <= double( cloudbox_limits[4] )  &&  ipos1 < ipos2 )
                    { ppath_set_background( ppath_step, 3 ); }

                  else if( ipos1 >= double( cloudbox_limits[5] )  &&
                           ipos1 > ipos2 )
                    { ppath_set_background( ppath_step, 3 ); }
                }
            }
        }
      //----------------------------------------------------------------------
      //---  End of boundary check
      //----------------------------------------------------------------------


      // Put new ppath_step in ppath_array
      ppath_array.push_back( ppath_step );

    } // End path steps


  // Combine all structures in ppath_array to form the return Ppath structure.
  //
  ppath_init_structure( ppath, atmosphere_dim, np );
  //
  const Index na = ppath_array.nelem();
  //
  if( na == 0 )    // No path, just the starting point
    {
      ppath_copy( ppath, ppath_step );
    }
 
  else   // Otherwise, merge the array elelments
    {
      np = 0;   // Now used as counter for points moved to ppath
      //
      for( Index i=0; i<na; i++ )
        {
          // For the first structure, the first point shall be included.
          // For later structures, the first point shall not be included, but
          // there will always be at least two points.
          
          Index n = ppath_array[i].np;

          // First index to include
          Index i1 = 1;
          if( i == 0 )
            { i1 = 0; }

          // Vectors and matrices that can be handled by ranges.
          ppath.r[ Range(np,n-i1) ] = ppath_array[i].r[ Range(i1,n-i1) ];
          ppath.pos( Range(np,n-i1), joker ) = 
                                   ppath_array[i].pos( Range(i1,n-i1), joker );
          ppath.los( Range(np,n-i1), joker ) = 
                                   ppath_array[i].los( Range(i1,n-i1), joker );
          ppath.nreal[ Range(np,n-i1) ] = 
                                        ppath_array[i].nreal[ Range(i1,n-i1) ];
          ppath.l_step[ Range(np-i1,n-1) ] = ppath_array[i].l_step; 

          // Grid positions must be handled by a loop
          for( Index j=i1; j<n; j++ )
            { ppath.gp_p[np+j-i1] = ppath_array[i].gp_p[j]; }
          if( atmosphere_dim >= 2 )
            {
              for( Index j=i1; j<n; j++ )
                { ppath.gp_lat[np+j-i1] = ppath_array[i].gp_lat[j]; }
            }
          if( atmosphere_dim == 3 )
            {
              for( Index j=i1; j<n; j++ )
                { ppath.gp_lon[np+j-i1] = ppath_array[i].gp_lon[j]; }
            }

          // Increase number of points done
          np += n - i1;
        }

      // Field just included once:
      // Constant and background can be taken from last path_step
      ppath.constant   = ppath_step.constant;
      ppath.background = ppath_step.background;
      // lspace extracted above
      ppath.lspace     = lspace;
    }
}
