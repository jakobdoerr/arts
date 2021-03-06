/* Copyright (C) 2002-2012
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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
  \file   m_rte.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2002-05-11 

  \brief  Workspace functions for solving clear sky radiative transfer.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "jacobian.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"

extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const String ABSSPECIES_MAINTAG;
extern const String TEMPERATURE_MAINTAG;
extern const String WIND_MAINTAG;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_P_GRID;
extern const Index GFIELD4_LAT_GRID;
extern const Index GFIELD4_LON_GRID;



/*===========================================================================
  === Workspace methods 
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void iyApplyUnit(
         Matrix&         iy,
   ArrayOfTensor4&       iy_aux,
   const Index&          stokes_dim,
   const Vector&         f_grid,
   const ArrayOfString&  iy_aux_vars,
   const String&         iy_unit,
   const Verbosity&)
{
  if( iy_unit == "1" )
    throw runtime_error( "No need to use this method with *iy_unit* = \"1\"." );

  if( max(iy(joker,0)) > 1e-3 )
    {
      ostringstream os;
      os << "The spectrum matrix *iy* is required to have original radiance\n"
         << "unit, but this seems not to be the case. This as a value above\n"
         << "1e-3 is found in *iy*.";
      throw runtime_error( os.str() );      
    }

  // Polarisation index variable
  ArrayOfIndex i_pol(stokes_dim);
  for( Index is=0; is<stokes_dim; is++ )
    { i_pol[is] = is + 1; }

  apply_iy_unit( iy, iy_unit, f_grid, 1, i_pol );
  
  for( Index i=0; i<iy_aux_vars.nelem(); i++ )
    {
      if( iy_aux_vars[i] == "iy"  ||  iy_aux_vars[i] == "Error"  || 
          iy_aux_vars[i] == "Error (uncorrelated)" )
        {
          if( iy_aux[i].nrows() > 1 )
            throw runtime_error( "Data marked as \"iy\" or \"Error\" "
                                 "have incorrect size." );
          for( Index j=0; j<iy_aux[i].ncols(); j++ )
            { apply_iy_unit( iy_aux[i](joker,joker,0,j), iy_unit, f_grid, 1, 
                                                                      i_pol ); }
        }
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyCalc(
         Workspace&        ws,
         Matrix&           iy,
         ArrayOfTensor4&   iy_aux,
         Ppath&            ppath,
   const Index&            atmfields_checked,
   const Index&            atmgeom_checked,
   const ArrayOfString&    iy_aux_vars,
   const Index&            iy_id,
   const Vector&           f_grid,
   const Tensor3&          t_field,
   const Tensor3&          z_field,
   const Tensor4&          vmr_field,
   const Index&            cloudbox_on,
   const Index&            cloudbox_checked,
   const Index&            scat_data_checked,
   const Vector&           rte_pos,
   const Vector&           rte_los,
   const Vector&           rte_pos2,
   const String&           iy_unit,  
   const Agenda&           iy_main_agenda,
   const Verbosity& )
{
  // Basics
  //
  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( atmgeom_checked != 1 )
    throw runtime_error( "The atmospheric geometry must be flagged to have "
                         "passed a consistency check (atmgeom_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );
  if( cloudbox_on )
    if( scat_data_checked != 1 )
      throw runtime_error( "The scattering data must be flagged to have "
                           "passed a consistency check (scat_data_checked=1)." );


  // iy_transmission is just input and can be left empty for first call
  Tensor3   iy_transmission(0,0,0);
  ArrayOfTensor3 diy_dx;
      
  iy_main_agendaExecute( ws, iy, iy_aux, ppath, diy_dx, 
                         1, iy_unit, iy_transmission, iy_aux_vars, 
                         iy_id, cloudbox_on, 0, t_field, 
                         z_field, vmr_field,
                         f_grid, rte_pos, rte_los, rte_pos2,
                         iy_main_agenda );
  
  // Don't allow NaNs (should suffice to check first stokes element)
  for( Index i=0; i<iy.nrows(); i++ )
    { 
      if( isnan(iy(i,0) ) )
        throw runtime_error( "One or several NaNs found in *iy*." );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandard2(
        Workspace&                  ws,
        Matrix&                     iy,
        ArrayOfMatrix&              iy_aux,
        ArrayOfTensor3&             diy_dx,
        Vector&                     ppvar_p,
        Vector&                     ppvar_t,
        Matrix&                     ppvar_t_nlte,
        Matrix&                     ppvar_vmr,
        Matrix&                     ppvar_wind,
        Matrix&                     ppvar_mag,         
        Matrix&                     ppvar_f,  
        Tensor3&                    ppvar_iy,  
  const Index&                      iy_id,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Tensor3&                    z_field,
  const Tensor3&                    t_field,
  const Tensor4&                    t_nlte_field,
  const Tensor4&                    vmr_field,
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Tensor3&                    wind_u_field,
  const Tensor3&                    wind_v_field,
  const Tensor3&                    wind_w_field,
  const Tensor3&                    mag_u_field,
  const Tensor3&                    mag_v_field,
  const Tensor3&                    mag_w_field,
  const Index&                      cloudbox_on,
  const String&                     iy_unit,
  const ArrayOfString&              iy_aux_vars,
  const Index&                      jacobian_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const Ppath&                      ppath,
  const Vector&                     rte_pos2, 
  const Agenda&                     propmat_clearsky_agenda,
  const Agenda&                     iy_main_agenda,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     iy_surface_agenda,
  const Agenda&                     iy_cloudbox_agenda,
  const Index&                      iy_agenda_call1,
  const Tensor3&                    iy_transmission,
  const Numeric&                    rte_alonglos_v,
  const Verbosity&                  verbosity )
{
  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;

  // Radiative background index
  const Index rbi = ppath_what_background( ppath );
  
  // Checks of input
  if( rbi < 1  ||  rbi > 9 )
    throw runtime_error( "ppath.background is invalid. Check your "
                         "calculation of *ppath*?" );
  if( !iy_agenda_call1  &&  np == 1  &&  rbi == 2  )
    throw runtime_error( "A secondary propagation path starting at the "
                         "surface and is going directly into the surface "
                         "is found. This is not allowed." );
  // iy_aux_vars checked below
  
  //  Init Jacobian quantities
  Index   j_analytical_do = 0;
  if( jacobian_do ) 
    { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  const Index     nq = j_analytical_do ? jacobian_quantities.nelem() : 0;
  ArrayOfTensor3  diy_dpath(nq); 
  ArrayOfIndex    jac_species_i(nq), jac_scat_i(nq), jac_is_t(nq), jac_wind_i(nq);
  ArrayOfIndex    jac_mag_i(nq), jac_other(nq);
  //
  // Flags for partial derivatives of propmat
  const PropmatPartialsData ppd(jacobian_quantities);
  //
  if( j_analytical_do )
    {
      const ArrayOfString  scat_species(0);
      const ArrayOfTensor4 dpnd_field_dx(nq);
      //
      rtmethods_jacobian_init( jac_species_i, jac_scat_i, jac_is_t, jac_wind_i,
                               jac_mag_i, jac_other, diy_dx,
                               diy_dpath,
                               ns, nf, np, nq, abs_species,
                               scat_species, dpnd_field_dx, ppd,
                               jacobian_quantities, iy_agenda_call1 );
    }
  
  // Init iy_aux and fill where possible
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize( naux );
  //
  for( Index i=0; i<naux; i++ )
    {
      iy_aux[i].resize(nf,ns); 
      
      if( iy_aux_vars[i] == "Transmission" )
        {} // Filled below
      else if( iy_aux_vars[i] == "Radiative background" )
        { iy_aux[i] = (Numeric)min( (Index)2, rbi-1 ); }
      else
        {
          ostringstream os;
          os << "The only allowed strings in *iy_aux_vars* are:\n"
             << "  \"Radiative background\"\n"
             << "  \"Optical depth\"\n"
             << "but you have selected: \"" << iy_aux_vars[i] << "\"";
          throw runtime_error( os.str() );      
        }
    }

  // Get atmospheric and radiative variables along the propagation path
  //
  Tensor3 J(np,nf,ns);
  Tensor4 trans_cumulat(np,nf,ns,ns), trans_partial(np,nf,ns,ns);
  Tensor5 dtrans_partial_dx_above(np,nq,nf,ns,ns);
  Tensor5 dtrans_partial_dx_below(np,nq,nf,ns,ns);
  Tensor4 dJ_dx(np,nq,nf,ns);
  //
  if( np == 1  &&  rbi == 1 )  // i.e. ppath is totally outside the atmosphere:
    {
      ppvar_p.resize(0);
      ppvar_t.resize(0);
      ppvar_vmr.resize(0,0);
      ppvar_t_nlte.resize(0,0);
      ppvar_wind.resize(0,0);
      ppvar_mag.resize(0,0);
      ppvar_f.resize(0,0);
      ppvar_iy.resize(0,0,0);
    }
  else
    {
      // ppvar_iy      
      ppvar_iy.resize(nf,ns,np);
      
      // Basic atmospheric variables
      get_ppath_atmvars( ppvar_p, ppvar_t, ppvar_t_nlte, ppvar_vmr,
                         ppvar_wind, ppvar_mag, 
                         ppath, atmosphere_dim, p_grid, t_field, t_nlte_field, 
                         vmr_field, wind_u_field, wind_v_field, wind_w_field,
                         mag_u_field, mag_v_field, mag_w_field );
      
      get_ppath_f( ppvar_f, ppath, f_grid,  atmosphere_dim, 
                   rte_alonglos_v, ppvar_wind );

      // Size radiative variables always used
      Vector B(nf);
      PropagationMatrix K_this(nf,ns), K_past(nf,ns), Kp(nf,ns);
      StokesVector a(nf,ns), S(nf,ns), Sp(nf,ns);
      ArrayOfIndex lte(np);

      // Init variables only used if analytical jacobians done
      Vector dB_dT(0);
      ArrayOfPropagationMatrix dK_this_dx(nq), dK_past_dx(nq), dKp_dx(nq);
      ArrayOfStokesVector da_dx(nq), dS_dx(nq), dSp_dx(nq);
      //
      if( j_analytical_do )
        {
          dB_dT.resize(nf);
          FOR_ANALYTICAL_JACOBIANS_DO
            (
              dK_this_dx[iq] = PropagationMatrix(nf,ns);
              dK_past_dx[iq] = PropagationMatrix(nf,ns);
              dKp_dx[iq]     = PropagationMatrix(nf,ns);
              da_dx[iq]      = StokesVector(nf,ns);
              dS_dx[iq]      = StokesVector(nf,ns);
              dSp_dx[iq]     = StokesVector(nf,ns);
            )
        }

      // Loop ppath points and determine radiative properties
      for( Index ip=0; ip<np; ip++ )
        {
          get_stepwise_blackbody_radiation( B,
                                            dB_dT,
                                            ppvar_f(joker,ip),
                                            ppvar_t[ip],
                                            ppd.do_temperature() );
          
          get_stepwise_clearsky_propmat( ws,
                                         K_this,
                                         S,
                                         lte[ip],
                                         dK_this_dx,
                                         dS_dx,
                                         propmat_clearsky_agenda,
                                         jacobian_quantities,
                                         ppd,
                                         ppvar_f(joker,ip),
                                         ppvar_mag(joker,ip),
                                         ppath.los(ip,joker),
                                         ppvar_t_nlte(joker,ip),
                                         ppvar_vmr(joker,ip),
                                         ppvar_t[ip],
                                         ppvar_p[ip],
                                         jac_species_i,
                                         j_analytical_do );

          if( j_analytical_do )
            {
              adapt_stepwise_partial_derivatives( dK_this_dx,
                                                  dS_dx,
                                                  jacobian_quantities,
                                                  ppvar_f(joker,ip),
                                                  ppath.los(ip,joker),
                                                  ppvar_vmr(joker,ip),
                                                  ppvar_t[ip],
                                                  ppvar_p[ip],
                                                  jac_species_i,
                                                  jac_wind_i,
                                                  lte[ip],
                                                  atmosphere_dim,
                                                  j_analytical_do );
            }

          // Here absorption equals extinction 
          a = K_this;
          if( j_analytical_do )
            {
              FOR_ANALYTICAL_JACOBIANS_DO          // we can only assign a PropMat to a
                (                                  // StokesVec, not full ArrayOf to
                  da_dx[iq] = dK_this_dx[iq];      // each other, hence loop over
                )                                  // jacobian quantities
            }

          get_stepwise_transmission_matrix(
                                 trans_cumulat(ip,joker,joker,joker),
                                 trans_partial(ip,joker,joker,joker),
                                 dtrans_partial_dx_above(ip,joker,joker,joker,joker),
                                 dtrans_partial_dx_below(ip,joker,joker,joker,joker),
                                 (ip>0)?
                                   trans_cumulat(ip-1,joker,joker,joker):
                                   Tensor3(0,0,0),
                                 K_past,
                                 K_this,
                                 dK_past_dx,
                                 dK_this_dx,
                                 (ip>0)?
                                   ppath.lstep[ip-1]:
                                   Numeric(1.0),
                                 ip==0 );

          get_stepwise_effective_source( J(ip,joker,joker),
                                         dJ_dx(ip,joker,joker,joker),
                                         K_this,
                                         a,
                                         S,
                                         dK_this_dx,
                                         da_dx,
                                         dS_dx,
                                         B,
                                         dB_dT,
                                         jacobian_quantities,
                                         j_analytical_do );
      
          swap( K_past, K_this );
          swap( dK_past_dx, dK_this_dx );
        }
    }

  // iy_transmission
  Tensor3 iy_trans_new;
  if( iy_agenda_call1 )
    { iy_trans_new = trans_cumulat(np-1,joker,joker,joker); }
  else
    { iy_transmission_mult( iy_trans_new, iy_transmission, 
                            trans_cumulat(np-1,joker,joker,joker) ); }

  // Copy transmission to iy_aux 
  for( Index i=0; i<naux; i++ )
    { if( iy_aux_vars[i] == "Transmission" )
        { for( Index iv=0; iv<nf; iv++ )
            { for( Index is=0; is<ns; is++ )
                { iy_aux[i](iv,is) = iy_trans_new(iv,is,is); }
    }   }   }

  // Radiative background
  get_iy_of_background( ws, iy, diy_dx, 
                        iy_trans_new, iy_id, jacobian_do, ppath, rte_pos2, 
                        atmosphere_dim, t_field, z_field, vmr_field, 
                        cloudbox_on, stokes_dim, f_grid, iy_unit,
                        iy_main_agenda, iy_space_agenda, iy_surface_agenda, 
                        iy_cloudbox_agenda, verbosity );
  //
  ppvar_iy(joker,joker,np-1) = iy;

  
  // Radiative transfer calculations
  if( np > 1 )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          Vector through_level(ns), from_level(ns);
          Vector dfrom_level_dx(ns);
          Matrix one_minus_transmission(ns,ns);
      
          for( Index ip=np-2; ip>=0; ip-- )
            {
        
              ConstMatrixView T = trans_partial(ip+1,iv,joker,joker);
        
              if( j_analytical_do )
                {
                  if( stokes_dim > 1 )
                    { id_mat(one_minus_transmission); }
                  else 
                    { one_minus_transmission = 1.; }
                  one_minus_transmission -= T;
                }

              from_level = J(ip,iv,joker);
              from_level += J(ip+1,iv,joker);
              from_level *= 0.5;
              through_level = iy(iv,joker);
              through_level -= from_level;

              if( j_analytical_do )
                {
                  FOR_ANALYTICAL_JACOBIANS_DO
                    (
                       get_diydx( diy_dpath[iq](ip,iv,joker), 
                                  diy_dpath[iq](ip+1,iv,joker), 
                                  one_minus_transmission,
                                  trans_cumulat(ip,iv,joker,joker), 
                                  dtrans_partial_dx_above(ip+1,iq,iv,joker,joker), 
                                  dtrans_partial_dx_below(ip+1,iq,iv,joker,joker), 
                                  through_level, 
                                  dJ_dx(ip,iq,iv,joker), 
                                  dJ_dx(ip+1,iq,iv,joker),
                                  stokes_dim );
                     )
                }
              
              // Equation is I1 = T (I0 - 0.5(J_1+J_2)) + 0.5(J_1+J_2)
              mult( iy(iv,joker), T, through_level );
              iy(iv,joker) += from_level;

              ppvar_iy(iv,joker,ip) = iy(iv,joker);
            }
        }
    }


  // Finalize analytical Jacobians
  if( j_analytical_do )
    {
      rtmethods_jacobian_finalisation( diy_dx, diy_dpath,
                                       ns, nf, np, atmosphere_dim, ppath,
                                       ppvar_p, ppvar_t, ppvar_vmr,
                                       iy_agenda_call1, iy_transmission,
                                       jacobian_quantities, jac_species_i,
                                       jac_is_t );
    }

  // Radiance unit conversions
  if( iy_agenda_call1 )
    {
      rtmethods_unit_conversion( iy, diy_dx, ppvar_iy,
                                 ns, np, f_grid, ppath, jacobian_quantities,
                                 j_analytical_do, iy_unit );
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandard(
         Workspace&                  ws,
         Matrix&                     iy,
         ArrayOfTensor4&             iy_aux,
         Ppath&                      ppath,
         ArrayOfTensor3&             diy_dx,
   const Index&                      iy_id,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    t_nlte_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Tensor3&                    wind_u_field,
   const Tensor3&                    wind_v_field,
   const Tensor3&                    wind_w_field,
   const Tensor3&                    mag_u_field,
   const Tensor3&                    mag_v_field,
   const Tensor3&                    mag_w_field,
   const Index&                      cloudbox_on,
   const String&                     iy_unit,
   const ArrayOfString&              iy_aux_vars,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const Agenda&                     ppath_agenda,
   const Agenda&                     propmat_clearsky_agenda,
   const Agenda&                     iy_main_agenda,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     iy_surface_agenda,
   const Agenda&                     iy_cloudbox_agenda,
   const Index&                      iy_agenda_call1,
   const Tensor3&                    iy_transmission,
   const Vector&                     rte_pos,      
   const Vector&                     rte_los,      
   const Vector&                     rte_pos2, 
   const Numeric&                    rte_alonglos_v,      
   const Numeric&                    ppath_lmax,
   const Numeric&                    ppath_lraytrace,
   const Verbosity&                  verbosity )
{
  // Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, ppath_lmax, ppath_lraytrace,
                       rte_pos, rte_los, rte_pos2, 
                       cloudbox_on, 0, t_field, z_field, vmr_field, f_grid, 
                       ppath_agenda );
  //
  if( !iy_agenda_call1  &&  ppath.np == 1  &&  ppath_what_background( ppath ) == 2  )
    throw runtime_error( "A secondary propagation path starting at the "
                         "surface and is going directly into the surface "
                         "is found. This is not allowed." );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index nq = jacobian_quantities.nelem();

  
  // Brief definition of some of the internal variables (and to a start also questions)
  //
  // abs_per_species
  //   The absorption for individual species. Used to fill iy_aux.
  // auxAbsSpecies,  auxAbsIsp, auxVmrSpecies, auxVmrIsp
  //   Coding of which iy_aux elements that cover particle absorption resp.
  //   VMR, and the corresponding index among abs_species. Left empty if iy_aux
  //   does not contain any quantity of concern.
  // auxAbsSum, absPressure ...
  //   Set to -1 if the variable of concern not is requested to be part of
  //   iy_aux. Otherwise set to the position inside iy_aux.
  // diy_dpath
  //   The derivate of iy with respect to changes at the ppath points.
  // dppath_ext_dx
  //   The derivatives of the extinction matrix at the ppath points.
  // dppath_nlte_source_dx
  //   The derivatives of the NLTE source term at the ppath points.
  // dtrans_partial_dx_above/below
  //   The derivatives of the transmission matrix at the ppath points.
  //   Below means the derivatives are from the path point itself and above
  //   means the derivatives are from the neighboring ppath point.  Together,
  //   Above and Below makes the layer through which the transfer is carried
  //   out.
  // iaps
  //   The index of species for which abs_per_species shall be filled. 
  // jac_mag_i: Works as jac_species_i, but uses JacobianType::MagField... and JacobianType::None
  //    for coding of each element.
  // jac_other
  //    Works as jac_species_i, but uses JacobianType::Other... and JacobianType::None
  //    for "other" propmat derivatives
  // jac_is_t
  //    Works as jac_species_i, but uses JacobianType::Temperature... and JacobianType::None
  //    for coding of each element.
  // jac_species_i
  //   An array of length nq. Elements are the species index where a retrieval 
  //   quantity is a species otherwise set to -1.  
  // jac_wind_i 
  //    Works as jac_species_i, but uses JacobianType::WindField... and JacobianType::None
  //    for coding of each element.
  // jac_to_integrate 
  //    Track keeping of integration variable
  // ppath_p, ppath_vmr, ppath_ext ...
  //    Holds the pressure, vmr, extinction matrix etc. at each point of ppath.
  // ppd
  //   Flags and indexing for the jacobian_quantities that works in propmat clearsky
  // scalar_tau
  //   The total optical thickness of the ppath (a scalar value).
  // trans_cumulat, trans_partial
  //   The transmission (as Mueller matrices) for each ppath step resp. from
  //   the end to each point. See further *get_ppath_trans*.

  
  //###### jacobian part #######################################################
  // Initialise analytical jacobians (diy_dx and help variables)
  //
  Index           j_analytical_do = 0;
  ArrayOfTensor3  diy_dpath; 
  ArrayOfIndex    jac_species_i(0), jac_scat_i(0), jac_is_t(0), jac_wind_i(0);
  ArrayOfIndex    jac_mag_i(0), jac_other(0); 
  // Flags for partial derivatives of propmat
  const PropmatPartialsData ppd(jacobian_quantities);
  //
  if( jacobian_do ) 
    { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( !j_analytical_do )
    { diy_dx.resize( 0 ); }
  else 
    {
      diy_dpath.resize( nq ); 
      jac_species_i.resize( nq ); 
      jac_scat_i.resize( nq ); 
      jac_is_t.resize( nq ); 
      jac_wind_i.resize( nq );  
      jac_mag_i.resize( nq ); 
      jac_other.resize(nq);
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dpath[iq].resize( np, nf, ns ); 
            diy_dpath[iq] = 0.0;
      )

      const ArrayOfString scat_species(0);
      
      get_pointers_for_analytical_jacobians( jac_species_i, jac_scat_i, jac_is_t, 
                                             jac_wind_i, jac_mag_i, 
                                             jacobian_quantities,
                                             abs_species, scat_species );

      FOR_ANALYTICAL_JACOBIANS_DO(
        jac_other[iq] = ppd.is_this_propmattype(iq)?Index(JacobianType::Other):Index(JacobianType::None);
      )
      
      if( iy_agenda_call1 )
        {
          diy_dx.resize( nq ); 
          //
          bool any_affine;
          ArrayOfArrayOfIndex jacobian_indices;
          jac_ranges_indices( jacobian_indices, any_affine,
                              jacobian_quantities, true );
          //
          FOR_ANALYTICAL_JACOBIANS_DO( diy_dx[iq].resize( 
            jacobian_indices[iq][1]-jacobian_indices[iq][0]+1, nf, ns ); 
            diy_dx[iq] = 0.0;
          )
        }
    } 
  //###########################################################################
  
  
  //=== iy_aux part ===========================================================
  Index auxPressure    = -1,
        auxTemperature = -1,
        auxAbsSum      = -1,
        auxBackground  = -1,
        auxIy          = -1,
        auxTrans       = -1,
        auxOptDepth    = -1;
  ArrayOfIndex iaps(0);
  ArrayOfIndex auxAbsSpecies(0), auxAbsIsp(0);
  ArrayOfIndex auxVmrSpecies(0), auxVmrIsp(0);
  //
  {
    const Index naux = iy_aux_vars.nelem();
    iy_aux.resize( naux );
    //
    for( Index i=0; i<naux; i++ )
      {
        if( iy_aux_vars[i] == "Pressure" )
          { auxPressure = i;      iy_aux[i].resize( 1, 1, 1, np ); }
        else if( iy_aux_vars[i] == "Temperature" )
          { auxTemperature = i;   iy_aux[i].resize( 1, 1, 1, np ); }
        else if( iy_aux_vars[i].substr(0,13) == "VMR, species " )
          { 
            Index ispecies;
            istringstream is(iy_aux_vars[i].substr(13,2));
            is >> ispecies;
            if( ispecies < 0  ||  ispecies>=abs_species.nelem() )
              {
                ostringstream os;
                os << "You have selected VMR of species with index "
                   << ispecies << ".\nThis species does not exist!";
                throw runtime_error( os.str() );
              }
            auxVmrSpecies.push_back(i);
            auxVmrIsp.push_back(ispecies);
            iy_aux[i].resize( 1, 1, 1, np );               
          }
        else if( iy_aux_vars[i] == "Absorption, summed" )
          { auxAbsSum = i;   iy_aux[i].resize( nf, ns, ns, np ); }
        else if( iy_aux_vars[i].substr(0,20) == "Absorption, species " )
          { 
            Index ispecies;
            istringstream is(iy_aux_vars[i].substr(20,2));
            is >> ispecies;
            if( ispecies < 0  ||  ispecies>=abs_species.nelem() )
              {
                ostringstream os;
                os << "You have selected absorption species with index "
                   << ispecies << ".\nThis species does not exist!";
                throw runtime_error( os.str() );
              }
            auxAbsSpecies.push_back(i);
            const Index ihit = find_first( iaps, ispecies );
            if( ihit >= 0 )
              { auxAbsIsp.push_back( ihit ); }
            else
              { 
                iaps.push_back(ispecies); 
                auxAbsIsp.push_back( iaps.nelem()-1 ); 
              }
            iy_aux[i].resize( nf, ns, ns, np );               
          }
        else if( iy_aux_vars[i] == "Radiative background" )
          { auxBackground = i;   iy_aux[i].resize( nf, 1, 1, 1 ); }
        else if( iy_aux_vars[i] == "iy"   &&  auxIy < 0 )
          { auxIy = i;           iy_aux[i].resize( nf, ns, 1, np ); }
        else if( iy_aux_vars[i] == "Transmission"   &&  auxTrans < 0 )
          { auxTrans = i;        iy_aux[i].resize( nf, ns, ns, np ); }
        else if( iy_aux_vars[i] == "Optical depth" )
          { auxOptDepth = i;     iy_aux[i].resize( nf, 1, 1, 1 ); }
        else if( iy_aux_vars[i].substr(0,14) == "Mass content, " )
          { iy_aux[i].resize( 0, 0, 0, 0 ); }
        else if( iy_aux_vars[i].substr(0,10) == "PND, type " )
          { iy_aux[i].resize( 0, 0, 0, 0 ); }
        else
          {
            ostringstream os;
            os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
               << "\"\nThis choice is not recognised.";
            throw runtime_error( os.str() );
          }
      }
  }
  //===========================================================================


  // Get atmospheric and attenuation quantities for each ppath point/step
  //
  // "atmvars"
  Vector    ppath_p, ppath_t;
  Matrix    ppath_vmr, ppath_wind, ppath_mag, ppath_f, ppath_t_nlte, ppath_pnd_dummy;
  // Attenuation vars
  ArrayOfPropagationMatrix ppath_ext;
  ArrayOfPropagationMatrix pnd_ext_mat_dummy;
  Tensor3                  pnd_abs_vec;
  ArrayOfArrayOfPropagationMatrix dppath_ext_dx;
  ArrayOfArrayOfPropagationMatrix abs_per_species;
  Tensor5   dtrans_partial_dx_above, dtrans_partial_dx_below;
  Tensor4   trans_partial, trans_cumulat;
  ArrayOfArrayOfStokesVector dppath_nlte_source_dx;
  ArrayOfStokesVector   ppath_nlte_source;
  Matrix    ppath_blackrad, dppath_blackrad_dt;
  Vector    scalar_tau;
  ArrayOfIndex         lte, clear2cloudy_dummy;
  ArrayOfArrayOfIndex  extmat_case;
  ArrayOfMatrix        dummy_ppath_dpnd_dx;
  ArrayOfTensor4       dummy_dpnd_field_dx;
  Array<ArrayOfArrayOfSingleScatteringData> scat_data_single_dummy;
  const ArrayOfArrayOfSingleScatteringData scat_data_dummy;
  const Tensor4 pnd_field_dummy;
  const ArrayOfIndex cloudbox_limits_dummy;
  //
  if( np > 1 )
    {
      get_ppath_atmvars( ppath_p, ppath_t, ppath_t_nlte, ppath_vmr,
                         ppath_wind, ppath_mag, 
                         ppath, atmosphere_dim, p_grid, t_field, t_nlte_field, 
                         vmr_field, wind_u_field, wind_v_field, wind_w_field,
                         mag_u_field, mag_v_field, mag_w_field );      

      get_ppath_f( ppath_f, ppath, f_grid,  atmosphere_dim, 
                   rte_alonglos_v, ppath_wind );

      get_ppath_pmat_and_tmat( ws, ppath_ext, ppath_nlte_source, lte, abs_per_species,
                               dppath_ext_dx, dppath_nlte_source_dx, trans_partial,
                               dtrans_partial_dx_above, dtrans_partial_dx_below,
                               extmat_case, clear2cloudy_dummy, trans_cumulat,
                               scalar_tau, pnd_ext_mat_dummy, pnd_abs_vec, ppath_pnd_dummy,
                               dummy_ppath_dpnd_dx, scat_data_single_dummy,
                               propmat_clearsky_agenda, jacobian_quantities,
                               ppd, ppath, ppath_p,  ppath_t, ppath_t_nlte, ppath_vmr, 
                               ppath_mag, ppath_f, f_grid, 
                               jac_species_i, jac_is_t, jac_wind_i, jac_mag_i,
                               jac_other, iaps,
                               scat_data_dummy,
                               pnd_field_dummy, dummy_dpnd_field_dx,
                               cloudbox_limits_dummy,
                               atmosphere_dim, stokes_dim,
                               jacobian_do, false, verbosity);
      
      get_ppath_blackrad( ppath_blackrad, ppath, ppath_t, ppath_f );

      get_dppath_blackrad_dt( dppath_blackrad_dt, ppath_t, ppath_f, jac_is_t, 
                              j_analytical_do );
    }
  else // For cases inside the cloudbox, or totally outside the atmosphere,
    {  // set zero optical thickness and unit transmission
      scalar_tau.resize( nf );
      scalar_tau = 0;
      trans_cumulat.resize( nf, ns, ns, np );
      for( Index iv=0; iv<nf; iv++ )
        { id_mat( trans_cumulat(iv,joker,joker,np-1) ); }
    }

  // iy_transmission
  //
  Tensor3 iy_trans_new;
  //
  if( iy_agenda_call1 )
    { iy_trans_new = trans_cumulat(joker,joker,joker,np-1); }
  else
    { iy_transmission_mult( iy_trans_new, iy_transmission, 
                            trans_cumulat(joker,joker,joker,np-1) ); }
  
  // Radiative background
  //
  get_iy_of_background( ws, iy, diy_dx, 
                        iy_trans_new, iy_id, jacobian_do, ppath, rte_pos2, 
                        atmosphere_dim, t_field, z_field, vmr_field, 
                        cloudbox_on, stokes_dim, f_grid, iy_unit,
                        iy_main_agenda, iy_space_agenda, iy_surface_agenda, 
                        iy_cloudbox_agenda, verbosity );


  //=== iy_aux part ===========================================================
  // Fill parts of iy_aux that are defined even for np=1.
  // Radiative background
  if( auxBackground >= 0 ) 
    { iy_aux[auxBackground](joker,0,0,0) = (Numeric)min( (Index)2,
                                              ppath_what_background(ppath)-1); }
  // Radiance 
  if( auxIy >= 0 ) 
    { iy_aux[auxIy](joker,joker,0,np-1) = iy; }
  // Transmission variables
  if( auxTrans >= 0 ) 
    { 
      if( np == 1 )
        { for( Index iv=0; iv<nf; iv++ ) {
            id_mat( iy_aux[auxTrans](iv,joker,joker,0) ); } }
      else
        { iy_aux[auxTrans] = trans_cumulat; }
    }
  if( auxOptDepth >= 0 ) 
    { iy_aux[auxOptDepth](joker,0,0,0) = scalar_tau; }
  //===========================================================================

  
  // Do RT calculations
  //
  if( np > 1 )
  {
    // Temperature disturbance, K   
    //
    // (This variable is used in some parts of the T-jacobian)
    //
    const Numeric   dt = 0.1;     


    //=== iy_aux part =======================================================
    // iy_aux for point np-1:
    
    // Pressure
    if( auxPressure >= 0 )  { iy_aux[auxPressure](0,0,0,np-1) = ppath_p[np-1]; }
    
    // Temperature
    if( auxTemperature >= 0 )  { iy_aux[auxTemperature](0,0,0,np-1) = ppath_t[np-1]; }
      
    // VMR
    for( Index j=0; j<auxVmrSpecies.nelem(); j++ ) { iy_aux[auxVmrSpecies[j]](0,0,0,np-1) = ppath_vmr(auxVmrIsp[j],np-1); }
    
    // Absorption
    if( auxAbsSum >= 0 ) 
    { 
      for( Index is1=0; is1<ns; is1++ )
      {
        iy_aux[auxAbsSum](joker,is1,is1,np-1) = ppath_ext[np-1].Kjj();
        if(is1 == 1)
        {
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) = iy_aux[auxAbsSum](joker,is1-1,is1,np-1) = ppath_ext[np-1].K12();
        }
        else if(is1 == 2)
        {
          iy_aux[auxAbsSum](joker,is1,is1-2,np-1) = iy_aux[auxAbsSum](joker,is1-2,is1,np-1) = ppath_ext[np-1].K13();
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) = iy_aux[auxAbsSum](joker,is1-1,is1,np-1) = ppath_ext[np-1].K23();
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) *= -1;
        }
        else if(is1 == 3)
        {
          iy_aux[auxAbsSum](joker,is1,is1-3,np-1) = iy_aux[auxAbsSum](joker,is1-3,is1,np-1) = ppath_ext[np-1].K14();
          iy_aux[auxAbsSum](joker,is1,is1-2,np-1) = iy_aux[auxAbsSum](joker,is1-2,is1,np-1) = ppath_ext[np-1].K24();
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) = iy_aux[auxAbsSum](joker,is1-1,is1,np-1) = ppath_ext[np-1].K34();
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) *= -1;
          iy_aux[auxAbsSum](joker,is1,is1-2,np-1) *= -1;
        }
      }
    }
    
    // Absorption by species
    for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
    {
      for( Index is1=0; is1<ns; is1++ )
      {
        iy_aux[auxAbsSum](joker,is1,is1,np-1) = abs_per_species[np-1][auxAbsIsp[j]].Kjj();
        if(is1 == 1)
        {
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) = iy_aux[auxAbsSum](joker,is1-1,is1,np-1) = abs_per_species[np-1][auxAbsIsp[j]].K12();
        }
        else if(is1 == 2)
        {
          iy_aux[auxAbsSum](joker,is1,is1-2,np-1) = iy_aux[auxAbsSum](joker,is1-2,is1,np-1) = abs_per_species[np-1][auxAbsIsp[j]].K13();
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) = iy_aux[auxAbsSum](joker,is1-1,is1,np-1) = abs_per_species[np-1][auxAbsIsp[j]].K23();
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) *= -1;
        }
        else if(is1 == 3)
        {
          iy_aux[auxAbsSum](joker,is1,is1-3,np-1) = iy_aux[auxAbsSum](joker,is1-3,is1,np-1) = abs_per_species[np-1][auxAbsIsp[j]].K14();
          iy_aux[auxAbsSum](joker,is1,is1-2,np-1) = iy_aux[auxAbsSum](joker,is1-2,is1,np-1) = abs_per_species[np-1][auxAbsIsp[j]].K24();
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) = iy_aux[auxAbsSum](joker,is1-1,is1,np-1) = abs_per_species[np-1][auxAbsIsp[j]].K34();
          iy_aux[auxAbsSum](joker,is1,is1-1,np-1) *= -1;
          iy_aux[auxAbsSum](joker,is1,is1-2,np-1) *= -1;
        }
      }
    }
    
    // Radiance for this point is handled above
    //=======================================================================


    //=======================================================================
    // Loop ppath steps
    for( Index ip=np-2; ip>=0; ip-- )
    {
      // Path step average of B: Bbar
      //
      Vector bbar(nf);
      //
      for( Index iv=0; iv<nf; iv++ )  
        { bbar[iv] = 0.5 * ( ppath_blackrad(iv,ip) +
                              ppath_blackrad(iv,ip+1) ); }

      // Extra variables for non-LTE
      //
      const bool nonlte = lte[ip] == 0 or lte[ip+1] == 0; 
      //
      StokesVector sourcebar;
      PropagationMatrix extbar;
      //
      if( nonlte )
      { 
        sourcebar = StokesVector(ppath_nlte_source[ip], ppath_nlte_source[ip+1], 0.5);
        extbar    = PropagationMatrix(ppath_ext[ip], ppath_ext[ip+1], 0.5);
      }

      
      //### jacobian part #################################################
      if( j_analytical_do )
      { 
        // The variable names used below refer largely to the
        // nomenclature used in the "Clear-sky Jacobians" chapter in AUG.
        // For example, si is the Stokes vector for position i.
        
        // Difference between local stokes and Planck (si-bi)
        Matrix sibi(nf,ns);
        
        // nlte terms
        Matrix nlte_inv(ns,ns);
        Vector nlte_sibi(ns,0.0);
        
        for( Index iv=0; iv<nf; iv++ )
        {
          if(nonlte) // Then sibi is difference between local Stokes and local Source
          {
            extbar.MatrixInverseAtPosition(nlte_inv, iv);
            mult(nlte_sibi, nlte_inv, sourcebar.VectorAtPosition(iv));
          }
          
          sibi(iv,0) = iy(iv,0) - bbar[iv] - nlte_sibi[0];
          for( Index is=1; is<ns; is++ )
          { 
              sibi(iv,is) = iy(iv,is) - nlte_sibi[is];
          }
        }

          
        for( Index iq=0; iq<nq; iq++ ) 
        {
          if( jacobian_quantities[iq].Analytical() )
          {
            if( jac_species_i[iq] >= 0 || jac_wind_i[iq] ||
              jac_mag_i[iq] || jac_other[iq] || jac_is_t[iq])
            {
              /* 
                *   We need to do blackbody derivation for temperature.
                *   Should we do this for wind (frequency) as well, then 
                *   add or statement here for jac_wind_i[iq] and turn
                *   dppath_blackrad_dt into a Tensor3 of size [nq,nf,np].
                */
              
              const bool this_is_t = jac_is_t[iq], 
              this_is_hse = this_is_t ? jacobian_quantities[iq].Subtag() == "HSE on" : false;
              
              get_diydx_replacement(diy_dpath[iq](ip, joker, joker),
                                    diy_dpath[iq](ip+1, joker, joker),
                                    iy,
                                    sibi,
                                    ppath_nlte_source[ip],
                                    ppath_nlte_source[ip+1],
                                    dppath_nlte_source_dx[ip][iq],
                                    dppath_nlte_source_dx[ip+1][iq],
                                    ppath_ext[ip],
                                    ppath_ext[ip+1],
                                    dppath_ext_dx[ip][iq],
                                    dppath_ext_dx[ip+1][iq],
                                    trans_partial(joker,joker,joker,ip),
                                    dtrans_partial_dx_below(iq,joker,joker,joker,ip),
                                    dtrans_partial_dx_above(iq,joker,joker,joker,ip),
                                    trans_cumulat(joker,joker,joker,ip  ),
                                    trans_cumulat(joker,joker,joker,ip+1),
                                    ppath_t[ip  ],
                                    ppath_t[ip+1],
                                    dt,
                                    dppath_blackrad_dt(joker,ip  ),
                                    dppath_blackrad_dt(joker,ip+1),
                                    ppath.lstep[ip],
                                    this_is_t,
                                    this_is_hse,
                                    nonlte);
            } // if this iq is analytical
          } // if this analytical
        } // for iq
      } // if any analytical
      //###################################################################


      // Spectrum at end of ppath step
      emission_rtstep_replacement(iy, stokes_dim, bbar, extmat_case[ip],
                                  trans_partial(joker,joker,joker,ip),
                                  nonlte, extbar, sourcebar );

      //=== iy_aux part ===================================================
      // Pressure
      if( auxPressure >= 0 )  { iy_aux[auxPressure](0,0,0,ip) = ppath_p[ip]; }
      
      // Temperature
      if( auxTemperature >= 0 )  { iy_aux[auxTemperature](0,0,0,ip) = ppath_t[ip]; }
      
      // VMR
      for( Index j=0; j<auxVmrSpecies.nelem(); j++ ) { iy_aux[auxVmrSpecies[j]](0,0,0,ip) = ppath_vmr(auxVmrIsp[j],ip); }
      
      // Absorption
      if( auxAbsSum >= 0 ) 
      { 
        for( Index is1=0; is1<ns; is1++ )
        {
          iy_aux[auxAbsSum](joker,is1,is1,ip) = ppath_ext[ip].Kjj();
          if(is1 == 1)
          {
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) = iy_aux[auxAbsSum](joker,is1-1,is1,ip) = ppath_ext[ip].K12();
          }
          else if(is1 == 2)
          {
            iy_aux[auxAbsSum](joker,is1,is1-2,ip) = iy_aux[auxAbsSum](joker,is1-2,is1,ip) = ppath_ext[ip].K13();
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) = iy_aux[auxAbsSum](joker,is1-1,is1,ip) = ppath_ext[ip].K23();
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) *= -1;
          }
          else if(is1 == 3)
          {
            iy_aux[auxAbsSum](joker,is1,is1-3,ip) = iy_aux[auxAbsSum](joker,is1-3,is1,ip) = ppath_ext[ip].K14();
            iy_aux[auxAbsSum](joker,is1,is1-2,ip) = iy_aux[auxAbsSum](joker,is1-2,is1,ip) = ppath_ext[ip].K24();
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) = iy_aux[auxAbsSum](joker,is1-1,is1,ip) = ppath_ext[ip].K34();
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) *= -1;
            iy_aux[auxAbsSum](joker,is1,is1-2,ip) *= -1;
          }
        }
      }
      
      // Absorption by species
      for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
      {
        for( Index is1=0; is1<ns; is1++ )
        {
          iy_aux[auxAbsSum](joker,is1,is1,ip) = abs_per_species[ip][auxAbsIsp[j]].Kjj();
          if(is1 == 1)
          {
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) = iy_aux[auxAbsSum](joker,is1-1,is1,ip) = abs_per_species[ip][auxAbsIsp[j]].K12();
          }
          else if(is1 == 2)
          {
            iy_aux[auxAbsSum](joker,is1,is1-2,ip) = iy_aux[auxAbsSum](joker,is1-2,is1,ip) = abs_per_species[ip][auxAbsIsp[j]].K13();
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) = iy_aux[auxAbsSum](joker,is1-1,is1,ip) = abs_per_species[ip][auxAbsIsp[j]].K23();
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) *= -1;
          }
          else if(is1 == 3)
          {
            iy_aux[auxAbsSum](joker,is1,is1-3,ip) = iy_aux[auxAbsSum](joker,is1-3,is1,ip) = abs_per_species[ip][auxAbsIsp[j]].K14();
            iy_aux[auxAbsSum](joker,is1,is1-2,ip) = iy_aux[auxAbsSum](joker,is1-2,is1,ip) = abs_per_species[ip][auxAbsIsp[j]].K24();
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) = iy_aux[auxAbsSum](joker,is1-1,is1,ip) = abs_per_species[ip][auxAbsIsp[j]].K34();
            iy_aux[auxAbsSum](joker,is1,is1-1,ip) *= -1;
            iy_aux[auxAbsSum](joker,is1,is1-2,ip) *= -1;
          }
        }
      }
      // Radiance 
      if( auxIy >= 0 ) 
        { iy_aux[auxIy](joker,joker,0,ip) = iy; }
    } // path point loop
    //=======================================================================


    //### jacobian part #####################################################
    // Map jacobians from ppath to retrieval grids
    // (this operation corresponds to the term Dx_i/Dx)
    if( j_analytical_do )
    { 
      // Weight with iy_transmission
      if( !iy_agenda_call1 )
      {
        Matrix X, Y; 
        //
        FOR_ANALYTICAL_JACOBIANS_DO
        ( 
          Y.resize(ns,diy_dpath[iq].npages());
          for( Index iv=0; iv<nf; iv++ )
          { 
                X = transpose( diy_dpath[iq](joker,iv,joker) );
                mult( Y, iy_transmission(iv,joker,joker), X );
                diy_dpath[iq](joker,iv,joker) = transpose( Y );
          }
        )
      }

      // Map to retrieval grids
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                            diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
      )
    }
    //#######################################################################
  } // if np>1


  // Unit conversions
  if( iy_agenda_call1 )
    {
      // Determine refractive index to use for the n2 radiance law
      Numeric n = 1.0; // First guess is that sensor is in space
      //
      if( ppath.end_lstep == 0 ) // If true, sensor inside the atmosphere
        { n = ppath.nreal[np-1]; }

      // Polarisation index variable
      ArrayOfIndex i_pol(ns);
      for( Index is=0; is<ns; is++ )
        { i_pol[is] = is + 1; }

      // Jacobian part (must be converted to Tb before iy for PlanckBT)
      // 
      if( j_analytical_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO( apply_iy_unit2( diy_dx[iq], iy, iy_unit,
                                                       f_grid, n, i_pol ); )
        } 

      // iy
      apply_iy_unit( iy, iy_unit, f_grid, n, i_pol );

      // iy_aux
      for( Index q=0; q<iy_aux.nelem(); q++ )
      {
        if( iy_aux_vars[q] == "iy")
        { 
          for( Index ip=0; ip<ppath.np; ip++ )
          { apply_iy_unit( iy_aux[q](joker,joker,0,ip), iy_unit, f_grid, 
            ppath.nreal[ip], i_pol ); }
        }
      }
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyLoopFrequencies(
         Workspace&        ws,
         Matrix&           iy,
         ArrayOfTensor4&   iy_aux,
         Ppath&            ppath,
         ArrayOfTensor3&   diy_dx,
   const ArrayOfString&    iy_aux_vars,
   const Index&            stokes_dim,
   const Vector&           f_grid,
   const Index&            atmosphere_dim,
   const Vector&           p_grid,         
   const Vector&           lat_grid,
   const Vector&           lon_grid,
   const Vector&           lat_true,
   const Vector&           lon_true,
   const Tensor3&          t_field,
   const Tensor3&          z_field,
   const Tensor4&          vmr_field,
   const Matrix&           z_surface,
   const Numeric&          ppath_lmax,
   const Numeric&          ppath_lraytrace,
   const Index&            cloudbox_on,
   const ArrayOfIndex&     cloudbox_limits, 
   const Tensor4&          pnd_field,
   const Index&            iy_agenda_call1,
   const String&           iy_unit,  
   const Tensor3&          iy_transmission,
   const Vector&           rte_pos,
   const Vector&           rte_los,
   const Vector&           rte_pos2,
   const Index&            jacobian_do,
   const Agenda&           iy_sub_agenda,
   const Verbosity& )
{
  // Throw error if unsupported features are requested
  if( !iy_agenda_call1 )
    throw runtime_error( "Recursive usage not possible (iy_agenda_call1 must be 1)." );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty." );

  const Index nf = f_grid.nelem();

  for( Index i=0; i<nf; i++ )
    {
      // Variables for 1 frequency
      Matrix         iy1;
      ArrayOfTensor4 iy_aux1; 
      ArrayOfTensor3 diy_dx1;
      
      iy_sub_agendaExecute( ws, iy1, iy_aux1, ppath, diy_dx1, iy_agenda_call1,
                            iy_unit, iy_transmission, iy_aux_vars, 0,
                            Vector(1,f_grid[i]), atmosphere_dim, p_grid,
                            lat_grid, lon_grid, lat_true, lon_true,
                            t_field, z_field, vmr_field, z_surface,
                            ppath_lmax, ppath_lraytrace,
                            cloudbox_on, cloudbox_limits, pnd_field,
                            jacobian_do, rte_pos, rte_los, rte_pos2, iy_sub_agenda );

      // After first frequency, give output its size
      if( i == 0 )
        {
          iy.resize( nf, stokes_dim );
          //
          iy_aux.resize( iy_aux1.nelem() );
          for( Index q=0; q<iy_aux1.nelem(); q++ )
            {
              if( iy_aux1[q].ncols() > 1 )
                throw runtime_error( "When using this method, *iy_aux_vars* "
                        "is not allowed to include along-the-path variables." );
              iy_aux[q].resize(nf,iy_aux1[q].npages(),iy_aux1[q].nrows(),1);
            }
          //
          diy_dx.resize( diy_dx1.nelem() );
          for( Index q=0; q<diy_dx1.nelem(); q++ )
            { diy_dx[q].resize( diy_dx1[q].npages(), nf, stokes_dim ); }
        }

      // Copy to output variables
      iy(i,joker) = iy1(0,joker);
      for( Index q=0; q<iy_aux1.nelem(); q++ )
        { iy_aux[q](i,joker,joker,0) = iy_aux1[q](0,joker,joker,0); }
      for( Index q=0; q<diy_dx1.nelem(); q++ )
        { diy_dx[q](joker,i,joker) = diy_dx1[q](joker,0,joker); }
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyMC(
         Workspace&                  ws,
         Matrix&                     iy,
         ArrayOfTensor4&             iy_aux,
         ArrayOfTensor3&             diy_dx,
   const Index&                      iy_agenda_call1,
   const Tensor3&                    iy_transmission,
   const Vector&                     rte_pos,      
   const Vector&                     rte_los,      
   const ArrayOfString&              iy_aux_vars,
   const Index&                      jacobian_do,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const Vector&                     refellipsoid,
   const Matrix&                     z_surface,
   const Index&                      cloudbox_on,
   const ArrayOfIndex&               cloudbox_limits,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const ArrayOfArrayOfSingleScatteringData&   scat_data,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     surface_rtprop_agenda,
   const Agenda&                     propmat_clearsky_agenda, 
   const Agenda&                     ppath_step_agenda, 
   const Numeric&                    ppath_lmax, 
   const Numeric&                    ppath_lraytrace, 
   const Tensor4&                    pnd_field,
   const String&                     iy_unit,
   const Numeric&                    mc_std_err,
   const Index&                      mc_max_time,
   const Index&                      mc_max_iter,
   const Index&                      mc_min_iter,
   const Numeric&                    mc_taustep_limit,
   const Verbosity&                  verbosity)
{
  // Throw error if unsupported features are requested
  if( atmosphere_dim != 3 )
    throw runtime_error( 
                "Only 3D atmospheres are allowed (atmosphere_dim must be 3)" );
  if( !cloudbox_on )
    throw runtime_error( 
                    "The cloudbox must be activated (cloudbox_on must be 1)" );
  if( jacobian_do )
    throw runtime_error( 
        "This method does not provide any jacobians (jacobian_do must be 0)" );
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty" );


  // Size output variables
  //
  const Index   nf = f_grid.nelem();
  //
  iy.resize( nf, stokes_dim );
  diy_dx.resize(0);
  //
  //=== iy_aux part ===========================================================
  Index auxError = -1;
  {
    const Index naux = iy_aux_vars.nelem();
    iy_aux.resize( naux );
    //
    for( Index i=0; i<naux; i++ )
      {
        if( iy_aux_vars[i] == "Error (uncorrelated)" )
          { auxError = i;      iy_aux[i].resize( nf, stokes_dim, 1, 1 ); }
        else
          {
            ostringstream os;
            os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
               << "\"\nThis choice is not recognised.";
            throw runtime_error( os.str() );
          }
      }
  }
  //===========================================================================

  // Some MC variables are only local here
  //
  MCAntenna mc_antenna;
  mc_antenna.set_pencil_beam();

  // Pos and los must be matrices 
  Matrix pos(1,3), los(1,2);
  //
  pos(0,joker) = rte_pos;
  los(0,joker) = rte_los;

  Workspace l_ws (ws);
  Agenda l_ppath_step_agenda (ppath_step_agenda);
  Agenda l_iy_space_agenda (iy_space_agenda);
  Agenda l_propmat_clearsky_agenda (propmat_clearsky_agenda);
  Agenda l_surface_rtprop_agenda (surface_rtprop_agenda);

  String fail_msg;
  bool failed = false;

  if (nf)
#pragma omp parallel for                   \
  if (!arts_omp_in_parallel() && nf > 1)   \
  firstprivate(l_ws, l_ppath_step_agenda, l_iy_space_agenda, \
               l_propmat_clearsky_agenda, l_surface_rtprop_agenda)
  for( Index f_index=0; f_index<nf; f_index++ )
    {
      if (failed) continue;

      try {
        // Seed reset for each loop. If not done, the errors
        // appear to be highly correlated.
        Index    mc_seed;
        MCSetSeedFromTime( mc_seed, verbosity );

        Vector       y, mc_error;
        Index        mc_iteration_count;
        Tensor3      mc_points;
        ArrayOfIndex mc_scat_order, mc_source_domain;

        MCGeneral( l_ws, y, mc_iteration_count, mc_error, mc_points, 
                   mc_scat_order, mc_source_domain, mc_antenna,
                   f_grid, f_index, pos, los, stokes_dim, atmosphere_dim,
                   l_ppath_step_agenda, ppath_lmax, ppath_lraytrace, l_iy_space_agenda, 
                   l_surface_rtprop_agenda, l_propmat_clearsky_agenda, 
                   p_grid, lat_grid, lon_grid, z_field, 
                   refellipsoid, z_surface, t_field, vmr_field,
                   cloudbox_on, cloudbox_limits,
                   pnd_field, scat_data, 1, 1, 1, 1, iy_unit,
                   mc_seed, mc_std_err, mc_max_time, mc_max_iter,
                   mc_min_iter, mc_taustep_limit, 1, verbosity);
          //cout << "Error: "      << mc_error << endl;
          //cout << "N photons: " << mc_iteration_count << endl;
        assert( y.nelem() == stokes_dim );

        iy(f_index,joker) = y;
          
        if( auxError >= 0 ) 
          { iy_aux[auxError](f_index,joker,0,0) = mc_error; }
      } catch (runtime_error e) {
        ostringstream os;
        os << "Error for f_index = " << f_index << " (" << f_grid[f_index] 
           << ")" << endl << e.what();
#pragma omp critical (iyMC_fail)
          { failed = true; fail_msg = os.str(); }
          continue;
      }
    }

  if (failed)
    throw runtime_error(fail_msg);
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyReplaceFromAux(
         Matrix&            iy,
   const ArrayOfTensor4&    iy_aux,
   const ArrayOfString&     iy_aux_vars,
   const Index&             jacobian_do,
   const String&            aux_var,
   const Verbosity& )
{
  if( iy_aux.nelem() != iy_aux_vars.nelem() )
    throw runtime_error( "*iy_aux* and *iy_aux_vars* must have the same "
                         "number of elements." );
  
  if( jacobian_do )
    throw runtime_error( "This method can not provide any jacobians and "
                         "*jacobian_do* must be 0." );

  bool ready = false;

  for( Index i=0; i<iy_aux.nelem() && !ready; i++ )
    {
      if( iy_aux_vars[i] == aux_var )
        {
          if( iy_aux[i].nrows() > 1  ||  iy_aux[i].ncols() > 1 )
            {
              throw runtime_error( "If an auxiliary variable shall be inserted "
                     "in *iy*, its row and page dimensions must have size 1." );
            }
          if( iy_aux[i].nbooks() != iy.nrows() )
            {
              throw runtime_error( "If an auxiliary variable shall be inserted "
                                   "in *iy*, its frequency dimension must match"
                                   "the length of existing *iy*." );
            }

          iy = 0;
          
          for( Index iv=0; iv<iy.nrows(); iv++ ) 
            {
              for( Index is=0; is<iy_aux[i].npages(); is++ )
                {
                  iy(iv,is) = iy_aux[i](iv,is,0,0);
                }
            }
          
          ready = true;
        }
    }

  if( !ready )
    throw runtime_error( "The selected auxiliary variable to insert in *iy* "
                         "is either not defined at all or is not set." );

}





/* Workspace method: Doxygen documentation will be auto-generated */
void iy_auxFillParticleVariables(
         ArrayOfTensor4&       iy_aux,
   const Index&                atmfields_checked,
   const Index&                cloudbox_checked,
   const Index&                atmosphere_dim,
   const Index&                cloudbox_on,
   const ArrayOfIndex&         cloudbox_limits,
   const Tensor4&              pnd_field,
   const Matrix&               particle_masses,
   const Ppath&                ppath,
   const ArrayOfString&        iy_aux_vars,
   const Verbosity& )
{
  // Some sizes
  const Index np = ppath.np; 
  const Index naux = iy_aux_vars.nelem();

  // Input checks
  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );
  if( !cloudbox_on )
    throw runtime_error( 
                    "The cloudbox must be activated (cloudbox_on must be 1)" );
  if( iy_aux.nelem() != naux )
    throw runtime_error( "*iy_aux_vars* and *iy_aux* must have the same array "
                         "length. (You can not call this WSM before the main "
                         "iy-WSM.)" );

  // Analayse iy_aux_vars
  ArrayOfIndex auxPartCont(0), auxPartContI(0);
  ArrayOfIndex auxPartField(0), auxPartFieldI(0);
  //
  for( Index i=0; i<naux; i++ )
    {
      if( iy_aux_vars[i].substr(0,14) == "Mass content, " )
        { 
          Index icont;
          istringstream is(iy_aux_vars[i].substr(14,2));
          is >> icont;
          if( icont < 0  ||  icont>=particle_masses.ncols() )
            {
              ostringstream os;
              os << "You have selected particle mass content category with "
                 << "index " << icont << ".\nThis category is not defined!";
              throw runtime_error( os.str() );
            }
          auxPartCont.push_back(i);
          auxPartContI.push_back(icont);
          iy_aux[i].resize( 1, 1, 1, np );
        }
      else if( iy_aux_vars[i].substr(0,10) == "PND, type " )
        { 
          Index ip;
          istringstream is(iy_aux_vars[i].substr(10,2));
          is >> ip;
          if( ip < 0  ||  ip>=pnd_field.nbooks() )
            {
              ostringstream os;
              os << "You have selected particle number density field with "
                 << "index " << ip << ".\nThis field is not defined!";
              throw runtime_error( os.str() );
            }
          auxPartField.push_back(i);
          auxPartFieldI.push_back(ip);
          iy_aux[i].resize( 1, 1, 1, np );
        }
    }

  if( auxPartCont.nelem() + auxPartField.nelem() > 0 )
    {
      // PND along the ppath
      Matrix ppath_pnd( pnd_field.nbooks(), np, 0 );
      //
      for( Index ip=0; ip<np; ip++ )
        {
          Matrix itw( 1, Index(pow(2.0,Numeric(atmosphere_dim))) );

          ArrayOfGridPos gpc_p(1), gpc_lat(1), gpc_lon(1);
          GridPos gp_lat, gp_lon;
          if( atmosphere_dim >= 2 ) { gridpos_copy( gp_lat, ppath.gp_lat[ip] );}
          if( atmosphere_dim == 3 ) { gridpos_copy( gp_lon, ppath.gp_lon[ip] );}
          if( is_gp_inside_cloudbox( ppath.gp_p[ip], gp_lat, gp_lon, 
                                     cloudbox_limits, true, atmosphere_dim ) )
            { 
              interp_cloudfield_gp2itw( itw(0,joker), 
                                        gpc_p[0], gpc_lat[0], gpc_lon[0], 
                                        ppath.gp_p[ip], gp_lat, gp_lon,
                                        atmosphere_dim, cloudbox_limits );
              for( Index i=0; i<pnd_field.nbooks(); i++ )
                {
                  interp_atmfield_by_itw( ppath_pnd(i,ip), atmosphere_dim,
                                          pnd_field(i,joker,joker,joker), 
                                          gpc_p, gpc_lat, gpc_lon, itw );
                }
            }
        }
      
      // Loop ppath steps
      for( Index ip=0; ip<np; ip++ )
        {
          // Particle mass content
          for( Index j=0; j<auxPartCont.nelem(); j++ )
            { iy_aux[auxPartCont[j]](0,0,0,ip) = ppath_pnd(joker,ip) *
                                      particle_masses(joker,auxPartContI[j]); }
          // Particle number density
          for( Index j=0; j<auxPartField.nelem(); j++ )
            { iy_aux[auxPartField[j]](0,0,0,ip) = 
                                              ppath_pnd(auxPartFieldI[j],ip); }
        }  
    }
}


void yCalc_mblock_loop_body(
         bool&                       failed,
         String&                     fail_msg,
         ArrayOfArrayOfVector&       iyb_aux_array,
         Workspace&                  ws,
         Vector&                     y,
         Vector&                     y_f,
         ArrayOfIndex&               y_pol,
         Matrix&                     y_pos,
         Matrix&                     y_los,
         Matrix&                     y_geo,
         Matrix&                     jacobian,
   const Index&                      atmosphere_dim,
   const Tensor3&                    t_field,
   const Tensor3&                    z_field,
   const Tensor4&                    vmr_field,
   const Index&                      cloudbox_on,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Matrix&                     transmitter_pos,
   const Matrix&                     mblock_dlos_grid,
   const Sparse&                     sensor_response,
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Matrix&                     sensor_response_dlos,
   const String&                     iy_unit,   
   const Agenda&                     iy_main_agenda,
   const Agenda&                     geo_pos_agenda,
   const Agenda&                     jacobian_agenda,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const ArrayOfString&              iy_aux_vars,
   const Verbosity&                  verbosity,
   const Index&                      mblock_index,
   const Index&                      n1y,
   const Index&                      j_analytical_do)
{
    try
    {
        // Calculate monochromatic pencil beam data for 1 measurement block
        //
        Vector          iyb, iyb_error, yb(n1y);
        ArrayOfMatrix   diyb_dx;
        Matrix          geo_pos_matrix;
        //
        iyb_calc(ws, iyb, iyb_aux_array[mblock_index], diyb_dx, geo_pos_matrix,
                 mblock_index, atmosphere_dim, t_field, z_field, vmr_field,
                 cloudbox_on, stokes_dim, f_grid, sensor_pos, sensor_los,
                 transmitter_pos, mblock_dlos_grid, 
                 iy_unit, iy_main_agenda, geo_pos_agenda,
                 j_analytical_do, jacobian_quantities,
                 jacobian_indices, iy_aux_vars, verbosity);

        // Apply sensor response matrix on iyb, and put into y
        //
        const Range rowind = get_rowindex_for_mblock(sensor_response,
                                                     mblock_index);
        const Index row0   = rowind.get_start();
        //
        mult( yb, sensor_response, iyb );
        //
        y[rowind] = yb;  // *yb* also used below, as input to jacobian_agenda

        // Fill information variables. And search for NaNs in *y*.
        //
        for( Index i=0; i<n1y; i++ )
          {
            const Index ii = row0 + i; 
            if( isnan(y[ii] ) )
              throw runtime_error( "One or several NaNs found in *y*." );
            y_f[ii]          = sensor_response_f[i];
            y_pol[ii]        = sensor_response_pol[i];
            y_pos(ii,joker)  = sensor_pos(mblock_index,joker);
            y_los(ii,joker)  = sensor_los(mblock_index,joker);
            y_los(ii,0)     += sensor_response_dlos(i,0);
            if( sensor_response_dlos.ncols() > 1 )
              { y_los(ii,1) += sensor_response_dlos(i,1); }
          }
        
        // Apply sensor response matrix on diyb_dx, and put into jacobian
        // (that is, analytical jacobian part)
        //
        if( j_analytical_do )
          {
            FOR_ANALYTICAL_JACOBIANS_DO(
              mult(jacobian(rowind,
                            Range(jacobian_indices[iq][0],
                                jacobian_indices[iq][1]-jacobian_indices[iq][0]+1)),
                                sensor_response, diyb_dx[iq] );
            )
          }

        // Calculate remaining parts of *jacobian*
        //
        if( jacobian_do )
          {
            jacobian_agendaExecute( ws, jacobian, mblock_index, iyb, yb,
                                    jacobian_agenda );
          }


        // Handle geo-positioning
        if( !isnan(geo_pos_matrix(0,0)) )  // No data are flagged as NaN
          {
            // Find bore sigtht direction be probing sensor_response
            const Index   nf   = f_grid.nelem();
            const Index   nlos = mblock_dlos_grid.nrows();
            const Index   niyb = nf * nlos * stokes_dim;
            ArrayOfIndex i_of_max( n1y );
            Vector max_contr( n1y, NAN );
            for( Index ilos=0; ilos<nlos; ilos++ )
              {
                Vector itry( niyb, 0 );
                itry[Range(ilos*nf*stokes_dim,nf,stokes_dim)] = 1;
                Vector ytry( n1y );
                mult( ytry, sensor_response, itry );
                for( Index i=0; i<n1y; i++ )
                  {
                    if( ytry[i] > max_contr[i] )
                      {
                        max_contr[i] = ytry[i];
                        i_of_max[i]  = ilos;
                      }
                  }
              }

            // Extract geo_pos_matrix for found bore-sights
            for( Index i=0; i<n1y; i++ )
              { y_geo(row0+i,joker) = geo_pos_matrix(i_of_max[i],joker); }
          }
    }

    catch (runtime_error e)
    {
#pragma omp critical (yCalc_fail)
        { fail_msg = e.what(); failed = true; }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void yCalc(
         Workspace&                  ws,
         Vector&                     y,
         Vector&                     y_f,
         ArrayOfIndex&               y_pol,
         Matrix&                     y_pos,
         Matrix&                     y_los,
         ArrayOfVector&              y_aux,
         Matrix&                     y_geo,
         Matrix&                     jacobian,
   const Index&                      atmgeom_checked,
   const Index&                      atmfields_checked,
   const Index&                      atmosphere_dim,
   const Tensor3&                    t_field,
   const Tensor3&                    z_field,
   const Tensor4&                    vmr_field,
   const Index&                      cloudbox_on,
   const Index&                      cloudbox_checked,
   const Index&                      scat_data_checked,
   const Index&                      sensor_checked,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Matrix&                     transmitter_pos,
   const Matrix&                     mblock_dlos_grid,
   const Sparse&                     sensor_response,
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Matrix&                     sensor_response_dlos,
   const String&                     iy_unit,   
   const Agenda&                     iy_main_agenda,
   const Agenda&                     geo_pos_agenda,
   const Agenda&                     jacobian_agenda,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfString&              iy_aux_vars,
   const Verbosity&                  verbosity )
{
  CREATE_OUT3;

  // Basics
  //
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  //
  if ( f_grid.empty() )
    { throw runtime_error ( "The frequency grid is empty." ); }
  chk_if_increasing ( "f_grid", f_grid );
  if( f_grid[0] <= 0) 
    { throw runtime_error( "All frequencies in *f_grid* must be > 0." ); }
  //
  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( atmgeom_checked != 1 )
    throw runtime_error( "The atmospheric geometry must be flagged to have "
                         "passed a consistency check (atmgeom_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );
  if( cloudbox_on )
    if( scat_data_checked != 1 )
      throw runtime_error( "The scattering data must be flagged to have "
                           "passed a consistency check (scat_data_checked=1)." );
  if( sensor_checked != 1 )
    throw runtime_error( "The sensor variables must be flagged to have "
                         "passed a consistency check (sensor_checked=1)." );


  // Some sizes
  const Index   nf      = f_grid.nelem();
  const Index   nlos    = mblock_dlos_grid.nrows();
  const Index   n1y     = sensor_response.nrows();
  const Index   nmblock = sensor_pos.nrows();
  const Index   niyb    = nf * nlos * stokes_dim;


  //---------------------------------------------------------------------------
  // Allocations and resizing
  //---------------------------------------------------------------------------

  // Resize *y* and *y_XXX*
  //
  y.resize( nmblock*n1y );
  y_f.resize( nmblock*n1y );
  y_pol.resize( nmblock*n1y );
  y_pos.resize( nmblock*n1y, sensor_pos.ncols() );
  y_los.resize( nmblock*n1y, sensor_los.ncols() );
  y_geo.resize( nmblock*n1y, atmosphere_dim );
  y_geo = NAN;   // Will be replaced if relavant data are provided (*geo_pos*)

  // For y_aux we don't know the number of quantities, and we need to 
  // store all output
  ArrayOfArrayOfVector  iyb_aux_array( nmblock );

  // Jacobian variables
  //
  Index               j_analytical_do = 0;
  ArrayOfArrayOfIndex jacobian_indices;
  //
  if( jacobian_do )
    {
      bool any_affine;
      jac_ranges_indices( jacobian_indices, any_affine,
                          jacobian_quantities, true );
      jacobian.resize( nmblock*n1y, 
                       jacobian_indices[jacobian_indices.nelem()-1][1]+1 );
      jacobian = 0;
      //
      FOR_ANALYTICAL_JACOBIANS_DO(
        j_analytical_do  = 1; 
      )
    }
  else
    { jacobian.resize( 0, 0 ); }


  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  String fail_msg;
  bool failed = false;

  if (nmblock >= arts_omp_get_max_threads()
      || (nf <= nmblock && nmblock >= nlos))
  {
      out3 << "  Parallelizing mblock loop (" << nmblock << " iterations)\n";

      // We have to make a local copy of the Workspace and the agendas because
      // only non-reference types can be declared firstprivate in OpenMP
      Workspace l_ws (ws);
      Agenda l_jacobian_agenda (jacobian_agenda);
      Agenda l_iy_main_agenda (iy_main_agenda);
      Agenda l_geo_pos_agenda (geo_pos_agenda);

#pragma omp parallel for                         \
  firstprivate(l_ws, l_jacobian_agenda, l_iy_main_agenda, l_geo_pos_agenda)
      for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
      {
          // Skip remaining iterations if an error occurred
          if (failed) continue;

          yCalc_mblock_loop_body( failed, fail_msg, iyb_aux_array, l_ws,
                                  y, y_f, y_pol, y_pos, y_los, y_geo, jacobian,
                                  atmosphere_dim, t_field, z_field,
                                  vmr_field,
                                  cloudbox_on, stokes_dim, f_grid, 
                                  sensor_pos, sensor_los, transmitter_pos,
                                  mblock_dlos_grid, sensor_response,
                                  sensor_response_f, sensor_response_pol,
                                  sensor_response_dlos, iy_unit,
                                  l_iy_main_agenda, l_geo_pos_agenda,
                                  l_jacobian_agenda,
                                  jacobian_do, jacobian_quantities,
                                  jacobian_indices, iy_aux_vars, verbosity,
                                  mblock_index, n1y, j_analytical_do );          
      }  // End mblock loop
  }
  else
  {
     out3 << "  Not parallelizing mblock loop (" << nmblock << " iterations)\n";

     for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
      {
          // Skip remaining iterations if an error occurred
          if (failed) continue;

          yCalc_mblock_loop_body( failed, fail_msg, iyb_aux_array, ws,
                                  y, y_f, y_pol, y_pos, y_los, y_geo, jacobian,
                                  atmosphere_dim, t_field, z_field,
                                  vmr_field,
                                  cloudbox_on, stokes_dim, f_grid, 
                                  sensor_pos, sensor_los, transmitter_pos,
                                  mblock_dlos_grid, sensor_response,
                                  sensor_response_f, sensor_response_pol,
                                  sensor_response_dlos, iy_unit,
                                  iy_main_agenda, geo_pos_agenda,
                                  jacobian_agenda,
                                  jacobian_do, jacobian_quantities,
                                  jacobian_indices, iy_aux_vars, verbosity,
                                  mblock_index, n1y, j_analytical_do ); 
      }  // End mblock loop
  }

  // Rethrow exception if a runtime error occurred in the mblock loop
  if (failed) throw runtime_error(fail_msg);

  // Compile y_aux
  //
  const Index nq = iyb_aux_array[0].nelem();
  y_aux.resize( nq );
  //
  for( Index q=0; q<nq; q++ )
    {
      y_aux[q].resize( nmblock*n1y );
      //
      for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
        {
          const Range rowind = get_rowindex_for_mblock( sensor_response, 
                                                                mblock_index );
          const Index row0   = rowind.get_start();

          // The sensor response must be applied in a special way for
          // uncorrelated errors. Schematically: sqrt( H.^2 * y.^2 )
          if( iy_aux_vars[q] == "Error (uncorrelated)" )
            {
              for( Index i=0; i<n1y; i++ )
                {
                  const Index row = row0+i;
                  y_aux[q][row] = 0;
                  for( Index j=0; j<niyb; j++ )
                    { y_aux[q][row] += pow( sensor_response(i,j) * 
                            iyb_aux_array[mblock_index][q][j], (Numeric)2.0 ); }
                  y_aux[q][row] = sqrt( y_aux[q][row] );              
                }
            }
          else
            { mult( y_aux[q][rowind], sensor_response,
                                      iyb_aux_array[mblock_index][q] ); }
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void yCalcAppend(
         Workspace&                  ws,
         Vector&                     y,
         Vector&                     y_f,
         ArrayOfIndex&               y_pol,
         Matrix&                     y_pos,
         Matrix&                     y_los,
         ArrayOfVector&              y_aux,
         Matrix&                     y_geo,
         Matrix&                     jacobian,
         ArrayOfRetrievalQuantity&   jacobian_quantities,
   const Index&                      atmfields_checked,
   const Index&                      atmgeom_checked,
   const Index&                      atmosphere_dim,
   const Tensor3&                    t_field,
   const Tensor3&                    z_field,
   const Tensor4&                    vmr_field,
   const Index&                      cloudbox_on,
   const Index&                      cloudbox_checked,
   const Index&                      scat_data_checked,
   const Index&                      sensor_checked,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Matrix&                     transmitter_pos,
   const Matrix&                     mblock_dlos_grid,
   const Sparse&                     sensor_response,
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Matrix&                     sensor_response_dlos,
   const String&                     iy_unit,   
   const Agenda&                     iy_main_agenda,
   const Agenda&                     geo_pos_agenda,
   const Agenda&                     jacobian_agenda,
   const Index&                      jacobian_do,
   const ArrayOfString&              iy_aux_vars,
   const ArrayOfRetrievalQuantity&   jacobian_quantities_copy,
   const Index&                      append_instrument_wfs,
   const Verbosity&                  verbosity )
{
  // The jacobian indices of old and new part (without transformations)
  ArrayOfArrayOfIndex jacobian_indices, jacobian_indices_copy;
  {
    bool any_affine; 
    jac_ranges_indices( jacobian_indices_copy, any_affine,
                        jacobian_quantities_copy, true );
    jac_ranges_indices( jacobian_indices, any_affine,
                        jacobian_quantities, true );
  }

  // Check consistency of data representing first measurement
  const Index n1   = y.nelem();
        Index nrq1 = 0; 
  if( y.empty() )
    throw runtime_error( "Input *y* is empty. Use *yCalc*" );
  if( y_f.nelem() != n1 )
    throw runtime_error( "Lengths of input *y* and *y_f* are inconsistent." );
  if( y_pol.nelem() != n1 )
    throw runtime_error( "Lengths of input *y* and *y_pol* are inconsistent." );
  if( y_pos.nrows() != n1 )
    throw runtime_error( "Sizes of input *y* and *y_pos* are inconsistent." );
  if( y_los.nrows() != n1 )
    throw runtime_error( "Sizes of input *y* and *y_los* are inconsistent." );
  if( y_geo.nrows() != n1 )
    throw runtime_error( "Sizes of input *y* and *y_geo* are inconsistent." );
  if( jacobian_do )
    {
      nrq1 = jacobian_quantities_copy.nelem();
      if( jacobian.nrows() != n1 )
        throw runtime_error( "Sizes of *y* and *jacobian* are inconsistent." );
      if( jacobian.ncols() != jacobian_indices_copy[nrq1-1][1]+1 )
        throw runtime_error( "Size of input *jacobian* and size implied " 
                             "*jacobian_quantities_copy* are inconsistent." );
    }

  // Calculate new measurement
  //
  Vector        y2, y_f2;
  Matrix        y_pos2, y_los2, y_geo2, jacobian2;
  ArrayOfIndex  y_pol2;
  ArrayOfVector y_aux2;
  //
  yCalc( ws, y2, y_f2, y_pol2, y_pos2, y_los2, y_aux2, y_geo2, jacobian2,
         atmfields_checked, atmgeom_checked, atmosphere_dim, t_field,
         z_field, vmr_field, cloudbox_on,
         cloudbox_checked, scat_data_checked, sensor_checked,
         stokes_dim, f_grid, sensor_pos, sensor_los, transmitter_pos,
         mblock_dlos_grid, sensor_response,
         sensor_response_f, sensor_response_pol, sensor_response_dlos, 
         iy_unit, iy_main_agenda, geo_pos_agenda, 
         jacobian_agenda, jacobian_do, jacobian_quantities,
         iy_aux_vars, verbosity );

  // Consistency checks
  if( y_pos.ncols() != y_pos2.ncols() )
    throw runtime_error( 
          "Different number of columns in *y_pos* between the measurements." );
  if( y_los.ncols() != y_los2.ncols() )
    throw runtime_error( 
          "Different number of columns in *y_los* between the measurements." );
  

  // y and y_XXX
  //
  const Index         n2 = y2.nelem(); 
  //
  {
    // Make copy of old measurement
    const Vector        y1=y, y_f1=y_f;
    const Matrix        y_pos1=y_pos, y_los1=y_los, y_geo1=y_geo;
    const ArrayOfIndex  y_pol1=y_pol;
    const ArrayOfVector y_aux1=y_aux;
    //
    y.resize( n1+n2 );
    y[Range(0,n1)] = y1;   y[Range(n1,n2)] = y2; 
    //
    y_f.resize( n1+n2 );
    y_f[Range(0,n1)] = y_f1;   y_f[Range(n1,n2)] = y_f2; 
    //
    y_pos.resize( n1+n2, y_pos1.ncols() );
    y_pos(Range(0,n1),joker) = y_pos1;   y_pos(Range(n1,n2),joker) = y_pos2; 
    //
    y_los.resize( n1+n2, y_los1.ncols() );
    y_los(Range(0,n1),joker) = y_los1;   y_los(Range(n1,n2),joker) = y_los2; 
    //
    y_geo.resize( n1+n2, y_geo1.ncols() );
    y_geo(Range(0,n1),joker) = y_geo1;   y_geo(Range(n1,n2),joker) = y_geo2; 
    //
    y_pol.resize( n1+n2 );
    for( Index i=0; i<n1; i++ )
      { y_pol[i] = y_pol1[i]; }
    for( Index i=0; i<n2; i++ )
      { y_pol[n1+i] = y_pol2[i]; }

    // y_aux
    const Index na1 = y_aux1.nelem();
    const Index na2 = y_aux2.nelem();
    const Index na  = max(na1,na2);
    //
    y_aux.resize( na );
    //
    for( Index a=0; a<na; a++ )
      {
        y_aux[a].resize( n1+n2 );        
        if( a < na1 )
          { y_aux[a][Range(0,n1)] = y_aux1[a]; }
        else
          { y_aux[a][Range(0,n1)] = 0; }
        if( a < na2 )
          { y_aux[a][Range(n1,n2)] = y_aux2[a]; }
        else
          { y_aux[a][Range(n1,n2)] = 0; }
      }
  }

  // Jacobian and friends
  if( jacobian_do )
    {
      // Put in old jacobian_quantities and jacobian_indices as first part in
      // new version of these variables
      ArrayOfRetrievalQuantity  jacobian_quantities2 = jacobian_quantities;
      ArrayOfArrayOfIndex       jacobian_indices2    = jacobian_indices;
      //
      jacobian_quantities = jacobian_quantities_copy;
      jacobian_indices    = jacobian_indices_copy;

      // Loop new jacobian_quantities to determine how new jacobian data shall
      // be inserted
      //
      const Index    nrq2 = jacobian_quantities2.nelem();
      ArrayOfIndex   map_table(nrq2);
      //
      for( Index q2=0; q2<nrq2; q2++ )
        {

          Index pos = -1;

          // Compare to old quantities, to determine if append shall be 
          // considered. Some special checks performed here, grids checked later
          if( jacobian_quantities2[q2].MainTag() == ABSSPECIES_MAINTAG   ||
              jacobian_quantities2[q2].MainTag() == TEMPERATURE_MAINTAG  ||
              jacobian_quantities2[q2].MainTag() == WIND_MAINTAG         ||
              append_instrument_wfs )
            {
              for( Index q1=0; q1<nrq1; q1++ && pos < 0 )
                {
                  if( jacobian_quantities2[q2].MainTag() ==
                      jacobian_quantities_copy[q1].MainTag() )
                    {
                      // Absorption species
                      if( jacobian_quantities2[q2].MainTag() == 
                                                           ABSSPECIES_MAINTAG )
                        {
                          if( jacobian_quantities2[q2].Subtag() ==
                              jacobian_quantities_copy[q1].Subtag() )
                            {
                              if( jacobian_quantities2[q2].Mode() ==
                                  jacobian_quantities_copy[q1].Mode() )
                                { pos = q1; }
                              else
                                {
                                  ostringstream os;
                                  os << "Jacobians for " 
                                     << jacobian_quantities2[q2].MainTag()<<"/" 
                                     << jacobian_quantities2[q2].Subtag() 
                                     << " shall be appended.\nThis requires "
                                     << "that the same retrieval unit is used "
                                     << "but it seems that this requirement is "
                                     << "not met.";
                                  throw runtime_error(os.str());
                                }
                            }
                        }
                      // Temperature
                      else if( jacobian_quantities2[q2].MainTag() == 
                                                          TEMPERATURE_MAINTAG )
                        {
                          if( jacobian_quantities2[q2].Subtag() ==
                              jacobian_quantities_copy[q1].Subtag() )
                            { pos = q1; }
                          else
                            {
                              ostringstream os;
                              os << "Jacobians for " 
                                 << jacobian_quantities2[q2].MainTag()<<"/" 
                                 << jacobian_quantities2[q2].Subtag() 
                                 << " shall be appended.\nThis requires "
                                 << "that HSE is either ON or OFF for both "
                                 << "parts but it seems that this requirement "
                                 << "is not met.";
                              throw runtime_error(os.str());
                            }
                        }
                      // Other
                      else if( jacobian_quantities2[q2].Subtag() ==
                               jacobian_quantities_copy[q1].Subtag() )
                        { pos = q1; }
                    }
                }
            }


          // New quantity
          if( pos < 0 )
            {
              map_table[q2] = jacobian_quantities.nelem();
              jacobian_quantities.push_back( jacobian_quantities2[q2] );
              ArrayOfIndex indices(2);
              indices[0] = jacobian_indices[jacobian_indices.nelem()-1][1]+1;
              indices[1] = indices[0] + jacobian_indices2[q2][1] -
                                        jacobian_indices2[q2][0];
              jacobian_indices.push_back( indices );
            }
          // Existing quantity
          else
            {
              map_table[q2] = pos;
              // Check if grids are equal
              ArrayOfVector grids1 = jacobian_quantities_copy[pos].Grids();
              ArrayOfVector grids2 = jacobian_quantities2[q2].Grids();
              bool any_wrong = false;
              if( grids1.nelem() != grids2.nelem() )
                { any_wrong = true; }
              else
                {
                  for( Index g=0; g<grids1.nelem(); g++ )
                    {
                      if( grids1[g].nelem() != grids2[g].nelem() )
                        { any_wrong = true; }
                      else
                        {
                          for( Index e=0; e<grids1[g].nelem(); e++ )
                            {
                              const Numeric v1 = grids1[g][e];
                              const Numeric v2 = grids2[g][e];
                              if( ( v1 == 0  &&  abs(v2) > 1e-9 ) || 
                                                 abs(v1-v2)/v1 > 1e-6 )
                                { any_wrong = true; }
                            }   
                        }   
                    }   
                }
              if( any_wrong )
                {
                  ostringstream os;
                  os << "Jacobians for " 
                     << jacobian_quantities2[q2].MainTag()
                     << "/" << jacobian_quantities2[q2].Subtag() 
                     << " shall be appended.\nThis requires that the "
                     << "same grids are used for both measurements,\nbut "
                     << "it seems that this requirement is not met.";
                  throw runtime_error(os.str());
                }
            }
        }

      // Create and fill *jacobian*
      //
      const Index  nrq = jacobian_quantities.nelem();
      const Matrix jacobian1 = jacobian;
      //
      jacobian.resize( n1+n2, jacobian_indices[nrq-1][1]+1 );
      jacobian = 0;
      //
      // Put in old part in top-left corner
      jacobian(Range(0,n1),Range(0,jacobian_indices_copy[nrq1-1][1]+1)) = jacobian1;
      // New parts
      for( Index q2=0; q2<nrq2; q2++ )
        {
          jacobian(Range(n1,n2),Range(jacobian_indices[map_table[q2]][0], 
                                      jacobian_indices[map_table[q2]][1]-
                                      jacobian_indices[map_table[q2]][0]+1) ) =
                jacobian2(joker,Range(jacobian_indices2[q2][0], 
                                      jacobian_indices2[q2][1]-
                                      jacobian_indices2[q2][0]+1) );
        } 
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void yApplyUnit(
         Vector&         y,
         Matrix&         jacobian,
   const Vector&         y_f,
   const ArrayOfIndex&   y_pol,
   const String&         iy_unit,
   const Verbosity&)
{
  if( iy_unit == "1" )
    { throw runtime_error(
        "No need to use this method with *iy_unit* = \"1\"." ); }

  if( max(y) > 1e-3 )
    {
      ostringstream os;
      os << "The spectrum vector *y* is required to have original radiance\n"
         << "unit, but this seems not to be the case. This as a value above\n"
         << "1e-3 is found in *y*.";
      throw runtime_error( os.str() );      
    }

  // Is jacobian set?
  //
  const Index ny = y.nelem();
  //
  const bool do_j = jacobian.nrows() == ny;

  // Some jacobian quantities can not be handled
  if( do_j  &&  max(jacobian) > 1e-3 )
    {
      ostringstream os;
      os << "The method can not be used with jacobian quantities that are not\n"
         << "obtained through radiative transfer calculations. One example on\n"
         << "quantity that can not be handled is *jacobianAddPolyfit*.\n"
         << "The maximum value of *jacobian* indicates that one or several\n"
         << "such jacobian quantities are included.";
      throw runtime_error( os.str() );      
    }

  // Planck-Tb 
  //-------------------------------------------------------------------------- 
  if( iy_unit == "PlanckBT" )
    {
      // Hard to use telescoping here as the data are sorted differently in y
      // and jacobian, than what is expected apply_iy_unit. Copy to temporary
      // variables instead.

      // Handle the elements in "frequency chunks"

      Index i0 = 0;           // Index of first element for present chunk
      //
      while( i0 < ny )
        { 
          // Find number of values for this chunk
          Index n = 1;
          //
          while( i0+n < ny &&  y_f[i0] == y_f[i0+n] ) 
            { n++; }                              

          Matrix yv(1,n);  
          ArrayOfIndex i_pol(n);
          bool any_quv = false;
          //
          for( Index i=0; i<n; i++ )
            { 
              const Index ix=i0+i;   
              yv(0,i) = y[ix];   
              i_pol[i] = y_pol[ix]; 
              if( i_pol[i] > 1  &&  i_pol[i] < 5 )
                { any_quv = true; }
            }

          // Index of elements to convert
          Range ii( i0, n );

          if( do_j )
            {
              if( any_quv  &&  i_pol[0] != 1 )
                {
                  ostringstream os;
                  os << "The conversion to PlanckBT, of the Jacobian and "
                     << "errors for Q, U and V, requires that I (first Stokes "
                     << "element) is at hand and that the data are sorted in "
                     << "such way that I comes first for each frequency.";
                  throw runtime_error( os.str() );      
                }

              // Jacobian
              if( do_j )
                {
                  Tensor3 J(jacobian.ncols(),1,n);
                  J(joker,0,joker) = transpose( jacobian(ii,joker) );
                  apply_iy_unit2( J, yv, iy_unit, y_f[i0], 1, i_pol ); 
                  jacobian(ii,joker) = transpose( J(joker,0,joker) );
                }
            }

          // y (must be done last)
          apply_iy_unit( yv, iy_unit, y_f[i0], 1, i_pol );
          y[ii] = yv(0,joker);

          i0 += n;
        }
    }


  // Other conversions
  //-------------------------------------------------------------------------- 
  else
    {
      // Here we take each element of y separately. 

      Matrix yv(1,1);
      ArrayOfIndex i_pol(1);

      for( Index i=0; i<ny; i++ )
        {
          yv(0,0)  = y[i];    
          i_pol[0] = y_pol[i];

          // Jacobian
          if( do_j )
            { apply_iy_unit2( MatrixView( jacobian(i,joker) ), yv, 
                              iy_unit, y_f[i], 1, i_pol ); }

          // y (must be done last)
          apply_iy_unit( yv, iy_unit, y_f[i], 1, i_pol );
          y[i] = yv(0,0);
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyIndependentBeamApproximation(
         Workspace&        ws,
         Matrix&           iy,
         ArrayOfTensor4&   iy_aux,
         Ppath&            ppath,
         ArrayOfTensor3&   diy_dx,
         GriddedField4&    atm_fields_compact,
   const Index&            iy_id,
   const Vector&           f_grid,
   const Index&            atmosphere_dim,
   const Vector&           p_grid,
   const Vector&           lat_grid,
   const Vector&           lon_grid,
   const Vector&           lat_true,
   const Vector&           lon_true,
   const Tensor3&          t_field,
   const Tensor3&          z_field,
   const Tensor4&          vmr_field,
   const Tensor4&          t_nlte_field,
   const Tensor3&          wind_u_field,
   const Tensor3&          wind_v_field,
   const Tensor3&          wind_w_field,
   const Tensor3&          mag_u_field,
   const Tensor3&          mag_v_field,
   const Tensor3&          mag_w_field,
   const Index&            cloudbox_on,
   const ArrayOfIndex&     cloudbox_limits, 
   const Tensor4&          pnd_field,
   const Matrix&           particle_masses,
   const Agenda&           ppath_agenda,
   const Numeric&          ppath_lmax,
   const Numeric&          ppath_lraytrace,
   const Index&            iy_agenda_call1,
   const String&           iy_unit,  
   const Tensor3&          iy_transmission,
   const Vector&           rte_pos,
   const Vector&           rte_los,
   const Vector&           rte_pos2,
   const Index&            jacobian_do,
   const ArrayOfString&    iy_aux_vars,
   const Agenda&           iy_sub_agenda,
   const Index&            return_atm1d,
   const Index&            skip_vmr,
   const Index&            skip_pnd,
   const Index&            return_masses,
   const Verbosity& )
{
  // Throw error if unsupported features are requested
  if( jacobian_do )
    throw runtime_error( "Jacobians not provided by the method, *jacobian_do* "
                         "must be 0." );
  if( !t_nlte_field.empty() )
    throw runtime_error( "This method does not yet support non-empty *t_nlte_field*." );
  if( !wind_u_field.empty() )
    throw runtime_error( "This method does not yet support non-empty *wind_u_field*." );
  if( !wind_v_field.empty() )
    throw runtime_error( "This method does not yet support non-empty *wind_v_field*." );
  if( !wind_w_field.empty() )
    throw runtime_error( "This method does not yet support non-empty *wind_w_field*." );
  if( !mag_u_field.empty() )
    throw runtime_error( "This method does not yet support non-empty *mag_u_field*." );
  if( !mag_v_field.empty() )
    throw runtime_error( "This method does not yet support non-empty *mag_v_field*." );
  if( !mag_w_field.empty() )
    throw runtime_error( "This method does not yet support non-empty *mag_w_field*." );
  //
  if( return_masses )
    {
      if( pnd_field.nbooks() != particle_masses.nrows() )
        throw runtime_error( "Sizes of *pnd_field* and *particle_masses* "
                               "are inconsistent." );
    }
  

  // Note that input 1D atmospheres are handled exactly as 2D and 3D, to make
  // the function totally general. And 1D must be handled for iterative calls.

  
  // Determine propagation path (with cloudbox deactivated) and check
  // that is OK for ICA
  //
  ppath_agendaExecute( ws, ppath, ppath_lmax, ppath_lraytrace,
                       rte_pos, rte_los, rte_pos2, 0, 0, t_field,
                       z_field, vmr_field, f_grid, ppath_agenda );
  // This check could be improved: 
  bool ppath_is_down = is_decreasing( ppath.r );
  if( !ppath_is_down  &&  !is_increasing( ppath.r ) )
    throw runtime_error( "A propagation path of limb character found. Such "
                         "viewing geometries are not supported by ICA. Propagation "
                         "paths must result in strictly increasing or decreasing "
                         "altitudes." );

  // If scattering and sensor inside atmosphere, we need a pseudo-ppath that
  // samples altitudes not covered by main ppath. We make this second path
  // strictly vertical.
  //
  Ppath   ppath2;
  //
  if( cloudbox_on  &&  ppath.end_lstep == 0 )
    {
      Vector los_tmp = rte_los;
      if( abs(rte_los[0]) < 90 )
        { los_tmp[0] = 180; }
      else
        { los_tmp[0] = 0; }
      //
      ppath_agendaExecute( ws, ppath2, ppath_lmax, ppath_lraytrace,
                           rte_pos, los_tmp, rte_pos2, 0, 0, t_field,
                           z_field, vmr_field, f_grid, ppath_agenda );
    }
  else
    { ppath2.np = 1; }

  // Grid positions, sorted correctly
  const Index     np = ppath.np + ppath2.np - 1;
  ArrayOfGridPos  gp_p(np), gp_lat(np), gp_lon(np);
  if( ppath_is_down )
    {
      // Copy ppath in reversed order
      for( Index i=0; i<ppath.np; i++ )
        {
          const Index ip = ppath.np - i - 1;
          gp_p[i]   = ppath.gp_p[ip];
          if( atmosphere_dim > 1 )
            {
              gp_lat[i] = ppath.gp_lat[ip];
              if( atmosphere_dim == 3 )
                { gp_lon[i] = ppath.gp_lon[ip]; }
            }
        }
      // Append ppath2, but skipping element [0]
      for( Index i=ppath.np; i<np; i++ )
        {
          const Index ip = i - ppath.np + 1;
          gp_p[i]   = ppath2.gp_p[ip];
          if( atmosphere_dim > 1 )
            {
              gp_lat[i] = ppath2.gp_lat[ip];
              if( atmosphere_dim == 3 )
                { gp_lon[i] = ppath2.gp_lon[ip]; }
            }
        }
    }
  else
    {
      // Copy ppath2 in reversed order, but skipping element [0]
      for( Index i=0; i<ppath2.np-1; i++ )
        {
          const Index ip = ppath2.np - i - 1;
          gp_p[i]   = ppath2.gp_p[ip];
          if( atmosphere_dim > 1 )
            {
              gp_lat[i] = ppath2.gp_lat[ip];
              if( atmosphere_dim == 3 )
                { gp_lon[i] = ppath2.gp_lon[ip]; }
            }
        }
      // Append ppath
      for( Index i=ppath2.np-1; i<np; i++ )
        {
          const Index ip = i - ppath2.np + 1;
          gp_p[i]   = ppath.gp_p[ip];
          if( atmosphere_dim > 1 )
            {
              gp_lat[i] = ppath.gp_lat[ip];
              if( atmosphere_dim == 3 )
                { gp_lon[i] = ppath.gp_lon[ip]; }
            }
        }
    }
  

  // 1D version of p_grid 
  Matrix         itw;
  Vector         p1( np );
  ArrayOfGridPos gp0(0), gp1(1);
  interp_atmfield_gp2itw( itw, 1, gp_p, gp0, gp0 );
  itw2p( p1, p_grid, gp_p, itw );

  // 1D version of lat and lon variables
  Vector lat1(0), lon1(0);
  Vector lat_true1(1), lon_true1(1);
  if( atmosphere_dim == 3 )
    {
      gp1[0] = gp_lat[0];
      interp_atmfield_gp2itw( itw, 1, gp1, gp0, gp0 );
      interp( lat_true1, itw, lat_grid, gp1 ); 
      gp1[0] = gp_lon[0];
      interp_atmfield_gp2itw( itw, 1, gp1, gp0, gp0 );
      interp( lon_true1, itw, lon_grid, gp1 );
    }
  else if( atmosphere_dim == 2 )
    {
      gp1[0] = gp_lat[0];
      interp_atmfield_gp2itw( itw, 1, gp1, gp0, gp0 );
      interp( lat_true1, itw, lat_true, gp1 ); 
      interp( lon_true1, itw, lon_true, gp1 );
    }
  else
    {
      lat_true1[0] = lat_true[0];
      lon_true1[0] = lon_true[0];
    }
  
  
  // 2D/3D interpolation weights
  interp_atmfield_gp2itw( itw, atmosphere_dim, gp_p, gp_lat, gp_lon );  
  
  // 1D temperature field
  Tensor3 t1( np, 1, 1 );  
  interp_atmfield_by_itw( t1(joker,0,0), atmosphere_dim, t_field,
                          gp_p, gp_lat, gp_lon, itw )  ;

  // 1D altitude field
  Tensor3 z1( np, 1, 1 );  
  interp_atmfield_by_itw( z1(joker,0,0), atmosphere_dim, z_field,
                          gp_p, gp_lat, gp_lon, itw )  ;
  
  // 1D VMR field
  Tensor4 vmr1( vmr_field.nbooks(), np, 1, 1 );  
  for( Index is=0; is<vmr_field.nbooks(); is++ )
    { interp_atmfield_by_itw( vmr1(is,joker,0,0), atmosphere_dim,
                              vmr_field( is, joker, joker, joker ), 
                              gp_p, gp_lat, gp_lon, itw ); }

  // 1D surface altitude
  Matrix zsurf1(1,1);
  zsurf1(0,0) = z1(0,0,0);
  
  // 1D version of rte_pos/los
  Vector pos1(1); pos1[0] = rte_pos[0]; 
  Vector los1(1); los1[0] = abs( rte_los[0] );
  Vector pos2(0); if( rte_pos2.nelem() ) { pos2 = rte_pos2[Range(0,rte_pos2.nelem())]; } 


  
  // Cloudbox variables
  //
  Index cbox_on1 = cloudbox_on;
  ArrayOfIndex cbox_lims1(0);
  Tensor4 pnd1(0,0,0,0);
  //
  if( cloudbox_on )
    {
      // Determine what p1-levels that are inside cloudbox
      Index ifirst = np;
      Index ilast  = -1;
      for( Index i=0; i<np; i++ )
        {
          if( is_gp_inside_cloudbox( gp_p[i], gp_lat[i], gp_lon[i], 
                                     cloudbox_limits, true, atmosphere_dim ) )
            {
              if( i < ifirst ) { ifirst = i; }
              ilast = i;
            }
        }

      // If no hit, deactive cloudbox
      if( ifirst == np )
        { cbox_on1 = 0; }

      // Otherwise set 1D cloud variables      
      else
        {
          // We can enter the cloudbox from the side, and we need to add 1
          // level on each side to be safe
          //
          const Index extra_bot = ifirst == 0 ? 0 : 1;
          const Index extra_top = ilast == np-1 ? 0 : 1;
          //
          cbox_lims1.resize(2);
          cbox_lims1[0] = ifirst - extra_bot;
          cbox_lims1[1] = ilast + extra_top;
          
          // pnd_field
          //
          pnd1.resize( pnd_field.nbooks(), cbox_lims1[1]-cbox_lims1[0]+1, 1, 1 );
          //
          itw.resize( 1, Index(pow(2.0,Numeric(atmosphere_dim))) );
          //
          for( Index i=extra_bot; i<pnd1.npages()-extra_top; i++ )
            {
              const Index i0 = cbox_lims1[0] + i;
              ArrayOfGridPos gpc_p(1), gpc_lat(1), gpc_lon(1);
              interp_cloudfield_gp2itw( itw(0,joker), 
                                            gpc_p[0], gpc_lat[0], gpc_lon[0], 
                                            gp_p[i0], gp_lat[i0], gp_lon[i0],
                                            atmosphere_dim, cloudbox_limits );
              for( Index p=0; p<pnd_field.nbooks(); p++ )
                {
                  interp_atmfield_by_itw( pnd1(p,i,0,0), atmosphere_dim,
                                          pnd_field(p,joker,joker,joker), 
                                          gpc_p, gpc_lat, gpc_lon, itw );
                }
            }
        }      
    }
  
  // Call sub agenda
  //
  {
    const Index adim1 = 1;
    const Numeric lmax1 = -1;
    Ppath ppath1d;   
    //
    iy_sub_agendaExecute( ws, iy, iy_aux, ppath1d, diy_dx, iy_agenda_call1,
                          iy_unit, iy_transmission, iy_aux_vars, iy_id,
                          f_grid, adim1, p1, lat1, lon1, lat_true1, lon_true1,
                          t1, z1, vmr1, zsurf1, lmax1, ppath_lraytrace,
                          cbox_on1, cbox_lims1, pnd1,
                          jacobian_do, pos1, los1, pos2, iy_sub_agenda );
  }


  // Fill *atm_fields_compact*?
  if( return_atm1d )
    {
      // Sizes and allocate memory
      const Index nvmr  = skip_vmr ? 0 : vmr1.nbooks();
      const Index npnd  = skip_pnd ? 0 : pnd1.nbooks();
      const Index nmass = return_masses ? particle_masses.ncols() : 0;
      const Index ntot  = 2 + nvmr + npnd + nmass;
      ArrayOfString field_names( ntot );
      atm_fields_compact.resize( ntot, np, 1, 1 );

      // Altitudes
      field_names[0] = "Geometric altitudes";
      atm_fields_compact.data(0,joker,0,0) = z1(joker,0,0);

      // Temperature
      field_names[1] = "Temperature";
      atm_fields_compact.data(1,joker,0,0) = t1(joker,0,0);

      // VMRs
      if( nvmr )
        {
          for( Index i=0; i<nvmr; i++ )
            {
              const Index iout = 2 + i;
              ostringstream sstr;
              sstr << "VMR species " << i;
              field_names[iout] = sstr.str();
              atm_fields_compact.data(iout,joker,0,0) = vmr1(i,joker,0,0);
            }
        }

      // PNDs
      if( npnd )
        {
          for( Index i=0; i<npnd; i++ )
            {
              const Index iout = 2 + nvmr + i;
              ostringstream sstr;
              sstr << "Scattering element " << i;
              field_names[iout] = sstr.str();
              atm_fields_compact.data(iout,joker,0,0) = 0;
              atm_fields_compact.data(iout,Range(cbox_lims1[0],pnd1.npages()),0,0) =
                pnd1(i,joker,0,0);
            }
        }

      // Masses
      if( nmass )
        {
          for( Index i=0; i<nmass; i++ )
            {
              const Index iout = 2 + nvmr + npnd + i;
              ostringstream sstr;
              sstr << "Mass category " << i;
              field_names[iout] = sstr.str();
              atm_fields_compact.data(iout,joker,0,0) = 0;
              for( Index ip=cbox_lims1[0]; ip<pnd1.npages(); ip++ )
                {
                  for( Index is=0; is<pnd1.nbooks(); is++ )
                    {
                      atm_fields_compact.data(iout,ip,0,0) +=
                        particle_masses(is,i) * pnd1(is,ip,0,0);
                    }
                }
            }
        }

      // Finally, set grids and names
      //
      atm_fields_compact.set_name( "Data created by *iyIndependentBeamApproximation*" );
      //
      atm_fields_compact.set_grid_name( GFIELD4_FIELD_NAMES, "Atmospheric quantity" );
      atm_fields_compact.set_grid( GFIELD4_FIELD_NAMES, field_names );  
      atm_fields_compact.set_grid_name( GFIELD4_P_GRID, "Pressure" );
      atm_fields_compact.set_grid( GFIELD4_P_GRID, p1 );
      atm_fields_compact.set_grid_name( GFIELD4_LAT_GRID, "Latitude" );
      atm_fields_compact.set_grid( GFIELD4_LAT_GRID, lat_true1 );
      atm_fields_compact.set_grid_name( GFIELD4_LON_GRID, "Longitude" );
      atm_fields_compact.set_grid( GFIELD4_LON_GRID, lon_true1 );
    }
}
