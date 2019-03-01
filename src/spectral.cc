/* Copyright (C) 2018-2020 Jakob Doerr <jakobdoerr@googlemail.com>

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
  \file   spectral.cc
  \author Jakob Doerr <jakobdoerr@googlemail.com>
  \date   Fri May 16 11:29:59 2018

  \brief  This file contains definitions and functions related to the
          spectral transformation.


 */

#include <cfloat>
#include <cmath>
#include <stdexcept>
#include "array.h"
#include "arts.h"
#include "check_input.h"
#include "interpolation.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "messages.h"
#include "optproperties.h"
#include "xml_io.h"
#include "spectral.h"
#include "shtns.h"
#include <fftw3.h>
#include <complex.h>

extern const Numeric DEG2RAD;
//! one-line descript
/*!
  Descript/Doc

  \param[out] name  desc.
  \param[in]  name  desc.

  \author Jana Mendrok
  \date   2018-01-15
*/
/*
void methodname(//Output
                type& name,
                //Input
                const type& name)
{
}


*/

//! Transforms spectral extinction matrix and absorption vector to grid

/*!
  Transforms given spectral extinction matrix and absorption vector onto
  a discrete grid provided as input using SHTns

  \param[out] ext_matrix  The extinction matrix on the discrete angular grid
  \param[out]  abs_vector  The absorption vector on the discrete angular grid

  \param[in]  ext_matrix_spectral  The extinction matric on the spectral grid
  \param[in]  abs_vector_spectral  The absorption vector on the spectral grid
  \param[in]  dir_array  The discrete angular grid (as pairs of theta and phi)
  \param[in]  any_m  Flag wether any m coefficients are present


  \author Jakob Doerr
  \date   2018-05-25
*/

void opt_prop_SpecToGrid(//Output
                Tensor5View ext_matrix,
                Tensor4View abs_vector,
                //Input
                const Tensor5& ext_matrix_spectral,
                const Tensor4& abs_vector_spectral,
                const Matrix& dir_array,
                const bool& any_m) {
  // Initialize gridded fields
  //cout << abs_vector_spectral(1,0,joker,0) << "\n";
  ext_matrix = 0.;
  abs_vector = 0.;
  shtns_cfg shtns;                // handle to a sht transform configuration
  complex<double> *Qlm;      // spherical harmonics coefficients (l,m space): complex numbers.
  int lmax = (int) ext_matrix_spectral.npages() - 1;
  int nlat = 40;
  int nphi = 10;
  int mres = 1;
  int mmax,NLM;
  if (any_m)
    mmax = lmax;
  else
    mmax = 0;
  shtns_verbose(0);                       // displays informations during initialization.
  shtns_use_threads(0);           // enable multi-threaded transforms (if supported).
  shtns = shtns_init(sht_gauss, lmax, mmax, mres, nlat, nphi);
  //shtns = shtns_create(lmax, mmax, mres, sht_orthonormal | SHT_REAL_NORM);
  //shtns_set_grid(shtns, sht_gauss, 0.0, nlat, nphi);
  NLM = shtns->nlm;
  // Memory allocation : the use of fftw_malloc is required because we need proper 16-byte alignement.
  // allocate spatial fields.
  // allocate SH representations.
  Qlm = (complex<double> *) fftw_malloc(NLM * sizeof(complex<double>));
  for (Index f_index = 0; f_index < ext_matrix_spectral.nshelves(); f_index++)
  {
    for (Index p_index = 0; p_index < ext_matrix_spectral.nbooks(); p_index++)
    {
      for (Index nst1 = 0; nst1 < ext_matrix_spectral.nrows(); nst1++)
      {
        for (Index nst2 = 0; nst2 < ext_matrix_spectral.ncols(); nst2++)
        {
          for (Index i_lm = 0; i_lm < NLM; i_lm++)
          {
            Qlm[i_lm] = (complex<double>) ext_matrix_spectral(
                           f_index, p_index, i_lm, nst1, nst2);
          }
          for (Index dind = 0; dind < dir_array.nrows(); dind++)
          {
            ext_matrix(f_index, p_index, dind, nst1, nst2) =
                        SH_to_point(shtns, Qlm, cos(dir_array(dind, 0)*DEG2RAD), dir_array(dind, 1)*DEG2RAD);
          }
        }
      }
      for (Index nst1 = 0; nst1 < abs_vector_spectral.ncols(); nst1++)
      {
        for (Index i_lm = 0; i_lm < lmax; i_lm++)
        {
          Qlm[i_lm] = (complex<double>) abs_vector_spectral(
                      f_index, p_index, i_lm, nst1);
        }
        for (Index dind = 0; dind < dir_array.nrows(); dind++)
        {
          abs_vector(f_index, p_index, dind, nst1) =
                    SH_to_point(shtns, Qlm, cos(dir_array(dind, 0)*DEG2RAD), dir_array(dind, 1)*DEG2RAD);
        }
      }
    }
  }
}


//! Transforms spectral phase matrix to grid

/*!
  Transforms given spectral phase matrix onto
  a discrete grid provided as input using SHTns

  \param[out] ext_matrix  The extinction matrix on the discrete angular grid
  \param[out]  abs_vector  The absorption vector on the discrete angular grid

  \param[in]  ext_matrix_spectral  The extinction matric on the spectral grid
  \param[in]  abs_vector_spectral  The absorption vector on the spectral grid
  \param[in]  dir_array  The discrete angular grid (as pairs of theta and phi)
  \param[in]  any_m  Flag wether any m coefficients are present


  \author Jakob Doerr
  \date   2018-07-25
*/

void pha_mat_SpecToGrid(//Output
        Tensor6View pha_matrix,
        //Input
        const Tensor6& pha_mat_real_spectral,
        const Tensor6& pha_mat_imag_spectral,
        const Matrix& pdir_array,
        const Matrix& idir_array,
        const Index& ptype,
        const bool& any_m_inc,
        const bool& any_m_sca)
{
  //if (pha_mat_real_spectral(0,0,0,0,0,0) == 0)
  //{
  //  ostringstream os;
  //  os << "Looks like you have all zeros in your data....dumbass!";
  //  throw runtime_error( os.str() );
  //}
  // Initialize gridded fields
  pha_matrix = 0.;
  shtns_cfg shtns_inc, shtns_sca;                // handle to a sht transform configuration
  complex<double> *Qlm_inc, *Qlm_sca;      // spherical harmonics coefficients (l,m space): complex numbers.

  int mmax_inc, mmax_sca, lmax_inc, lmax_sca, NLM_inc, NLM_sca;
  // Calculate the number of l and m for incoming and scattered directions from the data

  int ncoeff_inc = (int) pha_mat_real_spectral.npages();
  int ncoeff_sca = (int) pha_mat_real_spectral.nbooks();

  // Incoming
  if (any_m_sca)
  {
    lmax_sca = (int) (- 1.5 + sqrt(0.25 + 2 * ncoeff_sca));  // See doc.
    mmax_sca = lmax_sca;
  }
  else
  {
    lmax_sca = ncoeff_sca - 1;
    mmax_sca = 0;
  }

  // Scattered
  if (any_m_inc)
  {
    lmax_inc = (int) sqrt(ncoeff_inc) - 1;
    mmax_inc = lmax_inc;
  }
  else
  {
    lmax_inc = ncoeff_inc - 1;
    mmax_inc = 0;
  }

  int nlat_inc = lmax_inc + 2 + lmax_inc%2; // Nlat has to be bigger than lmax
  int nlat_sca = lmax_sca + 2 + lmax_sca%2;
  int nphi_inc = 2*mmax_inc + 2 + mmax_inc%2; // NPhi has to be bigger than 2*mmax
  int nphi_sca = 2*mmax_sca + 2 + mmax_sca%2;
  int mres = 1;
  shtns_verbose(0);                       // displays informations during initialization.
  shtns_use_threads(0);           // enable multi-threaded transforms (if supported).
  shtns_inc = shtns_init(sht_gauss, lmax_inc, mmax_inc, mres, nlat_inc, nphi_inc);
  shtns_sca = shtns_init(sht_gauss, lmax_sca, mmax_sca, mres, nlat_sca, nphi_sca);

  //shtns = shtns_create(lmax, mmax, mres, sht_orthonormal | SHT_REAL_NORM);
  //shtns_set_grid(shtns, sht_gauss, 0.0, nlat, nphi);
  NLM_inc = shtns_inc->nlm;
  NLM_sca = shtns_sca->nlm;
  // Memory allocation : the use of fftw_malloc is required because we need proper 16-byte alignement.
  // allocate spatial fields.
  // allocate SH representations.
  Qlm_inc = (complex<double> *) fftw_malloc(NLM_inc * sizeof(complex<double>));
  Qlm_sca = (complex<double> *) fftw_malloc(NLM_sca * sizeof(complex<double>));

  const Index nf = pha_matrix.nvitrines();
  const Index nT = pha_matrix.nshelves();
  const Index npDir = pdir_array.nrows();
  assert( pha_matrix.nbooks() == npDir );

  const Index stokes_dim = pha_matrix.ncols();
  assert( pha_matrix.nrows() == stokes_dim );

  const Index niDir = idir_array.nrows();
  assert( pha_matrix.npages() == niDir );

  // FIXME: Get rid of this assertion later, b\c it can handle any ptype..
  assert( ptype == PTYPE_TOTAL_RND or ptype == PTYPE_AZIMUTH_RND );

  if ( ptype == PTYPE_TOTAL_RND)
  {
    for( Index pdir=0; pdir<npDir; pdir++ )
    {
      for (Index idir = 0; idir < niDir; idir++)
      {
        // calc scat ang theta from incident and prop dirs
        Numeric theta = scat_angle(pdir_array(pdir,0), pdir_array(pdir,1),
                                   idir_array(idir,0), idir_array(idir,1));
        for (Index find = 0; find < nf; find++)
        {
          for (Index Tind = 0; Tind < nT; Tind++)
          {
            Vector pha_mat_int(pha_matrix.nrows()*pha_matrix.nrows());
            for (Index st1 = 0; st1 < pha_matrix.nrows(); st1++)
            {
              for (Index st2 = 0; st2 < pha_matrix.ncols(); st2++)
              {

                for (Index p_lm = 0; p_lm < NLM_sca; p_lm++)
                {
                  for (Index i_lm = 0; i_lm < NLM_inc; i_lm++)
                  {
                    // First spectral transformation (incoming anlges)
                    Qlm_inc[i_lm] = complex<double>(pha_mat_real_spectral(
                            find, Tind, p_lm, i_lm, st1, st2)
                                    ,pha_mat_imag_spectral(
                            find, Tind, p_lm, i_lm, st1, st2));
                  }
                  // Here the first trafo happens. We can take the angle
                  // theta=0, phi=0 since there is no incident angle
                  // dependence for totally random particles
                  Qlm_sca[p_lm] = SH_to_point(shtns_inc,Qlm_inc,
                            0,0);
                }
                // Second spectral transformation (scattered angles)
                pha_mat_int[st1*4+st2] = SH_to_point(shtns_sca,Qlm_sca,cos(theta*DEG2RAD),0);
              }
            }
            pha_mat_labCalc(pha_matrix(find,Tind,pdir,idir,joker,joker),pha_mat_int,
                            pdir_array(pdir,0), pdir_array(pdir,1),
                            idir_array(idir,0), idir_array(idir,1),
                            theta);
          }
        }
      }
    }
  } // FIXME: IS IT RIGHT LIKE THIS???                                                                                                                                                                                                                                                                                                           
  else // ptype azi random (and general)
  {
    Index nDir = npDir * niDir;
    Matrix delta_aa(npDir, niDir);

    for (Index pdir = 0; pdir < npDir; pdir++)
    {
      for (Index idir = 0; idir < niDir; idir++)
      {
        delta_aa(pdir, idir) = pdir_array(pdir, 1) - idir_array(idir, 1) +
                               (pdir_array(pdir, 1) - idir_array(idir, 1) < -180) * 360 -
                               (pdir_array(pdir, 1) - idir_array(idir, 1) > 180) * 360;
      }
    }
    for (Index find = 0; find < nf; find++)
    {
      for (Index Tind = 0; Tind < nT; Tind++)
      {
        for (Index pdir = 0; pdir < npDir; pdir++)
        {
          for (Index idir = 0; idir < niDir; idir++)
          {
            for (Index st1 = 0; st1 < pha_matrix.nrows(); st1++)
            {
              for (Index st2 = 0; st2 < pha_matrix.ncols(); st2++)
              {
                for (Index p_lm = 0; p_lm < NLM_sca; p_lm++)
                {
                  for (Index i_lm = 0; i_lm < NLM_inc; i_lm++)
                  {
                    // First spectral transformation (incoming angles)
                    Qlm_inc[i_lm] = complex<double>(pha_mat_real_spectral(
                                     find, Tind, p_lm, i_lm, st1, st2)
                                          ,pha_mat_imag_spectral(
                                     find, Tind, p_lm, i_lm, st1, st2));
                  }
                  // Here the first trafo happens.
                  Qlm_sca[p_lm] = SH_to_point(shtns_inc, Qlm_inc,
                                              cos(idir_array(idir, 0)*DEG2RAD),
                                              idir_array(idir, 1)*DEG2RAD);
                }
                pha_matrix(find, Tind, pdir, idir, st1, st2) =
                          SH_to_point(shtns_sca, Qlm_sca, cos(pdir_array(pdir, 0)*DEG2RAD),
                                       pdir_array(pdir, 1)*DEG2RAD);
              }
            }
            if( stokes_dim>2 )
            {
              if (delta_aa(pdir, idir) < 0.)
              {
                pha_matrix(find, joker, pdir, idir, 0, 2) *= -1;
                pha_matrix(find, joker, pdir, idir, 1, 2) *= -1;
                pha_matrix(find, joker, pdir, idir, 2, 0) *= -1;
                pha_matrix(find, joker, pdir, idir, 2, 1) *= -1;
              }
            }
            if( stokes_dim>3 )
            {
              if( delta_aa(pdir,idir)<0. )
              {
                pha_matrix(find,joker,pdir,idir,0,3) *= -1;
                pha_matrix(find,joker,pdir,idir,1,3) *= -1;
                pha_matrix(find,joker,pdir,idir,3,0) *= -1;
                pha_matrix(find,joker,pdir,idir,3,1) *= -1;
              }
            }
          }
        }
      }
    }
  }
}



