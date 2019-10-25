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

#include "spectral.h"
#include <complex.h>
#include <fftw3.h>
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
#include "shtns.h"
#include "xml_io.h"

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

void opt_prop_SpecToGrid(  //Output
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
  shtns_cfg shtns;  // handle to a sht transform configuration
  complex<double>*
      Qlm;  // spherical harmonics coefficients (l,m space): complex numbers.
  int lmax = (int)ext_matrix_spectral.npages() - 1;
  int nlat = 40;
  int nphi = 10;
  int mres = 1;
  int mmax, NLM;
  if (any_m)
    mmax = lmax;
  else
    mmax = 0;
  shtns_verbose(0);      // displays informations during initialization.
  shtns_use_threads(0);  // enable multi-threaded transforms (if supported).
  shtns = shtns_init(sht_gauss, lmax, mmax, mres, nlat, nphi);
  //shtns = shtns_create(lmax, mmax, mres, sht_orthonormal | SHT_REAL_NORM);
  //shtns_set_grid(shtns, sht_gauss, 0.0, nlat, nphi);
  NLM = shtns->nlm;
  // Memory allocation : the use of fftw_malloc is required because we need proper 16-byte alignement.
  // allocate spatial fields.
  // allocate SH representations.
  Qlm = (complex<double>*)fftw_malloc(NLM * sizeof(complex<double>));
  for (Index f_index = 0; f_index < ext_matrix_spectral.nshelves(); f_index++) {
    for (Index p_index = 0; p_index < ext_matrix_spectral.nbooks(); p_index++) {
      for (Index nst1 = 0; nst1 < ext_matrix_spectral.nrows(); nst1++) {
        for (Index nst2 = 0; nst2 < ext_matrix_spectral.ncols(); nst2++) {
          for (Index i_lm = 0; i_lm < NLM; i_lm++) {
            Qlm[i_lm] = (complex<double>)ext_matrix_spectral(
                f_index, p_index, i_lm, nst1, nst2);
          }
          for (Index dind = 0; dind < dir_array.nrows(); dind++) {
            ext_matrix(f_index, p_index, dind, nst1, nst2) =
                SH_to_point(shtns,
                            Qlm,
                            cos(dir_array(dind, 0) * DEG2RAD),
                            dir_array(dind, 1) * DEG2RAD);
          }
        }
      }
      for (Index nst1 = 0; nst1 < abs_vector_spectral.ncols(); nst1++) {
        for (Index i_lm = 0; i_lm < lmax; i_lm++) {
          Qlm[i_lm] = (complex<double>)abs_vector_spectral(
              f_index, p_index, i_lm, nst1);
        }
        for (Index dind = 0; dind < dir_array.nrows(); dind++) {
          abs_vector(f_index, p_index, dind, nst1) =
              SH_to_point(shtns,
                          Qlm,
                          cos(dir_array(dind, 0) * DEG2RAD),
                          dir_array(dind, 1) * DEG2RAD);
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
  \param[in]  any_m  Flag whether any m coefficients are present


  \author Jakob Doerr
  \date   2018-07-25
*/

void pha_mat_SpecToGrid(  //Output
    Tensor6View pha_matrix,
    //Input
    const Tensor6& pha_mat_real_spectral,
    const Tensor6& pha_mat_imag_spectral,
    const Matrix& pdir_array,
    const Matrix& idir_array,
    const Index& ptype,
    const bool& any_m_inc,
    const bool& any_m_sca) {
  // Initialize gridded fields
  pha_matrix = 0.;
  shtns_cfg shtns_inc, shtns_sca;  // handle to a sht transform configuration
  complex<double>*Qlm_inc, *Qlm_sca,
      *Sh;  // spherical harmonics coefficients (l,m space): complex numbers.
  double* Th;

  int mmax_inc, mmax_sca, lmax_inc, lmax_sca, NLM_inc, NLM_sca;
  // Calculate the number of l and m for incoming and scattered directions from the data

  int ncoeff_inc = (int)pha_mat_real_spectral.npages();
  int ncoeff_sca = (int)pha_mat_real_spectral.nbooks();

  // Incoming
  if (any_m_sca) {
    lmax_sca = (int)(-1.5 + sqrt(0.25 + 2 * ncoeff_sca));  // See doc.
    mmax_sca = lmax_sca;
  } else {
    lmax_sca = ncoeff_sca - 1;
    mmax_sca = 0;
  }

  // Scattered
  if (any_m_inc) {
    lmax_inc = (int)sqrt(ncoeff_inc) - 1;
    mmax_inc = lmax_inc;
  } else {
    lmax_inc = ncoeff_inc - 1;
    mmax_inc = 0;
  }

  int nlat_inc =
      2 * lmax_inc + 2;  //+ lmax_inc%2; // Nlat has to be bigger than lmax
  int nlat_sca = 2 * lmax_sca + 2;  //+ lmax_sca%2;
  int nphi_inc = 2 * mmax_inc + 1;  // NPhi has to be bigger than 2*mmax
  int nphi_sca = 2 * mmax_sca + 1;
  int mres = 1;
  shtns_verbose(0);      // displays information during initialization.
  shtns_use_threads(0);  // enable multi-threaded transforms (if supported).

  // Create the two shtns objects for the two transformations
  shtns_inc =
      shtns_init(sht_reg_poles, lmax_inc, mmax_inc, mres, nlat_inc, nphi_inc);
  shtns_sca =
      shtns_init(sht_reg_poles, lmax_sca, mmax_sca, mres, nlat_sca, nphi_sca);

  NLM_inc = shtns_inc->nlm;
  NLM_sca = shtns_sca->nlm;
  // Memory allocation : the use of fftw_malloc is required because we need proper 16-byte alignement.
  // allocate spatial fields.
  // allocate SH representations.
  Qlm_inc = (complex<double>*)fftw_malloc(NLM_inc * sizeof(complex<double>));
  Qlm_sca = (complex<double>*)fftw_malloc(NLM_sca * sizeof(complex<double>));
  // Allocate the variable for the product after the first transformation
  Sh = (complex<double>*)fftw_malloc(NSPAT_ALLOC(shtns_inc) *
                                     sizeof(complex<double>));

  // Allocate the variable for the product after the second transformation
  Th = (double*)fftw_malloc(NSPAT_ALLOC(shtns_sca) * sizeof(double));

  // Extract the native grids of the two shtns objects
  Vector za_grid_inc(nlat_inc, 0);
  for (Index lat_inc = 0; lat_inc < nlat_inc; lat_inc++) {
    za_grid_inc[lat_inc] = acos(shtns_inc->ct[lat_inc]) / DEG2RAD;
  }
  Vector za_grid_sca(nlat_sca, 0);
  for (Index lat_sca = 0; lat_sca < nlat_sca; lat_sca++) {
    za_grid_sca[lat_sca] = acos(shtns_sca->ct[lat_sca]) / DEG2RAD;
  }
  Vector aa_grid_sca(nphi_sca, 0);
  for (Index phi_sca = 0; phi_sca < nphi_sca; phi_sca++) {
    aa_grid_sca[phi_sca] = 360. / nphi_sca * (double)phi_sca;
  }

  const Index nf = pha_matrix.nvitrines();
  const Index nT = pha_matrix.nshelves();
  const Index npDir = pdir_array.nrows();
  assert(pha_matrix.nbooks() == npDir);

  const Index stokes_dim = pha_matrix.ncols();
  assert(pha_matrix.nrows() == stokes_dim);

  const Index niDir = idir_array.nrows();
  assert(pha_matrix.npages() == niDir);

  // FIXME: Get rid of this assertion later, b\c it can handle any ptype..
  assert(ptype == PTYPE_TOTAL_RND or ptype == PTYPE_AZIMUTH_RND);

  if (ptype == PTYPE_TOTAL_RND) {
    for (Index pdir = 0; pdir < npDir; pdir++) {
      for (Index idir = 0; idir < niDir; idir++) {
        // calc scat ang theta from incident and prop dirs
        Numeric theta = scat_angle(pdir_array(pdir, 0),
                                   pdir_array(pdir, 1),
                                   idir_array(idir, 0),
                                   idir_array(idir, 1));
        for (Index find = 0; find < nf; find++) {
          for (Index Tind = 0; Tind < nT; Tind++) {
            Vector pha_mat_int(pha_matrix.nrows() * pha_matrix.nrows());
            for (Index st1 = 0; st1 < pha_matrix.nrows(); st1++) {
              for (Index st2 = 0; st2 < pha_matrix.ncols(); st2++) {
                for (Index p_lm = 0; p_lm < NLM_sca; p_lm++) {
                  for (Index i_lm = 0; i_lm < NLM_inc; i_lm++) {
                    // Prepare first spectral transformation (incoming angles)
                    Qlm_inc[i_lm] = complex<double>(
                        pha_mat_real_spectral(find, Tind, p_lm, i_lm, st1, st2),
                        pha_mat_imag_spectral(
                            find, Tind, p_lm, i_lm, st1, st2));
                  }
                  // Here the first trafo happens. We can take the angle
                  // theta=0, phi=0 since there is no incident angle
                  // dependence for totally random particles
                  Qlm_sca[p_lm] = SH_to_point(shtns_inc, Qlm_inc, 0, 0);
                }
                // Second spectral transformation (scattered angles)
                pha_mat_int[st1 * 4 + st2] =
                    SH_to_point(shtns_sca, Qlm_sca, cos(theta * DEG2RAD), 0);
              }
            }
            pha_mat_labCalc(pha_matrix(find, Tind, pdir, idir, joker, joker),
                            pha_mat_int,
                            pdir_array(pdir, 0),
                            pdir_array(pdir, 1),
                            idir_array(idir, 0),
                            idir_array(idir, 1),
                            theta);
          }
        }
      }
    }
  }

  else if (ptype == PTYPE_AZIMUTH_RND)  // ptype azi random (and general)
  {
    ComplexMatrix pha_mat_singletrans(NLM_sca, nlat_inc * nphi_inc);
    Matrix pha_matrix_temp(nlat_sca * nphi_sca, nlat_inc * nphi_inc);
    Matrix delta_aa(npDir, niDir);

    // Variables for the trilinear interpolation (phi_sca angles)
    Index nDir = npDir * niDir;
    Tensor4 pha_matrix_full_temp(nlat_sca, nphi_sca, nlat_inc, nphi_inc);
    Vector pha_mat_int(nDir);
    Matrix dir_itw(nDir, 8);

    ArrayOfGridPos daa_gp(nDir), pza_gp(nDir), iza_gp(nDir);
    ArrayOfGridPos pza_gp_tmp(npDir), iza_gp_tmp(niDir);

    // Variables for the double linear interpolation (no phi_sca angles)
    Matrix pha_matrix_interp(npDir, nlat_inc);
    Matrix idir_itw(niDir, 2);
    Matrix pdir_itw(npDir, 2);

    ArrayOfGridPos idir_gp(niDir);
    ArrayOfGridPos pdir_gp(npDir);

    if (any_m_sca) {
      Vector adelta_aa(nDir);

      gridpos(pza_gp_tmp, za_grid_sca, pdir_array(joker, 0));
      gridpos(iza_gp_tmp, za_grid_inc, idir_array(joker, 0));

      Index j = 0;
      for (Index pdir = 0; pdir < npDir; pdir++) {
        for (Index idir = 0; idir < niDir; idir++) {
          delta_aa(pdir, idir) =
              pdir_array(pdir, 1) - idir_array(idir, 1) +
              (pdir_array(pdir, 1) - idir_array(idir, 1) < -180) * 360 -
              (pdir_array(pdir, 1) - idir_array(idir, 1) > 180) * 360;
          adelta_aa[j] = abs(delta_aa(pdir, idir));
          pza_gp[j] = pza_gp_tmp[pdir];
          iza_gp[j] = iza_gp_tmp[idir];
          j++;
        }
        gridpos(daa_gp, aa_grid_sca, adelta_aa);

        interpweights(dir_itw, pza_gp, daa_gp, iza_gp);
      }
    } else {
      // derive the direction interpolation weights for incoming angles.

      gridpos(idir_gp, za_grid_inc, idir_array(joker, 0));
      // only interpolating in za, ie 1D linear interpol

      interpweights(idir_itw, idir_gp);
      // derive the direction interpolation weights for scattered angles.

      gridpos(pdir_gp, za_grid_sca, pdir_array(joker, 0));
      // only interpolating in za, ie 1D linear interpol
      interpweights(pdir_itw, pdir_gp);

      for (Index pdir = 0; pdir < npDir; pdir++) {
        for (Index idir = 0; idir < niDir; idir++) {
          delta_aa(pdir, idir) =
              pdir_array(pdir, 1) - idir_array(idir, 1) +
              (pdir_array(pdir, 1) - idir_array(idir, 1) < -180) * 360 -
              (pdir_array(pdir, 1) - idir_array(idir, 1) > 180) * 360;
        }
      }
    }
    // --------

    for (Index find = 0; find < nf; find++) {
      for (Index Tind = 0; Tind < nT; Tind++) {
        for (Index st1 = 0; st1 < pha_matrix.nrows(); st1++) {
          for (Index st2 = 0; st2 < pha_matrix.ncols(); st2++) {
            for (Index p_lm = 0; p_lm < NLM_sca; p_lm++) {
              for (Index i_lm = 0; i_lm < NLM_inc; i_lm++) {
                // create the complex coefficients array from the real and
                // imaginary part
                Qlm_inc[i_lm] = complex<double>(
                    pha_mat_real_spectral(find, Tind, p_lm, i_lm, st1, st2),
                    pha_mat_imag_spectral(find, Tind, p_lm, i_lm, st1, st2));
              }
              // Here the first transformation happens.
              SH_to_spat_cplx(shtns_inc, Qlm_inc, Sh);
              for (Index idir = 0; idir < nlat_inc; idir++) {
                pha_mat_singletrans(p_lm, idir) = Sh[idir];
              }
            }
            for (Index idir = 0; idir < nlat_inc; idir++) {
              for (Index p_lm = 0; p_lm < NLM_sca; p_lm++) {
                Qlm_sca[p_lm] = pha_mat_singletrans(p_lm, idir);
              }
              // Here the second trafo happens
              SH_to_spat(shtns_sca, Qlm_sca, Th);
              for (Index pdir = 0; pdir < nlat_sca; pdir++) {
                pha_matrix_temp(pdir, idir) = Th[pdir];
              }
            }

            if (any_m_sca) {
              // Blow up the temporary phase matrix, so that it has 4 dimensions (za_sca, aa_sca, za_inc, aa_inc)
              for (Index idx_za_sca = 0; idx_za_sca < nlat_sca; idx_za_sca++) {
                for (Index idx_aa_sca = 0; idx_aa_sca < nphi_sca;
                     idx_aa_sca++) {
                  for (Index idx_za_inc = 0; idx_za_inc < nlat_inc;
                       idx_za_inc++) {
                    for (Index idx_aa_inc = 0; idx_aa_inc < nphi_inc;
                         idx_aa_inc++) {
                      pha_matrix_full_temp(
                          idx_za_sca, idx_aa_sca, idx_za_inc, idx_aa_inc) =
                          pha_matrix_temp(idx_aa_sca * nlat_sca + idx_za_sca,
                                          idx_aa_inc * nlat_inc + idx_za_inc);
                    }
                  }
                }
              }
              // perform the (tri-linear) angle interpolation. but only for the
              // pha_mat elements that we actually need.

              interp(pha_mat_int,
                     dir_itw,
                     pha_matrix_full_temp(joker, joker, joker, 0),
                     pza_gp,
                     daa_gp,
                     iza_gp);

              // sort direction-combined 1D-array back into prop and incident
              // direction matrix
              Index i = 0;
              for (Index pdir = 0; pdir < npDir; pdir++)
                for (Index idir = 0; idir < niDir; idir++) {
                  pha_matrix(find, Tind, pdir, idir, st1, st2) = pha_mat_int[i];
                  i++;
                }
            } else {
              for (Index idir = 0; idir < nlat_inc; idir++) {
                // Interpolate the temporary phase matrix from the native SHTNS
                // angular grids to the grid specified in the input of this method
                interp(pha_matrix_interp(joker, idir),
                       pdir_itw,
                       pha_matrix_temp(joker, idir),
                       pdir_gp);
              }
              for (Index pdir = 0; pdir < npDir; pdir++) {
                // Interpolate the temporary phase matrix from the native SHTNS
                // angular grids to the grid specified in the input of this method
                interp(pha_matrix(find, Tind, pdir, joker, st1, st2),
                       idir_itw,
                       pha_matrix_interp(pdir, joker),
                       idir_gp);
              }
            }
          }
        }
      }
      for (Index pdir = 0; pdir < npDir; pdir++) {
        for (Index idir = 0; idir < niDir; idir++) {
          if (stokes_dim > 2) {
            if (delta_aa(pdir, idir) < 0.) {
              pha_matrix(find, joker, pdir, idir, 0, 2) *= -1;
              pha_matrix(find, joker, pdir, idir, 1, 2) *= -1;
              pha_matrix(find, joker, pdir, idir, 2, 0) *= -1;
              pha_matrix(find, joker, pdir, idir, 2, 1) *= -1;
            }
          }
          if (stokes_dim > 3) {
            if (delta_aa(pdir, idir) < 0.) {
              pha_matrix(find, joker, pdir, idir, 0, 3) *= -1;
              pha_matrix(find, joker, pdir, idir, 1, 3) *= -1;
              pha_matrix(find, joker, pdir, idir, 3, 0) *= -1;
              pha_matrix(find, joker, pdir, idir, 3, 1) *= -1;
            }
          }
        }
      }
    }
  }
}
