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
  Descript/Doc

  \param[out] name  desc.
  \param[in]  name  desc.

  \author Jakob Doerr
  \date   2018-05-25
*/

void par_optpropSpecToGrid(//Output
                Tensor6View ext_matrix,
                Tensor5View abs_vector,
                //Input
                const Tensor5& ext_matrix_spectral,
                const Tensor4& abs_vector_spectral,
                const Vector& za_grid,
                const Vector& aa_grid)
{
    for (Index f_index = 0; f_index < ext_matrix_spectral.nshelves(); f_index++)
    {
        for (Index p_index = 0; p_index < ext_matrix_spectral.nbooks(); p_index++)
        {
            for (Index nst1 = 0; nst1 < ext_matrix_spectral.ncols(); nst1++)
            {
                spec_to_grid(abs_vector(f_index,p_index,joker,joker,nst1),
                             abs_vector_spectral(f_index,p_index,joker,nst1),
                             za_grid, aa_grid);

                for (Index nst2 = 0; nst2 < ext_matrix_spectral.nrows(); nst2++)
                {
                    spec_to_grid(ext_matrix(f_index,p_index,joker,joker,nst1,nst2),
                                 ext_matrix_spectral(f_index,p_index,joker,nst1,nst2),
                                 za_grid, aa_grid);

                }
            }
        }
    }
}


//! Spectral transformation Spectral ---> Grid
/*!
  Descript/Doc

  \param[out] name  desc.
  \param[in]  name  desc.

  \author Jakob Doerr
  \date   2018-05-25
*/

void spec_to_grid(//Output
                MatrixView grid_field,
                //Input
                const Vector& spec_field,
                const Vector& za_grid,
                const Vector& aa_grid)
{
    shtns_cfg shtns;                // handle to a sht transform configuration
    long int lmax,mmax,nlat,nphi,mres, NLM;
    //complex double *Slm, *Tlm;      // spherical harmonics coefficients (l,m space): complex numbers.
    double *Sh, *Th;                // real space : theta,phi
    long int i,im,lm;
    double t;
    lmax = spec_field.nelem();       nlat = 32;
    mmax = 0;       nphi = 10;
    mres = 1;
    shtns_verbose(0);                       // displays informations during initialization.
    shtns_use_threads(0);           // enable multi-threaded transforms (if supported).
    //shtns = shtns_init( sht_gauss, lmax, mmax, mres, nlat, nphi );
      shtns = shtns_create(lmax, mmax, mres, sht_orthonormal | SHT_REAL_NORM);
      shtns_set_grid(shtns, za_grid, 0.0, nlat, nphi);
    NLM = shtns->nlm;

    // Memory allocation : the use of fftw_malloc is required because we need proper 16-byte alignement.
// allocate spatial fields.
    //Sh = (double *) fftw_malloc( NSPAT_ALLOC(shtns) * sizeof(double));

// allocate SH representations.
    //Slm = (complex double *) fftw_malloc( NLM * sizeof(complex double));
    //Tlm = (complex double *) fftw_malloc( NLM * sizeof(complex double));

    //SH_to_spat(shtns, spec_field,Sh);
}




