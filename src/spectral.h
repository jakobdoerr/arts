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
  ===  File description
  ===========================================================================*/

/*!
  \file   spectral.h
  \author Jakob Doerr <jakobdoerr@googlemail.com>
  \date   2018-05-18

  \brief  Spectral transformation function definitions.

   This file contains the definition of the  functions in spectral.cc
   that are of interest elsewhere.
*/

#ifndef spectral_h
#define spectral_h

#include "matpackVII.h"
#include "mystring.h"
#include "messages.h"
#include "gridded_fields.h"
#include "propagationmatrix.h"


void opt_prop_SpecToGrid(//Output
        Tensor5View ext_matrix,
        Tensor4View abs_vector,
        //Input
        const Tensor5& ext_matrix_spectral,
        const Tensor4& abs_vector_spectral,
        const Matrix& dir_array,
        const bool& any_m);

void pha_mat_SpecToGrid(//Output
        Tensor6View pha_matrix,
        //Input
        const Tensor6& pha_mat_real_spectral,
        const Tensor6& pha_mat_imag_spectral,
        const Matrix& pdir_array,
        const Matrix& idir_array,
        const Index& ptype,
        const bool& any_m_inc,
        const bool& any_m_sca);

#endif /* spectral_h */
