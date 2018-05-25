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

#include "matpackVII.h"
#include "mystring.h"
#include "messages.h"
#include "gridded_fields.h"
#include "propagationmatrix.h"

void spec_to_grid(//Output
        MatrixView grid_field,
        //Input
        const Vector& spec_field,
        const Vector& za_grid,
        const Vector& aa_grid);

void par_optpropSpecToGrid(//Output
        Tensor6View ext_matrix,
        Tensor5View abs_vector,
        //Input
        const Tensor5& ext_matrix_spectral,
        const Tensor4& abs_vector_spectral,
        const Vector& za_grid,
        const Vector& aa_grid);
