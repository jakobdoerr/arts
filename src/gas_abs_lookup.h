/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>

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

/*!
  \file   gas_abs_lookup.h
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Sep 19 16:49:00 2002
  
  \brief  Declarations for the gas absorption lookup table.  
*/

#ifndef gas_abs_lookup_h
#define gas_abs_lookup_h

#include "matpackIV.h"
#include "absorption.h"

//! An absorption lookup table.
/*! This class holds an absorption lookup table, as well as all
    information that is necessary to use the table to extract
    absorption. Extraction routines are implemented as member functions. */
struct GasAbsLookup {
public:
  // Documentation is with the implementation!
  void Adapt( const ArrayOfArrayOfSpeciesTag& current_species,
              ConstVectorView current_f_grid );

  // Documentation is with the implementation!
  void Extract( MatrixView sga,
                Index      f_index,
                Numeric    p,
                Numeric    T,
                VectorView vmrs );

  // Documentation is with the implementation!
  void Extract( Tensor3View     sga,
                ConstVectorView p,
                ConstVectorView T);
//private:

  //! The species tags for which the table is valid.
  ArrayOfArrayOfSpeciesTag species; 

  //! The species tags with non-linear treatment.
  /*! This must be inside the range of species. If nonlinear_species
    is an empty vector, it means that all species should be treated
    linearly. (No absorption for perturbed species profiles is stored.) */
  ArrayOfIndex nonlinear_species; 

  //! The frequency grid [Hz].
  /*! Must be sorted in ascending order. */
  Vector    f_grid;

  //! The pressure grid for the table [Pa].
  /*! Must be sorted in decreasing order. */
  Vector    p_grid;  

  //! The base 10 logarithm of the pressure grid.
  /*! This is not stored along with the table, but calculated when the
    table is initialized with Adapt. The reason to have this is that
    vertical interpolation should be linear in log(p). */
  Vector   log_p_grid;  

  //! The reference VMR profiles.
  /*! The VMRs for all species, associated with p_grid. Dimension:
    [n_species, n_p_grid]. These VMRs are needed to scale the
    absorption coefficient to other VMRs. We are never working with
    "absorption cross-sections", always with real absorption coefficients,
    so we have to remember the associated VMR values. 

    Physical unit: Absolute value. */
  Matrix    vmrs_ref;

  //! The reference temperature profile [K].
  /*! This is a temperature profile. The dimension must be the same as
    p_grid. */
  Vector    t_ref;

  //! The vector of temperature perturbations [K].
  /*! This can have any number of elements. Example:
    [-20,-10,0,10,20]. The actual temperatures for which absorption is
    stored are t_ref + t_pert for each level. The reference
    temperature itself should normally also be included, hence t_pert should
    always include 0. Must be sorted in ascending order!

    The vector t_pert may be an empty vector (nelem()=0), which marks
    the special case that no interpolation in temperature should be
    done. If t_pert is not empty, you will get an error message if you
    try to extract absorption for temperatures outside the range of
    t_pert. */
  Vector    t_pert;

  //! The vector of perturbations for the VMRs of the nonlinear species.
  /*!
    These apply to all the species that have been set as
    nonlinear_species.

    Fractional units are used! Example: [0,.5,1,10,100],
    meaning from VMR 0 to 100 times the profile given in
    vmrs. The reference value should normally be included, hence
    nls_pert should always include the value 1.

    If nonlinear_species is an empty vector, it means that there are
    no nonlinear species. Then nls_pert must also be an empty vector.
  */
  Vector    nls_pert;

  //! Absorption coefficients.
  /*!
    Physical unit: 1/m

    Dimension: [ a, b, c, d ]

    Simplest case (no temperature perturbations, no vmr perturbations): <br>
    a = 1 <br>
    b = n_species <br>
    c = n_f_grid <br>
    d = n_p_grid <br>

    Standard case (temperature perturbations, but no vmr perturbations): <br>
    a = n_t_pert <br>
    b = n_species <br>
    c = n_f_grid <br>
    d = n_p_grid <br>

    Full case (with temperature perturbations and vmr perturbations): <br>
    a = n_t_pert <br>
    b = n_species + n_nonlinear_species * ( n_nls_pert - 1 ) <br>
    c = n_f_grid <br>
    d = n_p_grid <br>

    Note that the last three dimensions are identical to the
    dimensions of abs_per_tg in ARTS-1-0. This should simplify
    computation of the lookup table with this old ARTS version.
  */
  Tensor4 abs;

};

#endif //  gas_abs_lookup_h
