/* Copyright (C) 2018 Oliver Lemke <olemke@core-dump.info>

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
  \file   hitran_xsec.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2018-01-08

  \brief  Methods and classes for HITRAN absorption cross section data.
*/


#include "hitran_xsec.h"
#include "absorption.h"
#include "check_input.h"

extern const Numeric PI;

Numeric func_2straights(const Numeric x, const Vector& coeffs)
{
    return (x <= coeffs[0]) ? coeffs[1] * x
                            : coeffs[2] * (x - coeffs[0]) +
                              coeffs[1] * coeffs[0];
}

Numeric lorentz_pdf(const Numeric x, const Numeric x0, const Numeric gamma)
{
    const Numeric xx0 = x - x0;
    return gamma / PI / (xx0 * xx0 + gamma * gamma);

}

String XsecRecord::SpeciesName() const
{
    // The function species_name_from_species_index internally does an assertion
    // that the species with this index really exists.
    return species_name_from_species_index(mspecies);
}


void XsecRecord::Extract(VectorView result,
                         ConstVectorView f_grid,
                         const Numeric& pressure,
                         const Verbosity& verbosity) const
{
    CREATE_OUTS;

    const Index nf = f_grid.nelem();

    // Assert that result vector has right size:
    assert(result.nelem() == nf);

    // Initialize result to zero (important for those frequencies outside the data grid).
    result = 0.;

    const Index ndatasets = mxsecs.nelem();
    for (Index idataset = 0; idataset < ndatasets; idataset++)
    {
        const Vector& data_f_grid = mfgrids[idataset];
        const Numeric fmin = data_f_grid[0];
        const Numeric fmax = data_f_grid[mfgrids[idataset].nelem() - 1];

        if (out3.sufficient_priority())
        {
            // Some detailed information to the most verbose output stream:
            ostringstream os;
            os << "    f_grid:      " << f_grid[0] << " - "
               << f_grid[nf - 1] << " Hz\n"
               << "    data_f_grid: " << fmin << " - " << fmax << " Hz\n"
               << "    pressure: " << pressure << " K\n";
            out3 << os.str();
        }

        // We want to return result zero for all f_grid points that are outside the
        // data_f_grid, because xsec datasets are defined only where the absorption
        // was measured. So, we have to find out which part of f_grid is inside
        // data_f_grid.
        Index i_fstart, i_fstop;

        for (i_fstart = 0; i_fstart < nf; ++i_fstart)
            if (f_grid[i_fstart] >= fmin) break;

        // Return directly if all frequencies are below data_f_grid:
        if (i_fstart == nf) return;

        for (i_fstop = nf - 1; i_fstop >= 0; --i_fstop)
            if (f_grid[i_fstop] <= fmax) break;

        // Return directly if all frequencies are above data_f_grid:
        if (i_fstop == -1) return;

        // Extent for active frequency vector:
        const Index f_extent = i_fstop - i_fstart + 1;

        if (out3.sufficient_priority())
        {
            ostringstream os;
            os << "    " << f_extent
               << " frequency extraction points starting at "
               << "frequency index " << i_fstart << ".\n";
            out3 << os.str();
        }

        // If f_extent is less than one, then the entire data_f_grid is between two
        // grid points of f_grid. (So that we do not have any f_grid points inside
        // data_f_grid.) Return also in this case.
        if (f_extent < 1) return;


        // This is the part of f_grid for which we have to do the interpolation.
        ConstVectorView f_grid_active = f_grid[Range(i_fstart, f_extent)];

        // We have to create a matching view on the result vector:
        VectorView result_active = result[Range(i_fstart, f_extent)];
        Vector xsec_interp(f_extent);


        // Decide on interpolation orders:
        const Index f_order = 3;

        // The frequency grid has to have enough points for this interpolation
        // order, otherwise throw a runtime error.
        if (data_f_grid.nelem() < f_order + 1)
        {
            ostringstream os;
            os << "Not enough frequency grid points in Hitran Xsec data.\n"
               << "You have only " << data_f_grid.nelem() << " grid points.\n"
               << "But need at least " << f_order + 1 << ".";
            throw runtime_error(os.str());
        }


        // Check if frequency is inside the range covered by the data:
        chk_interpolation_grids("Frequency interpolation for cross sections",
                                data_f_grid,
                                f_grid_active,
                                f_order);


        // Find frequency grid positions:
        ArrayOfGridPosPoly f_gp(f_grid_active.nelem()), T_gp(1);
        gridpos_poly(f_gp, data_f_grid, f_grid_active, f_order);

        {
            Matrix itw(f_gp.nelem(), f_order + 1);
            interpweights(itw, f_gp);
            interp(xsec_interp, itw, mxsecs[idataset], f_gp);
        }

        // Apply pressure dependent broadening and set negative values to zero.
        // (These could happen due to overshooting of the higher order interpolation.)
        const Numeric pdiff = pressure - mrefpressure[idataset];
        const Numeric fwhm = func_2straights(pdiff, mcoeffs);
//        std::cout << mcoeffs << " - ";
//        std::cout << "pdiff: " << pdiff << " - fwhm: " << fwhm << " - ";
//        std::cout << "lorentz: " << lorentz_pdf(24e12, 24e12, fwhm / 2)
//                  << std::endl;
        for (Index i = 0; i < result_active.nelem(); ++i)
        {
            const Numeric x0 = f_grid[i_fstart + i];
            Numeric lsum = 0.;
            Numeric pdf;
            for (Index j = 0; j < result_active.nelem(); j++)
            {
                pdf = lorentz_pdf(f_grid[i_fstart + j], x0, fwhm / 2);
                lsum += pdf;
                result_active[i] += xsec_interp[j] * pdf;
            }
            // Normalize integral to 1
            result_active[i] /= lsum;

            if (result_active[i] < 0)
                result_active[i] = 0;
        }
    }
}


/** Get the index in hitran_xsec_data for the given species.

 \param[in] hitran_xsec_data Hitran Xsec data array
 \param[in] species Species name

 \returns Index of this species in hitran_xsec_data. -1 if not found.
 */
Index hitran_xsec_get_index(const ArrayOfXsecRecord& xsec_data,
                            const Index species)
{
    for (Index i = 0; i < xsec_data.nelem(); i++)
        if (xsec_data[i].Species() == species)
            return i;

    return -1;
}

std::ostream& operator<<(std::ostream& os, const XsecRecord& xd)
{
    os << "Species: " << xd.Species() << std::endl;
    return os;
}
