/* Copyright (C) 2017 Oliver Lemke <olemke@core-dump.info> and Stefan
   Buehler <sbuehler@ltu.se>. 

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
  \file   m_hitran_xsec.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \author Stefan Buehler
  \date   2017-09-28

  \brief  Workspace methods for HITRAN absorption cross section data.

*/

#include "arts.h"
#include "absorption.h"
#include "messages.h"
#include "physics_funcs.h"
#include "auto_md.h"
#include "xml_io.h"


extern const Numeric SPEED_OF_LIGHT;

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddHitranXsec(// WS Output:
        ArrayOfMatrix& abs_xsec_per_species,
        ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
        // WS Input:
        const ArrayOfArrayOfSpeciesTag& abs_species,
        const ArrayOfRetrievalQuantity& jacobian_quantities,
        const ArrayOfIndex& abs_species_active,
        const Vector& f_grid,
        const Vector& abs_p,
        const Vector& abs_t,
        const Matrix& hitran_xsec_data,
        // WS Generic Input:
        const Numeric& T_extrapolfac,
        const Index& robust,
        // Verbosity object:
        const Verbosity& verbosity)
{
    CREATE_OUTS;

    {
        // Check that all parameters that should have the number of tag
        // groups as a dimension are consistent:
        const Index n_tgs = abs_species.nelem();
        const Index n_xsec = abs_xsec_per_species.nelem();

        if (n_tgs != n_xsec)
        {
            ostringstream os;
            os << "The following variables must all have the same dimension:\n"
               << "abs_species:          " << n_tgs << "\n"
               << "abs_xsec_per_species: " << n_xsec;
            throw runtime_error(os.str());
        }
    }

    // Jacobian overhead START
    /* NOTE:  The calculations below are inefficient and could
              be made much better by using interp in Extract to
              return the derivatives as well. */
    const PropmatPartialsData ppd(jacobian_quantities);
    const bool do_jac = ppd.supportsHitranXsec();
    const bool do_freq_jac = ppd.do_frequency();
    const bool do_temp_jac = ppd.do_temperature();
    Vector dfreq, dabs_t;
    const Numeric df = ppd.Frequency_Perturbation();
    const Numeric dt = ppd.Temperature_Perturbation();
    if (do_freq_jac)
    {
        dfreq.resize(f_grid.nelem());
        dfreq = f_grid;
        dfreq += df;
    }
    if (do_temp_jac)
    {
        dabs_t.resize(abs_t.nelem());
        dabs_t = abs_t;
        dabs_t += dt;
    }
    // Jacobian overhead END

    // Useful if there is no Jacobian to calculate
    ArrayOfMatrix empty;

    {
        // Check that all parameters that should have the the dimension of p_grid
        // are consistent:
        const Index n_p = abs_p.nelem();
        const Index n_t = abs_t.nelem();

        if (n_p != n_t)
        {
            ostringstream os;
            os << "The following variables must all have the same dimension:\n"
               << "abs_p:          " << n_p << "\n"
               << "abs_t:          " << n_t;
            throw runtime_error(os.str());
        }
    }

    // Allocate a vector with dimension frequencies for constructing our
    // cross-sections before adding them (more efficient to allocate this here
    // outside of the loops)
    Vector xsec_temp(f_grid.nelem(), 0.);

    // Jacobian vectors START
    Vector dxsec_temp_dT;
    Vector dxsec_temp_dF;
    if (do_freq_jac)
        dxsec_temp_dF.resize(f_grid.nelem());
    if (do_temp_jac)
        dxsec_temp_dT.resize(f_grid.nelem());
    // Jacobian vectors END

    // Loop over CIA data sets.
    // Index ii loops through the outer array (different tag groups),
    // index s through the inner array (different tags within each goup).
    for (Index ii = 0; ii < abs_species_active.nelem(); ii++)
    {
        const Index i = abs_species_active[ii];

        for (Index s = 0; s < abs_species[i].nelem(); s++)
        {
            const SpeciesTag& this_species = abs_species[i][s];

            // Check if this is a CIA tag
            if (this_species.Type() != SpeciesTag::TYPE_HITRAN_XSEC)
                continue;

//            const CIARecord& this_cia = abs_cia_data[this_cia_index];
            Matrix& this_xsec = abs_xsec_per_species[i];
            ArrayOfMatrix& this_dxsec = do_jac ? dabs_xsec_per_species_dx[i]
                                               : empty;

            // Check that the dimension of this_xsec is
            // consistent with abs_p and f_grid.
            if (this_xsec.nrows() != f_grid.nelem())
            {
                ostringstream os;
                os
                        << "Wrong dimension of abs_xsec_per_species.nrows for species "
                        << i << ":\n"
                        << "should match f_grid (" << f_grid.nelem()
                        << ") but is "
                        << this_xsec.nrows() << ".";
                throw runtime_error(os.str());
            }
            if (this_xsec.ncols() != abs_p.nelem())
            {
                ostringstream os;
                os
                        << "Wrong dimension of abs_xsec_per_species.ncols for species "
                        << i << ":\n"
                        << "should match abs_p (" << abs_p.nelem()
                        << ") but is "
                        << this_xsec.ncols() << ".";
                throw runtime_error(os.str());
            }

            // Loop over pressure:
            for (Index ip = 0; ip < abs_p.nelem(); ip++)
            {
                // Get the binary absorption cross sections from the HITRAN data:
                try
                {
                    // TODO OLE: CALCULATE CROSS SECTION VALUE HERE!!!!!!!!!
                    xsec_temp = 1e-17;
//                    this_cia.Extract(xsec_temp, f_grid, abs_t[ip],
//                                     this_species.CIADataset(),
//                                     T_extrapolfac, robust, verbosity);
                    //                    if(do_freq_jac)
                    //                        this_cia.Extract(dxsec_temp_dF, dfreq, abs_t[ip],
                    //                                         this_species.CIADataset(),
                    //                                         T_extrapolfac, robust, verbosity);
                    //                    if(do_temp_jac)
                    //                        this_cia.Extract(dxsec_temp_dT, f_grid, dabs_t[ip],
                    //                                         this_species.CIADataset(),
                    //                                         T_extrapolfac, robust, verbosity);
                } catch (runtime_error e)
                {
                    ostringstream os;
                    os << "Problem with HITRAN cross section species " << this_species.Name()
                       << ":\n"
                       << e.what();
                    throw runtime_error(os.str());
                }


                if (!do_jac)
                {
                    // Add to result variable:
                    this_xsec(joker, ip) += xsec_temp;
                }
                else
                {
                    for (Index iv = 0; iv < xsec_temp.nelem(); iv++)
                    {
                        this_xsec(iv, ip) += xsec_temp[iv];
                        for (Index iq = 0; iq < ppd.nelem(); iq++)
                        {
                            if (ppd(iq) == JQT_frequency ||
                                ppd(iq) == JQT_wind_magnitude ||
                                ppd(iq) == JQT_wind_u ||
                                ppd(iq) == JQT_wind_v || ppd(iq) == JQT_wind_w)
                                this_dxsec[iq](iv, ip) += (dxsec_temp_dF[iv] -
                                                           xsec_temp[iv]) / df;
                            else if (ppd(iq) == JQT_temperature)
                                this_dxsec[iq](iv, ip) += (dxsec_temp_dT[iv] -
                                                           xsec_temp[iv]) / dt;
                            else if (ppd(iq) == JQT_VMR)
                            {
                                if (abs_species[i][0].Species() ==
                                    ppd.species(iq))
                                {
                                    this_dxsec[iq](iv, ip) += xsec_temp[iv];
                                }
                            }
                            // Note for coef that d/dt(a*n*n) = da/dt * n1*n2 + a * dn1/dt * n2 + a * n1 * dn2/dt, 
                            // we now output da/dt*n2 + a*dn2/dt and coef conversion then have to do 
                            // dxsec*n1 + xsec*dn1/dt, which is what it has to do anyways, so no problems!
                            // Also note that d/dvmr gives a factor for other species absorption.
                            // even if own absorption is zero...
                        }
                    }
                }
            }
        }
    }
}

