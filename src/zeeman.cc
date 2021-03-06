/* Copyright (C) 2014
   Richard Larsson <ric.larsson@gmail.com>

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


#include "zeeman.h"
#include "global_data.h"
#include "linefunctions.h"
#include "linescaling.h"

/*!
 *  Defines the phase of the propagation matrix.
 *  Read Larsson et al. (2013) for an explanation.
 * 
 *  \param  K       Out:    The rotation extinction matrix.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2012-08-03
 */
void phase_matrix(MatrixView K, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(K.nrows() == 4 );
    assert(K.ncols() == 4 );
    
    const Numeric
    S2T = sin(theta)*sin(theta),  
    CT  = cos(theta), 
    CE2 = cos(2*eta), 
    SE2 = sin(2*eta);
    
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            K(0,0) =   0;  K(0,1) =           0;  K(0,2) =           0;  K(0,3) =           0;
            K(1,0) =   0;  K(1,1) =           0;  K(1,2) =      2 * CT;  K(1,3) =   S2T * SE2;
            K(2,0) =   0;  K(2,1) =    - 2 * CT;  K(2,2) =           0;  K(2,3) = - S2T * CE2;
            K(3,0) =   0;  K(3,1) = - S2T * SE2;  K(3,2) =   S2T * CE2;  K(3,3) =           0;
            break;
        case  1: // Transitions parallel to the magnetic field
            K(0,0) =   0;  K(0,1) =           0;  K(0,2) =           0;  K(0,3) =           0;
            K(1,0) =   0;  K(1,1) =           0;  K(1,2) =    - 2 * CT;  K(1,3) =   S2T * SE2;
            K(2,0) =   0;  K(2,1) =      2 * CT;  K(2,2) =           0;  K(2,3) = - S2T * CE2;
            K(3,0) =   0;  K(3,1) = - S2T * SE2;  K(3,2) =   S2T * CE2;  K(3,3) =           0;
            break;
        case  0:// Transitions perpendicular to the magnetic field
            K(0,0) =   0;  K(0,1) =           0;  K(0,2) =           0;  K(0,3) =           0;
            K(1,0) =   0;  K(1,1) =           0;  K(1,2) =           0;  K(1,3) = - S2T * SE2;
            K(2,0) =   0;  K(2,1) =           0;  K(2,2) =           0;  K(2,3) =   S2T * CE2;
            K(3,0) =   0;  K(3,1) =   S2T * SE2;  K(3,2) = - S2T * CE2;  K(3,3) =           0;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to phase matrix");
            break;
    };
};


/*!
 *  Defines the attenuation of the propagation matrix.
 *  Read Larsson et al. (2013) for an explanation.
 * 
 *  \param  K       Out:    The rotation extinction matrix.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2012-08-03
 */
void attenuation_matrix(MatrixView K, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(K.nrows() == 4 );
    assert(K.ncols() == 4 );
    
    const Numeric 
    S2T  = sin(theta)*sin(theta),  
    C2T  = cos(theta)*cos(theta),  
    CT   = cos(theta), 
    CE2  = cos(2*eta), 
    SE2  = sin(2*eta);
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            K(0,0) =     1 + C2T;  K(0,1) =   S2T * CE2;  K(0,2) =   S2T * SE2;  K(0,3) =    2 * CT;
            K(1,0) =   S2T * CE2;  K(1,1) =     1 + C2T;  K(1,2) =           0;  K(1,3) =         0;
            K(2,0) =   S2T * SE2;  K(2,1) =           0;  K(2,2) =     1 + C2T;  K(2,3) =         0;
            K(3,0) =      2 * CT;  K(3,1) =           0;  K(3,2) =           0;  K(3,3) =   1 + C2T;
            break;
        case  1: // Transitions parallel to the magnetic field
            K(0,0) =     1 + C2T;  K(0,1) =   S2T * CE2;  K(0,2) =   S2T * SE2;  K(0,3) =  - 2 * CT;
            K(1,0) =   S2T * CE2;  K(1,1) =     1 + C2T;  K(1,2) =           0;  K(1,3) =         0;
            K(2,0) =   S2T * SE2;  K(2,1) =           0;  K(2,2) =     1 + C2T;  K(2,3) =         0;
            K(3,0) =    - 2 * CT;  K(3,1) =           0;  K(3,2) =           0;  K(3,3) =   1 + C2T;
            break;
        case  0: // Transitions perpendicular to the magnetic field
            K(0,0) =         S2T;  K(0,1) = - S2T * CE2;  K(0,2) = - S2T * SE2;  K(0,3) =         0;
            K(1,0) = - S2T * CE2;  K(1,1) =         S2T;  K(1,2) =           0;  K(1,3) =         0;
            K(2,0) = - S2T * SE2;  K(2,1) =           0;  K(2,2) =         S2T;  K(2,3) =         0;
            K(3,0) =           0;  K(3,1) =           0;  K(3,2) =           0;  K(3,3) =       S2T;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to attenuation matrix");
            break;
    };
};


/*!
 *  Defines the derivative of the phase of the propagation matrix with regards to theta
 * 
 *  \param  dK      Out:    The rotation extinction matrix derivative with regards to theta.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2015-08-14
 */
void dphase_matrix_dtheta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(dK.nrows() == 4 );
    assert(dK.ncols() == 4 );
    
    const Numeric 
    ST   = sin(theta),
    CTST = sin(theta)*cos(theta), 
    CE2  = cos(2*eta), 
    SE2  = sin(2*eta);
    
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            dK(0,0) =   0;  dK(0,1) =                0;  dK(0,2) =                0;  dK(0,3) =                0;
            dK(1,0) =   0;  dK(1,1) =                0;  dK(1,2) =         - 2 * ST;  dK(1,3) =   2 * CTST * SE2;
            dK(2,0) =   0;  dK(2,1) =           2 * ST;  dK(2,2) =                0;  dK(2,3) = - 2 * CTST * CE2;
            dK(3,0) =   0;  dK(3,1) = - 2 * CTST * SE2;  dK(3,2) =   2 * CTST * CE2;  dK(3,3) =                0;
            break;
        case  1: // Transitions parallel to the magnetic field
            dK(0,0) =   0;  dK(0,1) =                0;  dK(0,2) =                0;  dK(0,3) =                0;
            dK(1,0) =   0;  dK(1,1) =                0;  dK(1,2) =           2 * ST;  dK(1,3) =   2 * CTST * SE2;
            dK(2,0) =   0;  dK(2,1) =        -  2 * ST;  dK(2,2) =                0;  dK(2,3) = - 2 * CTST * CE2;
            dK(3,0) =   0;  dK(3,1) = - 2 * CTST * SE2;  dK(3,2) =   2 * CTST * CE2;  dK(3,3) =                0;
            break;
        case  0:// Transitions perpendicular to the magnetic field
            dK(0,0) =   0;  dK(0,1) =                0;  dK(0,2) =                0;  dK(0,3) =                0;
            dK(1,0) =   0;  dK(1,1) =                0;  dK(1,2) =                0;  dK(1,3) = - 2 * CTST * SE2;
            dK(2,0) =   0;  dK(2,1) =                0;  dK(2,2) =                0;  dK(2,3) =   2 * CTST * CE2;
            dK(3,0) =   0;  dK(3,1) =   2 * CTST * SE2;  dK(3,2) = - 2 * CTST * CE2;  dK(3,3) =                0;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to phase matrix");
            break;
    };
};


/*!
 *  Defines the derivative of the phase of the propagation matrix with regards to eta
 * 
 *  \param  dK      Out:    The rotation extinction matrix derivative with regards to eta.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2015-08-14
 */
void dphase_matrix_deta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(dK.nrows() == 4 );
    assert(dK.ncols() == 4 );
    
    const Numeric 
    S2T = sin(theta)*sin(theta),  
    CE2 = cos(2*eta),
    SE2 = sin(2*eta);
    
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            dK(0,0) =   0;  dK(0,1) =               0;  dK(0,2) =               0;  dK(0,3) =               0;
            dK(1,0) =   0;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   2 * S2T * CE2;
            dK(2,0) =   0;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   2 * S2T * SE2;
            dK(3,0) =   0;  dK(3,1) = - 2 * S2T * CE2;  dK(3,2) = - 2 * S2T * SE2;  dK(3,3) =               0;
            break;
        case  1: // Transitions parallel to the magnetic field
            dK(0,0) =   0;  dK(0,1) =               0;  dK(0,2) =               0;  dK(0,3) =               0;
            dK(1,0) =   0;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   2 * S2T * CE2;
            dK(2,0) =   0;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   2 * S2T * SE2;
            dK(3,0) =   0;  dK(3,1) = - 2 * S2T * CE2;  dK(3,2) = - 2 * S2T * SE2;  dK(3,3) =               0;
            break;
        case  0:// Transitions perpendicular to the magnetic field
            dK(0,0) =   0;  dK(0,1) =               0;  dK(0,2) =               0;  dK(0,3) =               0;
            dK(1,0) =   0;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) = - 2 * S2T * CE2;
            dK(2,0) =   0;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) = - 2 * S2T * SE2;
            dK(3,0) =   0;  dK(3,1) =   2 * S2T * CE2;  dK(3,2) =   2 * S2T * SE2;  dK(3,3) =               0;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to phase matrix");
            break;
    };
};


/*!
 *  Defines the derivative of the attenuation of the propagation matrix with regards to theta.
 * 
 *  \param  dK      Out:    The rotation extinction matrix derivative with regards to theta.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2015-08-14
 */
void dattenuation_matrix_dtheta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(dK.nrows() == 4 );
    assert(dK.ncols() == 4 );
    
    const Numeric 
    CTST = cos(theta)*sin(theta), 
    ST   = sin(theta),
    CE2  = cos(2*eta), 
    SE2  = sin(2*eta);
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            dK(0,0) =       - 2 * CTST;  dK(0,1) =   2 * CTST * CE2;  dK(0,2) =   2 * CTST * SE2;  dK(0,3) =   - 2 * ST;
            dK(1,0) =   2 * CTST * CE2;  dK(1,1) =       - 2 * CTST;  dK(1,2) =                0;  dK(1,3) =          0;
            dK(2,0) =   2 * CTST * SE2;  dK(2,1) =                0;  dK(2,2) =       - 2 * CTST;  dK(2,3) =          0;
            dK(3,0) =         - 2 * ST;  dK(3,1) =                0;  dK(3,2) =                0;  dK(3,3) = - 2 * CTST;
            break;
        case  1: // Transitions parallel to the magnetic field
            dK(0,0) =       - 2 * CTST;  dK(0,1) =   2 * CTST * CE2;  dK(0,2) =   2 * CTST * SE2;  dK(0,3) =     2 * ST;
            dK(1,0) =   2 * CTST * CE2;  dK(1,1) =       - 2 * CTST;  dK(1,2) =                0;  dK(1,3) =          0;
            dK(2,0) =   2 * CTST * SE2;  dK(2,1) =                0;  dK(2,2) =       - 2 * CTST;  dK(2,3) =          0;
            dK(3,0) =           2 * ST;  dK(3,1) =                0;  dK(3,2) =                0;  dK(3,3) = - 2 * CTST;
            break;
        case  0: // Transitions perpendicular to the magnetic field
            dK(0,0) =         2 * CTST;  dK(0,1) = - 2 * CTST * CE2;  dK(0,2) = - 2 * CTST * SE2;  dK(0,3) =          0;
            dK(1,0) = - 2 * CTST * CE2;  dK(1,1) =         2 * CTST;  dK(1,2) =                0;  dK(1,3) =          0;
            dK(2,0) = - 2 * CTST * SE2;  dK(2,1) =                0;  dK(2,2) =         2 * CTST;  dK(2,3) =          0;
            dK(3,0) =                0;  dK(3,1) =                0;  dK(3,2) =                0;  dK(3,3) =   2 * CTST;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to attenuation matrix");
            break;
    };
};


/*!
 *  Defines the derivative of the attenuation of the propagation matrix with regards to eta.
 * 
 *  \param  dK      Out:    The rotation extinction matrix derivative with regards to eta.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2015-08-14
 */
void dattenuation_matrix_deta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(dK.nrows() == 4 );
    assert(dK.ncols() == 4 );
    
    const Numeric
    SE2 = sin(2*eta),
    CE2 = cos(2*eta),
    ST2 = sin(theta)*sin(theta);
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            dK(0,0) =               0;  dK(0,1) = - 2 * ST2 * SE2;  dK(0,2) =   2 * ST2 * CE2;  dK(0,3) =   0;
            dK(1,0) = - 2 * ST2 * SE2;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   0;
            dK(2,0) =   2 * ST2 * CE2;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   0;
            dK(3,0) =               0;  dK(3,1) =               0;  dK(3,2) =               0;  dK(3,3) =   0;
            break;
        case  1: // Transitions parallel to the magnetic field
            dK(0,0) =               0;  dK(0,1) = - 2 * ST2 * SE2;  dK(0,2) =   2 * ST2 * CE2;  dK(0,3) =   0;
            dK(1,0) = - 2 * ST2 * SE2;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   0;
            dK(2,0) =   2 * ST2 * CE2;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   0;
            dK(3,0) =               0;  dK(3,1) =               0;  dK(3,2) =               0;  dK(3,3) =   0;
            break;
        case  0: // Transitions perpendicular to the magnetic field
            dK(0,0) =               0;  dK(0,1) =   2 * ST2 * SE2;  dK(0,2) = - 2 * ST2 * CE2;  dK(0,3) =   0;
            dK(1,0) =   2 * ST2 * SE2;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   0;
            dK(2,0) = - 2 * ST2 * CE2;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   0;
            dK(3,0) =               0;  dK(3,1) =               0;  dK(3,2) =               0;  dK(3,3) =   0;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to attenuation matrix");
            break;
    };
};


Numeric gs_caseb(const Rational& N, const Rational& J, const Rational& S, const Numeric& GS) { return (GS*((J*(J+1)+S*(S+1)-N*(N+1))/(J*(J+1))/2).toNumeric()); }
Numeric gs_casea(const Rational& Omega, const Rational& J, const Rational& Lambda, const Rational& Sigma, const Numeric& GS) { return GS*(Omega/2/J/(J+1)*(Lambda+2*Sigma)).toNumeric(); }


/*!
 * Return the relative strength of the split Zeeman line parts as found in
 * Berdyugina and Solnaki (2002). Note that this is the same as the general case
 * of Schadee (1978).
 * 
 * \param  __U__   Void.
 * \param  m       In:     Secondary rotational quantum number.
 * \param  j       In:     Spin-Orbit Coupling number.
 * \param  DJ      In:     Change in the main rotational quantum number.
 * \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 * \author Richard Larsson
 * \date   2012-10-26
 */
Numeric relative_strength(const Rational& M, const Rational& J, const Index& dj, const Index& dm)
{
    
    // Variable to be returned.
    Rational ret_val;
    
    static const Rational three_fourths(3,4),one_and_a_half(3,2);
    
    switch ( dj )
    {
        case -1:
            switch ( dm )
            {
                case -1: // Transitions anti-parallel to the magnetic field
                    ret_val = three_fourths*(J+M)*(J-1+M)/(2*J*(2*J-1)*(2*J+1));
                    break;
                case  0: // Transitions perpendicular to the magnetic field
                    ret_val = one_and_a_half*(J*J-M*M)/(J*(2*J-1)*(2*J+1));
                    break;
                case +1: // Transitions parallel to the magnetic field
                    ret_val = three_fourths*(J-M)*(J-1-M)/(2*J*(2*J-1)*(2*J+1));
                    break;
                default:
                    throw std::runtime_error("Something is extremely wrong.");
                    break;
            }
            break;
        case  0:
            switch ( dm )
            {
                case -1: // Transitions anti-parallel to the magnetic field
                    ret_val = three_fourths*(J+M)*(J+1-M)/(2*J*(J+1)*(2*J+1));
                    break;
                case  0: // Transitions perpendicular to the magnetic field
                    ret_val = one_and_a_half*M*M/(J*(J+1)*(2*J+1));
                    break;
                case +1: // Transitions parallel to the magnetic field
                    ret_val = three_fourths*(J-M)*(J+1+M)/(2*J*(J+1)*(2*J+1));
                    break;
                default:
                    throw std::runtime_error("Something is extremely wrong.");
                    break;
            }
            break;
        case +1:
            switch ( dm )
            {
                case -1: // Transitions anti-parallel to the magnetic field
                    ret_val = three_fourths*(J+1-M)*(J+2-M)/(2*(J+1)*(2*J+1)*(2*J+3));
                    break;
                case  0: // Transitions perpendicular to the magnetic field
                    ret_val = one_and_a_half*((J+1)*(J+1)-M*M)/((J+1)*(2*J+1)*(2*J+3));
                    break;
                case +1: // Transitions parallel to the magnetic field
                    ret_val = three_fourths*(J+1+M)*(J+2+M)/(2*(J+1)*(2*J+1)*(2*J+3));
                    break;
                default:
                    throw std::runtime_error("Something is extremely wrong.");
                    break;
            }
            break;
        default:
            throw std::runtime_error("Something is extremely wrong.");
            break;
    }

    return ret_val.toNumeric();
}


/*!
 * Return the frequency change of the split Zeeman lines following
 * Berdyugina and Solnaki (2002) and Lenoir (1967).
 * 
 * \param  lr      In:     The line record with quantum numbers.
 * \param  H_mag   In:     The magnitude of the magnetic field.
 * \param  GS      In:     The Lande factor.
 * 
 * \author Richard Larsson
 * \date   2016-05-19 (combination of two precious versions from 2012)
 */
Numeric frequency_change(const LineRecord& lr,
                         const Numeric& H_mag, 
                         const Numeric& GS)
{
    Numeric Upper_E_part, Lower_E_part;
    
    // Lower:
    assert(abs(lr.QuantumNumbers().Lower()[QN_M])<=lr.QuantumNumbers().Lower()[QN_J]);
    switch(lr.QuantumNumbers().Lower()[QN_Hund].toIndex())
    {
        case Hund_Case_A:
            // This follows Berdyugina and Solnaki
            Lower_E_part = lr.QuantumNumbers().Lower()[QN_M].toNumeric() * 
            gs_casea(lr.QuantumNumbers().Lower()[QN_Omega],
                     lr.QuantumNumbers().Lower()[QN_J],
                     lr.QuantumNumbers().Lower()[QN_S],
                     lr.QuantumNumbers().Lower()[QN_Lambda],
                     GS);
            break;
        case Hund_Case_B:
            // This follows Lenoir
            if( lr.QuantumNumbers().Lower()[QN_J] == 0 )
            {
                Lower_E_part = 0;
            }
            else
            {
                Lower_E_part = lr.QuantumNumbers().Lower()[QN_M].toNumeric() * 
                gs_caseb(lr.QuantumNumbers().Lower()[QN_N],
                         lr.QuantumNumbers().Lower()[QN_J],
                         lr.QuantumNumbers().Lower()[QN_S],
                         GS);
            }
            break;
        default:
            throw std::runtime_error("Does not recognize Hund case.\n");
            break;
    }
    
    // Upper:
    assert(abs(lr.QuantumNumbers().Upper()[QN_M])<=lr.QuantumNumbers().Upper()[QN_J]);
    switch(lr.QuantumNumbers().Upper()[QN_Hund].toIndex())
    {
        case Hund_Case_A:
            // This follows Berdyugina and Solnaki
            Upper_E_part = lr.QuantumNumbers().Upper()[QN_M].toNumeric() * 
            gs_casea(lr.QuantumNumbers().Upper()[QN_Omega],
                     lr.QuantumNumbers().Upper()[QN_J],
                     lr.QuantumNumbers().Upper()[QN_S],
                     lr.QuantumNumbers().Upper()[QN_Lambda],
                     GS);
            break;
        case Hund_Case_B:
            // This follows Lenoir
            Upper_E_part = lr.QuantumNumbers().Upper()[QN_M].toNumeric() * 
            gs_caseb(lr.QuantumNumbers().Upper()[QN_N],
                     lr.QuantumNumbers().Upper()[QN_J],
                     lr.QuantumNumbers().Upper()[QN_S],
                     GS);
            break;
        default:
            throw std::runtime_error("Does not recognize Hund case.\n");
            break;
    }
    
    // convert from energy state to frequency and be done with it
    return H_mag * (Lower_E_part-Upper_E_part) / PLANCK_CONST * BOHR_MAGNETON;
    
}


void xsec_species_line_mixing_wrapper_with_zeeman(  ArrayOfPropagationMatrix& propmat_clearsky, 
                                                    ArrayOfStokesVector& nlte_source,
                                                    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                                    ArrayOfStokesVector& dnlte_dx_source,
                                                    ArrayOfStokesVector& nlte_dsource_dx,
                                                    const ArrayOfArrayOfSpeciesTag& abs_species, 
                                                    const PropmatPartialsData& flag_partials,
                                                    const Index& abs_lineshape_ls, 
                                                    const Index& abs_lineshape_lsn, 
                                                    const Numeric& abs_lineshape_cutoff, 
                                                    const ArrayOfLineRecord& lr, 
                                                    const Vector&  Zeeman_DF,
                                                    const Vector&  planck_BT,
                                                    const Matrix&  dplanck_BT,
                                                    const SpeciesAuxData& isotopologue_ratios, 
                                                    const SpeciesAuxData& partition_functions,
                                                    const Matrix& abs_t_nlte, 
                                                    const Matrix& abs_vmrs, 
                                                    const Vector& abs_p,
                                                    const Vector& abs_t, 
                                                    const Vector& f_grid,
                                                    const Vector& rtp_mag,
                                                    const Vector& r_path_los,
                                                    const Numeric& lm_p_lim,
                                                    const Numeric& theta, 
                                                    const Numeric& eta, 
                                                    const Numeric& H_mag, 
                                                    const Index& DM, 
                                                    const Index& this_species,
                                                    const Verbosity& verbosity )
                                                
                                                
{
  const bool do_src  =  !nlte_source.empty();
  const Index nq = flag_partials.nelem();
  const Index nf = f_grid.nelem();
  const Numeric n = abs_vmrs(this_species, 0)*number_density( abs_p[0],abs_t[0]);
  const Numeric dn_dT = abs_vmrs(this_species, 0)*dnumber_density_dt( abs_p[0],abs_t[0]);
  
  // Setting up variables
  Matrix attenuation(nf, 1),phase(nf, 1),source(do_src?nf:0,do_src?1:0),
  attenuation_du(nf, 1),phase_du(nf, 1),attenuation_dv(nf, 1),phase_dv(nf, 1),attenuation_dw(nf, 1),phase_dw(nf, 1),
  source_du(do_src?nf:0,do_src?1:0),source_dv(do_src?nf:0,do_src?1:0),source_dw(do_src?nf:0,do_src?1:0);
  
  
  ArrayOfMatrix partial_attenuation(nq), partial_phase(nq), partial_source(do_src?nq:0);
  for(Index iq = 0; iq<nq; iq++)
  {
    if(flag_partials(iq)!=JQT_NOT_JQT)
    {
      partial_attenuation[iq].resize(nf, 1);
      partial_phase[iq].resize(nf, 1);
      if(do_src)
        partial_source[iq].resize(nf, 1);
    }
  }
  // FIXME: Have to perturb for magnetic u, v, and w since it is too complicated otherwise
  Numeric dB=0.0;
  if(nq)
    dB = flag_partials.Magnetic_Field_Perturbation();
  
  Numeric H_dummy, deta_du, dtheta_du, deta_dv, dtheta_dv, deta_dw, dtheta_dw;
  Vector dmag = rtp_mag;
  
  if(flag_partials.do_zeeman_u())
  {
    dmag[0]+=dB;
    set_magnetic_parameters(H_dummy,deta_du,dtheta_du,0,0,0,0,dmag,r_path_los);
    deta_du -= eta;
    deta_du /= dB;
    dtheta_du -= theta;
    dtheta_du /= dB;
    dmag[0]-=dB;
  }
  
  if(flag_partials.do_zeeman_v())
  {
    dmag[1]+=dB;
    set_magnetic_parameters(H_dummy,deta_dv,dtheta_dv,0,0,0,0,dmag,r_path_los);
    deta_dv -= eta;
    deta_dv /= dB;
    dtheta_dv -= theta;
    dtheta_dv /= dB;
    dmag[1]-=dB;
  }
  
  if(flag_partials.do_zeeman_w())
  {
    dmag[2]+=dB;
    set_magnetic_parameters(H_dummy,deta_dw,dtheta_dw,0,0,0,0,dmag,r_path_los);
    deta_dw -= eta;
    deta_dw /= dB;
    dtheta_dw -= theta;
    dtheta_dw /= dB;
    dmag[2]-=dB;
  }
  // JACOBIAN SETUP END
  
  attenuation(joker, 0) = 0.;
  phase(joker, 0) = 0.;
  if(flag_partials.do_zeeman_u())
  {
    attenuation_du(joker,0)=0.;
    phase_du(joker,0)=0.;
  }
  if(flag_partials.do_zeeman_v())
  {
    attenuation_dv(joker,0)=0.;
    phase_dv(joker,0)=0.;
  }
  if(flag_partials.do_zeeman_w())
  {
    attenuation_dw(joker,0)=0.;
    phase_dw(joker,0)=0.;
  }
  if(do_src)
    source(joker, 0)=0.;
  for(Index iq = 0; iq < nq; iq++)
  {
    if(flag_partials(iq)!=JQT_NOT_JQT)
    {
      partial_attenuation[iq](joker,0)=0.;
      partial_phase[iq](joker,0)=0.;
      if(do_src)
        partial_source[iq](joker,0)=0.;
    }
  }
  
  xsec_species_line_mixing_wrapper(   attenuation,         source,         phase, 
                                      partial_attenuation, partial_source, partial_phase, flag_partials,
                                      f_grid, abs_p, abs_t, abs_t_nlte, abs_vmrs, abs_species, 
                                      this_species, lr, Zeeman_DF, H_mag,
                                      abs_lineshape_ls,abs_lineshape_lsn,lm_p_lim,abs_lineshape_cutoff,
                                      isotopologue_ratios, partition_functions, verbosity ); // Now in cross section
  if(flag_partials.do_zeeman_u())
  {
    dmag[0]+=dB;
    
    xsec_species_line_mixing_wrapper(         attenuation_du,         source_du,         phase_du, 
                                              partial_attenuation, partial_source, partial_phase, PropmatPartialsData(ArrayOfRetrievalQuantity(0)),
                                              f_grid, abs_p, abs_t, abs_t_nlte, abs_vmrs, abs_species, 
                                              this_species, lr, Zeeman_DF, sqrt(dmag*dmag),
                                              abs_lineshape_ls,abs_lineshape_lsn,lm_p_lim,abs_lineshape_cutoff,
                                              isotopologue_ratios, partition_functions, verbosity ); // Now in cross section
    dmag[0]-=dB;
  }
  if(flag_partials.do_zeeman_v())
  {
    dmag[1]+=dB;
    
    xsec_species_line_mixing_wrapper(         attenuation_dv,         source_dv,         phase_dv, 
                                              partial_attenuation, partial_source, partial_phase, PropmatPartialsData(ArrayOfRetrievalQuantity(0)),
                                              f_grid, abs_p, abs_t, abs_t_nlte, abs_vmrs, abs_species, 
                                              this_species, lr, Zeeman_DF, sqrt(dmag*dmag),
                                              abs_lineshape_ls,abs_lineshape_lsn,lm_p_lim,abs_lineshape_cutoff,
                                              isotopologue_ratios, partition_functions, verbosity ); // Now in cross section
    dmag[1]-=dB;
  }
  if(flag_partials.do_zeeman_w())
  {
    dmag[2]+=dB;
    
    xsec_species_line_mixing_wrapper(         attenuation_dw,         source_dw,         phase_dw, 
                                              partial_attenuation, partial_source, partial_phase, PropmatPartialsData(ArrayOfRetrievalQuantity(0)),
                                              f_grid, abs_p, abs_t, abs_t_nlte, abs_vmrs, abs_species, 
                                              this_species, lr, Zeeman_DF, sqrt(dmag*dmag),
                                              abs_lineshape_ls,abs_lineshape_lsn,lm_p_lim,abs_lineshape_cutoff,
                                              isotopologue_ratios, partition_functions, verbosity ); // Now in cross section
    dmag[2]-=dB;
  }
  
  if(DM == 0)
  {
    propmat_clearsky[this_species].AddZeemanPiComponent(attenuation(joker, 0), phase(joker, 0), n, theta*DEG2RAD, 
                                                        eta*DEG2RAD);
    
    if(do_src)
      nlte_source[this_species].AddZeemanPiComponent(source(joker,0), phase(joker, 0), n, theta*DEG2RAD, 
                                                     eta*DEG2RAD, planck_BT);
    
    for(Index iq = 0; iq < nq; iq++)
    {
      if(flag_partials(iq)==JQT_magnetic_u)
      {
        attenuation_du -= attenuation;
        attenuation_du /= dB;
        phase_du -= phase;
        phase_du /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentDerivative(attenuation(joker, 0), attenuation_du(joker, 0), 
                                                                phase(joker, 0), phase_du(joker, 0), n, theta*DEG2RAD, 
                                                                dtheta_du*DEG2RAD, eta*DEG2RAD, deta_du*DEG2RAD);
        
        if(do_src)
        {
          source_du -= source;
          source_du /= dB;
          dnlte_dx_source[iq].AddZeemanPiComponentDerivative(source(joker, 0), source_du(joker, 0), phase(joker, 0), 
                                                             phase_du(joker, 0), n, theta*DEG2RAD, dtheta_du*DEG2RAD, 
                                                             eta*DEG2RAD, deta_du*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_v)
      {
        attenuation_dv -= attenuation;
        attenuation_dv /= dB;
        phase_dv -= phase;
        phase_dv /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentDerivative(attenuation(joker, 0), attenuation_dv(joker, 0), 
                                                                phase(joker, 0), phase_dv(joker, 0), n, theta*DEG2RAD, 
                                                                dtheta_dv*DEG2RAD, eta*DEG2RAD, deta_dv*DEG2RAD);
        
        if(do_src)
        {
          source_dv -= source;
          source_dv /= dB;
          dnlte_dx_source[iq].AddZeemanPiComponentDerivative(source(joker, 0), source_dv(joker, 0), phase(joker, 0),
                                                             phase_dv(joker, 0), n, theta*DEG2RAD, dtheta_dv*DEG2RAD, 
                                                             eta*DEG2RAD, deta_dv*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_w)
      {
        attenuation_dw -= attenuation;
        attenuation_dw /= dB;
        phase_dw -= phase;
        phase_dw /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentDerivative(attenuation(joker, 0), attenuation_dw(joker, 0), 
                                                                phase(joker, 0), phase_dw(joker, 0), n, theta*DEG2RAD,
                                                                dtheta_dw*DEG2RAD, eta*DEG2RAD, deta_dw*DEG2RAD);
        
        if(do_src)
        {
          source_dw -= source;
          source_dw /= dB;
          dnlte_dx_source[iq].AddZeemanPiComponentDerivative(source(joker, 0), source_dw(joker, 0), phase(joker, 0),
                                                             phase_dw(joker, 0), n, theta*DEG2RAD, dtheta_dw*DEG2RAD,
                                                             eta*DEG2RAD, deta_dw*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_theta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentThetaDerivative(attenuation(joker, 0), phase(joker, 0), n, 
                                                                     theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanPiComponentThetaDerivative(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                                  eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_eta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentEtaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                   theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dpropmat_clearsky_dx[iq].AddZeemanPiComponentEtaDerivative(source(joker, 0), phase(joker, 0), n, 
                                                                     theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_temperature)
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponent(attenuation(joker, 0), phase(joker, 0), dn_dT, theta*DEG2RAD, 
                                                      eta*DEG2RAD);
        dpropmat_clearsky_dx[iq].AddZeemanPiComponent(partial_attenuation[iq](joker, 0), partial_phase[iq](joker, 0), n, 
                                                      theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanPiComponent(source(joker, 0), phase(joker, 0), dn_dT, theta*DEG2RAD, eta*DEG2RAD, 
                                                   planck_BT);
          dnlte_dx_source[iq].AddZeemanPiComponent(partial_source[iq](joker, 0), phase(joker, 0), n, theta*DEG2RAD, 
                                                   eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanPiComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD, eta*DEG2RAD, 
                                                   dplanck_BT(0, joker));
        }
      }
      else if(flag_partials(iq)==JQT_frequency)
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponent(partial_attenuation[iq](joker, 0), partial_phase[iq](joker, 0), n, 
                                                      theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanPiComponent(partial_source[iq](joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                   eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanPiComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD, eta*DEG2RAD,
                                                   dplanck_BT(1, joker));
        }
      }
      else if(flag_partials(iq)==JQT_VMR)
      {  
        dpropmat_clearsky_dx[iq].AddZeemanPiComponent(attenuation(joker, 0), phase(joker, 0), 
                                                      n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanPiComponent(source(joker, 0), phase(joker, 0), n/abs_vmrs(this_species, 0),
                                                   theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq) not_eq JQT_NOT_JQT)
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponent(partial_attenuation[iq](joker, 0), partial_phase[iq](joker, 0), n,
                                                      theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanPiComponent(partial_source[iq](joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                   eta*DEG2RAD, planck_BT);
        }
      }
    }
  }
  else if(DM == -1)
  {
    propmat_clearsky[this_species].AddZeemanSigmaMinusComponent(attenuation(joker, 0), phase(joker, 0), n,
                                                                theta*DEG2RAD, eta*DEG2RAD);
    
    if(do_src)
      nlte_source[this_species].AddZeemanSigmaMinusComponent(source(joker,0), phase(joker, 0), n, theta*DEG2RAD,
                                                             eta*DEG2RAD, planck_BT);
    
    for(Index iq = 0; iq < nq; iq++)
    {
      if(flag_partials(iq)==JQT_magnetic_u)
      {
        attenuation_du -= attenuation;
        attenuation_du /= dB;
        phase_du -= phase;
        phase_du /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentDerivative(attenuation(joker, 0), attenuation_du(joker, 0),
                                                                        phase(joker, 0), phase_du(joker, 0), n, 
                                                                        theta*DEG2RAD, dtheta_du*DEG2RAD, eta*DEG2RAD,
                                                                        deta_du*DEG2RAD);
        
        if(do_src)
        {
          source_du -= source;
          source_du /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponentDerivative(source(joker, 0), source_du(joker, 0),
                                                                     phase(joker, 0), phase_du(joker, 0), n, 
                                                                     theta*DEG2RAD, dtheta_du*DEG2RAD, eta*DEG2RAD,
                                                                     deta_du*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_v)
      {
        attenuation_dv -= attenuation;
        attenuation_dv /= dB;
        phase_dv -= phase;
        phase_dv /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentDerivative(attenuation(joker, 0), attenuation_dv(joker, 0),
                                                                        phase(joker, 0), phase_dv(joker, 0), n,
                                                                        theta*DEG2RAD, dtheta_dv*DEG2RAD, eta*DEG2RAD,
                                                                        deta_dv*DEG2RAD);
        
        if(do_src)
        {
          source_dv -= source;
          source_dv /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponentDerivative(source(joker, 0), source_dv(joker, 0), 
                                                                     phase(joker, 0), phase_dv(joker, 0), n, 
                                                                     theta*DEG2RAD, dtheta_dv*DEG2RAD, eta*DEG2RAD,
                                                                     deta_dv*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_w)
      {
        attenuation_dw -= attenuation;
        attenuation_dw /= dB;
        phase_dw -= phase;
        phase_dw /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentDerivative(attenuation(joker, 0), attenuation_dw(joker, 0),
                                                                        phase(joker, 0), phase_dw(joker, 0), n,
                                                                        theta*DEG2RAD, dtheta_dw*DEG2RAD, eta*DEG2RAD,
                                                                        deta_dw*DEG2RAD);
        
        if(do_src)
        {
          source_dw -= source;
          source_dw /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponentDerivative(source(joker, 0), source_dw(joker, 0),
                                                                     phase(joker, 0), phase_dw(joker, 0), n,
                                                                     theta*DEG2RAD, dtheta_dw*DEG2RAD, eta*DEG2RAD,
                                                                     deta_dw*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_theta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentThetaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                             theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponentThetaDerivative(source(joker, 0), phase(joker, 0), n,
                                                                          theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_eta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentEtaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                           theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentEtaDerivative(source(joker, 0), phase(joker, 0), n,
                                                                             theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_temperature)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(attenuation(joker, 0), phase(joker, 0), dn_dT,
                                                              theta*DEG2RAD, eta*DEG2RAD);
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(partial_attenuation[iq](joker, 0), 
                                                              partial_phase[iq](joker, 0), n, theta*DEG2RAD, 
                                                              eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(source(joker, 0), phase(joker, 0), dn_dT, theta*DEG2RAD,
                                                           eta*DEG2RAD, planck_BT);
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                           theta*DEG2RAD, eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanSigmaMinusComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                           eta*DEG2RAD, dplanck_BT(0, joker));
        }
      }
      else if(flag_partials(iq)==JQT_frequency)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(partial_attenuation[iq](joker, 0), 
                                                              partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                              eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                           theta*DEG2RAD, eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanSigmaMinusComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                           eta*DEG2RAD, dplanck_BT(1, joker));
        }
      }
      else if(flag_partials(iq)==JQT_VMR)
      {  
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(attenuation(joker, 0), phase(joker, 0),
                                                              n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(source(joker, 0), phase(joker, 0),
                                                           n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD,
                                                           planck_BT);
        }
      }
      else if(flag_partials(iq) not_eq JQT_NOT_JQT)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(partial_attenuation[iq](joker, 0),
                                                              partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                              eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                           theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
    }
  }
  else if(DM == 1)
  {
    propmat_clearsky[this_species].AddZeemanSigmaPlusComponent(attenuation(joker, 0), phase(joker, 0), n,
                                                               theta*DEG2RAD, eta*DEG2RAD);
    
    if(do_src)
      nlte_source[this_species].AddZeemanSigmaPlusComponent(source(joker,0), phase(joker, 0), n, theta*DEG2RAD,
                                                            eta*DEG2RAD, planck_BT);
    
    for(Index iq = 0; iq < nq; iq++)
    {
      if(flag_partials(iq)==JQT_magnetic_u)
      {
        attenuation_du -= attenuation;
        attenuation_du /= dB;
        phase_du -= phase;
        phase_du /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentDerivative(attenuation(joker, 0), attenuation_du(joker, 0),
                                                                       phase(joker, 0), phase_du(joker, 0), n, 
                                                                       theta*DEG2RAD, dtheta_du*DEG2RAD, eta*DEG2RAD, 
                                                                       deta_du*DEG2RAD);
        
        if(do_src)
        {
          source_du -= source;
          source_du /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponentDerivative(source(joker, 0), source_du(joker, 0), 
                                                                    phase(joker, 0), phase_du(joker, 0), n, 
                                                                    theta*DEG2RAD, dtheta_du*DEG2RAD, eta*DEG2RAD, 
                                                                    deta_du*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_v)
      {
        attenuation_dv -= attenuation;
        attenuation_dv /= dB;
        phase_dv -= phase;
        phase_dv /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentDerivative(attenuation(joker, 0), attenuation_dv(joker, 0),
                                                                       phase(joker, 0), phase_dv(joker, 0), n,
                                                                       theta*DEG2RAD, dtheta_dv*DEG2RAD, eta*DEG2RAD,
                                                                       deta_dv*DEG2RAD);
        
        if(do_src)
        {
          source_dv -= source;
          source_dv /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponentDerivative(source(joker, 0), source_dv(joker, 0),
                                                                    phase(joker, 0), phase_dv(joker, 0), n,
                                                                    theta*DEG2RAD, dtheta_dv*DEG2RAD, eta*DEG2RAD,
                                                                    deta_dv*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_w)
      {
        attenuation_dw -= attenuation;
        attenuation_dw /= dB;
        phase_dw -= phase;
        phase_dw /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentDerivative(attenuation(joker, 0), attenuation_dw(joker, 0),
                                                                       phase(joker, 0), phase_dw(joker, 0), n, 
                                                                       theta*DEG2RAD, dtheta_dw*DEG2RAD, eta*DEG2RAD,
                                                                       deta_dw*DEG2RAD);
        
        if(do_src)
        {
          source_dw -= source;
          source_dw /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponentDerivative(source(joker, 0), source_dw(joker, 0), 
                                                                    phase(joker, 0), phase_dw(joker, 0), n,
                                                                    theta*DEG2RAD, dtheta_dw*DEG2RAD, eta*DEG2RAD,
                                                                    deta_dw*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_theta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentThetaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                            theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponentThetaDerivative(source(joker, 0), phase(joker, 0), n,
                                                                         theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_magnetic_eta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentEtaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                          theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentEtaDerivative(source(joker, 0), phase(joker, 0), n,
                                                                            theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials(iq)==JQT_temperature)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(attenuation(joker, 0), phase(joker, 0), dn_dT,
                                                             theta*DEG2RAD, eta*DEG2RAD);
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(partial_attenuation[iq](joker, 0),
                                                             partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                             eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(source(joker, 0), phase(joker, 0), dn_dT, theta*DEG2RAD,
                                                          eta*DEG2RAD, planck_BT);
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                          theta*DEG2RAD, eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanSigmaPlusComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                          eta*DEG2RAD, dplanck_BT(0, joker));
        }
      }
      else if(flag_partials(iq)==JQT_frequency)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(partial_attenuation[iq](joker, 0),
                                                             partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                             eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                          theta*DEG2RAD, eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanSigmaPlusComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                          eta*DEG2RAD, dplanck_BT(1, joker));
        }
      }
      else if(flag_partials(iq)==JQT_VMR)
      {  
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(attenuation(joker, 0), phase(joker, 0),
                                                             n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(source(joker, 0), phase(joker, 0),
                                                          n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD,
                                                          planck_BT);
        }
      }
      else if(flag_partials(iq) not_eq JQT_NOT_JQT)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(partial_attenuation[iq](joker, 0),
                                                             partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                             eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                          theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
    }
  }
}

/*
void xsec_species_line_mixing_wrapper_with_zeeman_replacement(PropagationMatrix& propmat_clearsky, 
                                                              StokesVector& nlte_source,
                                                              ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                                              ArrayOfStokesVector& dnlte_dx_source,
                                                              ArrayOfStokesVector& nlte_dsource_dx,
                                                              const ArrayOfArrayOfSpeciesTag& abs_species, 
                                                              const PropmatPartialsData& flag_partials,
                                                              const Index& abs_lineshape_ls, 
                                                              const Index& abs_lineshape_lsn, 
                                                              const Numeric& abs_lineshape_cutoff, 
                                                              const ArrayOfLineRecord& lr, 
                                                              const Vector&  Zeeman_DF,
                                                              const Vector&  planck_BT,
                                                              const Matrix&  dplanck_BT,
                                                              const SpeciesAuxData& isotopologue_ratios, 
                                                              const SpeciesAuxData& partition_functions,
                                                              const Vector& T_nlte, 
                                                              const Vector& abs_vmrs, 
                                                              const Numeric& P,
                                                              const Numeric& T, 
                                                              const Vector& f_grid,
                                                              const Vector& rtp_mag,
                                                              const Vector& r_path_los,
                                                              const Numeric& lm_p_lim,
                                                              const Numeric& theta, 
                                                              const Numeric& eta, 
                                                              const Numeric& H_mag, 
                                                              const Index& DM, 
                                                              const Index& this_species,
                                                              const Verbosity& verbosity )


{
  // Setup the size of the problem
  const Index nf = f_grid.nelem();
  const Index nq = flag_partials.nelem();
  const Index nl = lr.nelem();
  const Index nnlte = T_nlte.nelem();
  
  if(0 == nf or 0 == nl)
    return;
  
  // Compute variables for pressure broadening
  Numeric G0, G2, L0, L2, e, FVC;
  
  // Compute variables for pressure broadening derivatives
  Numeric dG0_dT, dG2_dT, dL0_dT, dL2_dT, de_dT, dFVC_dT;
  
  // Compute variables for line mixing
  Numeric Y, G, DV;
  
  // Compute variables for line mixing derivatives
  Numeric dY_dT, dG_dT, dDV_dT;
  
  // Compute variables for Line strength
  Numeric QT, QT0, K3, K4, Tu, Tl, r_low;
  
  // Compute variables for Line strength derivatives
  Numeric dQT_dT, dK1_dT, dK2_dT, dK2_dF0, dK3_dT, dK4_dT, 
    dK3_dF0, dK3_dTl, dK3_dTu, dK4_dTu;
  
  // Compute vectors attenuation/phase
  ComplexVector F, F_single_line(nf), N, N_single_line;  // Only single line lineshape needs to be initialized
  
  // Compute vectors derivative of attenuation/phase
  ArrayOfComplexVector dF, dF_single_line(nq);
  for(auto& cv : dF_single_line) 
    cv.resize(nf);
  
  ArrayOfComplexVector dN, dN_single_line(nq);
  if(nnlte)
  {
    N_single_line.resize(nf);
    
    for(auto& cv : dN_single_line) 
      cv.resize(nf);
  }
  
  // All should be of same species here
  const Index spec = lr[0].Species();
  const Index isop = lr[0].Isotopologue();
  const Numeric T0 = lr[0].Ti0();
  const Numeric isotopic_ratio = isotopologue_ratios.getParam(spec, isop)[0].data[0];
  
  // Partial pressure of this species
  const Numeric P_partial = abs_vmrs[this_species] * P;
                          
  // Number density of molecules of the level
  const Numeric nd = abs_vmrs[this_species] * number_density(P, T); // NB: OF SPECIES NOT ISOTOPOLOGUE
  
  // Doppler broadening parameters
  const Numeric GD_div_F0 = Linefunctions::DopplerConstant(T, lr[0].IsotopologueData().Mass());
  const Numeric dGD_div_F0_dT = Linefunctions::dDopplerConstant_dT(T, lr[0].IsotopologueData().Mass());
  
  // Broadening species
  ArrayOfIndex broad_spec_locations;
  find_broad_spec_locations(broad_spec_locations, abs_species, this_species);
  
  // Water index
  const Index h2o_index = find_first_species_tg(abs_species, species_index_from_species_name("H2O"));
  
  // Pressure broadening derivatives
  ComplexVector dgamma;
  
  partition_function(QT0, QT, lr[0].Ti0(), T, 
                     partition_functions.getParamType(spec, isop), 
                     partition_functions.getParam(spec, isop));
  
  if(flag_partials.do_temperature())
  {
    dpartition_function_dT( dQT_dT,
                            QT,
                            T,
                            flag_partials.Temperature_Perturbation(), 
                            partition_functions.getParamType(spec, isop), 
                            partition_functions.getParam(spec, isop));
  }
  
  for(Index il = 0; il < nl; il++)
  {
    const LineRecord& l = lr[il];
    
    // Need to check that Ti0, Species, and Isotopologue of line has not changed...
    assert(T0 == l.Ti0());
    assert(spec == l.Species());
    assert(isop == l.Isotopologue());
    
    const Numeric& F0 = l.F();
    const Numeric& S0 = l.I0();
    const Numeric& E0 = l.Elow();
    
    const Numeric gamma = stimulated_emission(T, F0);
    const Numeric gamma_ref = stimulated_emission(T0, F0);
    const Numeric K1 = boltzman_ratio(T, T0, E0);
    const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
    
    // QuantumIdentifier
    const QuantumIdentifier QI = lr[il].QuantumIdentity(); 
    
    l.PressureBroadening().GetPressureBroadeningParams(G0, G2, e, L0, L2, FVC,
                                                       l.Ti0()/T, P, P_partial, 
                                                       this_species, h2o_index,
                                                       broad_spec_locations,
                                                       abs_vmrs, verbosity);
    
    // Get line mixing parameters of this line in this
    l.LineMixing().GetLineMixingParams(Y,  G,  DV, T, P, lm_p_lim, 1);
    
    // Treat line mixing as a line center shift
    L0 += DV;
    
    if(flag_partials.do_temperature())
    {
      l.PressureBroadening().GetPressureBroadeningParams_dT(dG0_dT, dL0_dT,
                                                            T, l.Ti0(), P, P_partial, 
                                                            this_species, h2o_index,
                                                            broad_spec_locations,
                                                            abs_vmrs, verbosity);
      
      if(flag_partials.get_first_pressure_term() > -1)
        l.PressureBroadening().SetInternalDerivatives(dgamma, flag_partials, QI, l.Ti0()/T,
                                                      P, P_partial, this_species, h2o_index,
                                                      abs_vmrs, verbosity);
      
      l.LineMixing().GetLineMixingParams_dT(dY_dT,  dG_dT,  dDV_dT, T, 
                                            flag_partials.Temperature_Perturbation(),
                                            P, lm_p_lim, 1);
      
      dL0_dT += dDV_dT;
      
      dK1_dT = dboltzman_ratio_dT(K1, T, E0);
      dK2_dT = dstimulated_relative_emission_dT(gamma, gamma_ref, F0);
    }
    
    if(flag_partials.do_line_center())
    {
      dK2_dF0 = dstimulated_relative_emission_dF0(gamma, gamma_ref, T);
    }
      
    switch(l.PressureBroadening().Type())
    {
      case PressureBroadeningData::PB_SD_AIR_VOLUME:
        Linefunctions::set_htp(F_single_line, dF_single_line,
                               f_grid, Zeeman_DF[il], H_mag, F0,
                               GD_div_F0, G0, L0, L2, G2, e, FVC,
                               flag_partials, QI, 
                               dGD_div_F0_dT, dG0_dT, 
                               dL0_dT, dG2_dT, dL2_dT, de_dT, dFVC_dT);
        break;
      case PressureBroadeningData::PB_AIR_AND_WATER_BROADENING:
      case PressureBroadeningData::PB_PLANETARY_BROADENING:
      case PressureBroadeningData::PB_AIR_BROADENING:
        Linefunctions::set_faddeeva_algorithm916(F_single_line, dF_single_line, 
                                                 f_grid, Zeeman_DF[il], H_mag, F0,
                                                 GD_div_F0, G0, L0, 
                                                 flag_partials, QI,
                                                 dGD_div_F0_dT, dG0_dT, dL0_dT);
        break;
      default:
        throw std::runtime_error("Developer messed up");
    }
    
    
    Linefunctions::apply_linemixing(F_single_line, dF_single_line,
                                    Y, G, flag_partials, QI, dY_dT, dG_dT);
    
    Linefunctions::apply_VVH(F_single_line, dF_single_line, f_grid, F0, T, flag_partials);
    
    Linefunctions::apply_linestrength(F_single_line, dF_single_line, 
                                      S0, isotopic_ratio,
                                      QT, QT0, K1, K2,
                                      flag_partials, QI,
                                      dQT_dT, dK1_dT, dK2_dT, dK2_dF0);
    
    
    if(flag_partials.get_first_pressure_term() > -1)
    {
      Linefunctions::apply_pressurebroadening_jacobian(dF_single_line, flag_partials, QI, dgamma);
    }
    
    if(nnlte)
    {
      const Index evlow_index = l.EvlowIndex();
      const Index evupp_index = l.EvuppIndex();
      const Numeric El = l.Evlow();
      const Numeric Eu = l.Evupp();
      
      if(evlow_index > -1)
      {
        Tu = T_nlte[evupp_index];
        K4 = boltzman_ratio(Tu, T, Eu);
      }
      else
      {
        Tu = T;
        K4 = 1.0;
      }
      
      if(evlow_index > -1)
      {
        Tl = T_nlte[evlow_index];
        r_low = boltzman_ratio(Tl, T, El);
      }
      else 
      {
        Tl = T;
        r_low = 1.0;
      }
      
      K3 = absorption_nlte_ratio(gamma, K4, r_low);
      
      if(flag_partials.do_frequency())
      {
        dK3_dF0 = dabsorption_nlte_rate_dF0(gamma, T, K4, r_low);
      }
      
      if(flag_partials.do_temperature())
      {
        dK3_dT = dabsorption_nlte_rate_dT(gamma, T, F0, El, Eu, K4, r_low);
      }
      
      if(El > 0)
      {
        dK3_dTl = dabsorption_nlte_rate_dTl(gamma, T, Tl, El, r_low);
      }
      
      if(Eu > 0)
      {
        dK3_dTu = dabsorption_nlte_rate_dTu(gamma, T, Tu, Eu, K4);
        dK4_dTu = dboltzman_ratio_dT(K4, Tu, Eu);
      }
      
      Linefunctions::apply_nonlte(F_single_line, dF_single_line, N_single_line, dN_single_line, K3, K4, 
                                  flag_partials, QI, dK3_dT, dK4_dT, dK3_dF0, dK3_dTl, dK3_dTu, dK4_dTu);
    }
    
    if(il)
    {
      F += F_single_line;
      
      if(nnlte)
        N += N_single_line;
      
      for(Index iq = 0; iq < nq; iq++)
      {
        dF[iq] += dF_single_line[iq];
        
        if(nnlte)
          dN[iq] += dN_single_line[iq];
      }
    }
    else
    {
      F = F_single_line;
      
      dF = dF_single_line;
      
      if(nnlte)
        N = N_single_line;
      
      if(nnlte)
        dN = dN_single_line;
    }
  }
  
  
}
*/

void set_magnetic_parameters(Numeric& H_mag,
                             Numeric& eta,
                             Numeric& theta,
                             const Index manual_zeeman_tag,
                             const Numeric& manual_zeeman_eta,
                             const Numeric& manual_zeeman_theta,
                             const Numeric& manual_zeeman_magnetic_field_strength,
                             ConstVectorView rtp_mag,
                             ConstVectorView r_path_los)
{
//Get the magnitude of the magnetic field and store a local unit Vector for simplified angle calculations.
  H_mag = manual_zeeman_tag != 0?manual_zeeman_magnetic_field_strength:sqrt( rtp_mag * rtp_mag );

  if(manual_zeeman_tag!=0)
  { // Leaving it up to the user to manually tag the angles for simplified magnetic fields.
    eta   = manual_zeeman_eta;
    theta = manual_zeeman_theta;
  }
  else if(H_mag==0.0)
  {
      eta = 0.0;
      theta = 0.0;
  }
  else
  { 
    const Numeric 
    aa=DEG2RAD*r_path_los[1], 
    za=DEG2RAD*r_path_los[0], 
    Bu = rtp_mag[0], 
    Bv = rtp_mag[1], 
    Bw = rtp_mag[2],
    cosaa=cos(aa),
    cosza=cos(za),
    sinaa=sin(aa),
    sinza=sin(za),
    cosaa2 = cosaa*cosaa,
    //cosza2 = cosza*cosza,
    sinza2 = sinza*sinza,
    sinaa2 = sinaa*sinaa,
    Bu2 = Bu*Bu,
    Bv2 = Bv*Bv,
    Bw2 = Bw*Bw,
    H2  = H_mag*H_mag,
    term1=(Bu*(cosaa2*sinza2 - 1.) + Bw*cosaa*cosza*sinza + Bv*cosaa*sinaa*sinza2),
    term2=(Bv*(sinaa2*sinza2 - 1.) + Bw*sinaa*cosza*sinza + Bu*cosaa*sinaa*sinza2),
    term3=( - Bw*sinza2            + Bu*cosaa*cosza*sinza + Bv*sinaa*cosza*sinza),
    eta_test=(PI - acos(((Bu*cosaa*cosza - Bw*sinza + Bv*sinaa*cosza)*sqrt(term1*term1 + term2*term2 +term3*term3))/(H2)))*RAD2DEG,
    x1 = (Bv*cosaa - Bu*sinaa),
    x2 = -((Bu2*cosaa-Bv2*cosaa+2.0*Bu*Bv*sinaa)*cosaa*sinza2 - Bu2 - Bw2 + 2.0*Bu*Bw*cosaa*cosza*sinza + (Bw2*cosza - Bv2*cosza+ 2.0*Bv*Bw*sinaa*sinza)*cosza),
    fx2=sqrt(x2),
    x3 = 1.0/H2,
    x=x1*fx2*x3;
    
    theta = acos((Bw*cosza + Bu*cosaa*sinza + Bv*sinaa*sinza)/H_mag) * RAD2DEG;
    
    if ((abs(x)-(Numeric)1.0)>0.0) // Numerical drifts can cause this...
      eta = 0.0;
    else 
    {
      eta=acos(x)*RAD2DEG;
      
      if(eta_test>90.0) eta*=-1.0;
    }
  }
}


void set_magnetic_parameters_derivative(
    Numeric& dH_du,
    Numeric& dH_dv,
    Numeric& dH_dw,
    Numeric& deta_du,
    Numeric& deta_dv,
    Numeric& deta_dw,
    Numeric& dtheta_du,
    Numeric& dtheta_dv,
    Numeric& dtheta_dw,
    ConstVectorView rtp_mag,
    ConstVectorView r_path_los)
{
    const Numeric& Bu = rtp_mag[0],Bv = rtp_mag[1],Bw = rtp_mag[2];
    
    const Numeric
    Bu2=Bu*Bu, Bv2 = Bv*Bv, Bw2 = Bw*Bw, 
    H_mag = sqrt(Bu2+Bv2+Bw2), H2 = H_mag*H_mag,
    aa=DEG2RAD*r_path_los[1], za=DEG2RAD*r_path_los[0], 
    cosaa=cos(aa), cosza=cos(za), sinaa=sin(aa), sinza=sin(za),
    sinza2=sinza*sinza,cosza2=cosza*cosza,cosaa2=cosaa*cosaa,sinaa2=sinaa*sinaa,
    a=(Bw*cos(za) + Bu*cos(aa)*sin(za) + Bv*sin(aa)*sin(za)),
    b=-cosaa*sinza*Bv2 + Bu*sinaa*sinza*Bv - cosaa*sinza*Bw2 + Bu*cosza*Bw,
    term1=(Bu*(cosaa2*sinza2 - 1.) + Bw*cosaa*cosza*sinza + Bv*cosaa*sinaa*sinza2),
    term2=(Bv*(sinaa2*sinza2 - 1.) + Bw*sinaa*cosza*sinza + Bu*cosaa*sinaa*sinza2),
    term3=(Bw*(cosza2 - 1.) + Bu*cosaa*cosza*sinza + Bv*sinaa*cosza*sinza),
    eta_test=(PI - acos(((Bu*cosaa*cosza - Bw*sinza + 
    Bv*sinaa*cosza)*sqrt(term1*term1/(H_mag*H_mag) + term2*term2/(H_mag*H_mag) +
    term3*term3/(H_mag*H_mag)))/(H_mag)))*RAD2DEG; 
    
    const Numeric
    x1 = (Bv*cosaa - Bu*sinaa),
    x2 = -((Bu2*cosaa-Bv2*cosaa+2.0*Bu*Bv*sinaa)*cosaa*sinza2
    - Bu2 - Bw2 + 2.0*Bu*Bw*cosaa*cosza*sinza +
    (Bw2*cosza - Bv2*cosza+ 2.0*Bv*Bw*sinaa*sinza)*cosza),
    fx2=sqrt(x2),
    x3 = 1.0/H2,
    x=x1*fx2*x3, // eta = acos(x) * RAD2DEG;
    dx1_dBu = -sinaa,
    dx1_dBv = cosaa,
    dx2_dBu = 2.0*Bu - cosaa*sinza2*(2.0*Bu*cosaa + 2.0*Bv*sinaa) - 2.0*Bw*cosaa*cosza*sinza,
    dx2_dBv = cosza*(2.0*Bv*cosza - 2.0*Bw*sinaa*sinza) + cosaa*sinza2*(2.0*Bv*cosaa - 2.0*Bu*sinaa),
    dx2_dBw = 2.0*Bw - cosza*(2.0*Bw*cosza + 2.0*Bv*sinaa*sinza) - 2.0*Bu*cosaa*cosza*sinza,
    dx3_dBu = -2.0*Bu*x3*x3,
    dx3_dBv = -2.0*Bv*x3*x3,
    dx3_dBw = -2.0*Bw*x3*x3,
    c1 = fx2*x3, c2 = x1*0.5*x3/fx2,c3=x1*fx2,
    d_acos  = -1/sqrt(1-x*x) * (eta_test>90.0?-RAD2DEG:RAD2DEG); //RAD2DEG is here for convenience...
    
    dH_du = Bu/H_mag;
    dH_dv = Bv/H_mag;
    dH_dw = Bw/H_mag;
    
    dtheta_du = b/sqrt(H_mag*H_mag-a*a)*H_mag*H_mag * RAD2DEG;
    dtheta_dv = dtheta_du * (- sinaa*sinza*Bu2 + Bv*cosaa*sinza*Bu - sinaa*sinza*Bw2 + Bv*cosza*Bw)/b;
    dtheta_dw = dtheta_du * (-(cosza*Bu2 - Bw*cosaa*sinza*Bu + cosza*Bv2 - Bw*sinaa*sinza*Bv)/b);
    
    deta_du = d_acos * (dx1_dBu*c1 + dx2_dBu*c2 +  dx3_dBu*c3);
    deta_dv = d_acos * (dx1_dBv*c1 + dx2_dBv*c2 +  dx3_dBv*c3);
    deta_dw = d_acos * (/* empty */  dx2_dBw*c2 +  dx3_dBw*c3);
}




void set_quantumnumbers( LineRecord& this_LR,
                         const Rational& hund,
                         const SpeciesAuxData& isotopologue_quantum)
{
  Rational Lambda, S;
  
  S = Rational((Index) (2*isotopologue_quantum.getParam(this_LR.Species(), 
                                                        this_LR.Isotopologue())[0].data[AuxIndex_S]),2);
  S.Simplify();
  
  Lambda = (Index) isotopologue_quantum.getParam(this_LR.Species(), 
                                         this_LR.Isotopologue())[0].data[AuxIndex_Lambda];
  
  this_LR.SetQuantumNumberLower(QN_S, S);
  this_LR.SetQuantumNumberUpper(QN_S, S);
  
  this_LR.SetQuantumNumberLower(QN_Lambda, Lambda);
  this_LR.SetQuantumNumberUpper(QN_Lambda, Lambda);
  
  // Flag quantum number
  this_LR.SetQuantumNumberLower(QN_Hund, hund);
  this_LR.SetQuantumNumberUpper(QN_Hund, hund);
}

void alter_linerecord( LineRecord& new_LR,
                       Numeric& Test_RS,
                       const Numeric& old_LS,
                       const Rational& J_up,
                       const Rational& J_lo,
                       const Rational& M_up,
                       const Rational& M_lo)
{
    
    // Test that we did not mess up somewhere
    assert(abs(M_lo)<=J_lo);
    assert(abs(M_up)<=J_up);
    
    // Find the relative strength
    const Numeric RS = relative_strength(M_lo, J_lo, (J_up-J_lo).toIndex(), (M_up-M_lo).toIndex());
    
    // Setup a test that the relative strength is reasonable
    Test_RS += RS;
    new_LR.setI0( old_LS * RS );
    
    // Set quantum numbers
    new_LR.SetQuantumNumberLower(QN_M, M_lo);
    new_LR.SetQuantumNumberUpper(QN_M, M_up);
}



void create_Zeeman_linerecordarrays(
        ArrayOfArrayOfLineRecord& aoaol,
        ArrayOfVector& z1_frequencyshift,
        const ArrayOfArrayOfSpeciesTag& abs_species,
        const ArrayOfArrayOfLineRecord& abs_lines_per_species,
        const SpeciesAuxData& isotopologue_quantum,
        const Verbosity& verbosity)
{
    CREATE_OUT3;
    // Note that this function assumes that all tests that are not line specifice are done elsewhere
    
    using global_data::species_data;
    
    // This margin is for relative strength and can perhaps be lowered by returning RS as Rational?
    const Numeric margin    = 1e-4; 

      // For all species
      for(Index II = 0; II<abs_species.nelem(); II++)
      {
          // If the species isn't Zeeman, look at the next species
          if(!is_zeeman(abs_species[II])) continue;
          
          // If there are no lines give up on this species
          if(!abs_lines_per_species[II].nelem()) continue;
          
          const Index nis = isotopologue_quantum.nisotopologues(abs_lines_per_species[II][0].Species());
          
          const SpeciesRecord& sr = species_data[abs_lines_per_species[II][0].Species()];
          
          for(Index is=0;is<nis;is++)
          {
              const IsotopologueRecord& ir = sr.Isotopologue()[is];
              
              if(ir.isContinuum())
                  continue;
              
              const ArrayOfGriddedField1& aogf1 =
              isotopologue_quantum.getParam(abs_lines_per_species[II][0].Species(),is);
              if(!aogf1.nelem())
              {
                  std::ostringstream os;
                  os<<"No data in isotopologue_quantum for species "<< sr.Name();
              }
              
              const GriddedField1& gf1 = aogf1[0];
              if(gf1.data.nelem()!=AuxIndex_TotalCount)
              {
                  std::ostringstream os;
                  os << "There are undefined isotopologues in *isotopologue_quantum* for species " <<sr.Name() <<"\n";
                  os << "All isotopologues must be defined for a species.\n";
                  throw std::runtime_error(os.str());
              }
          }
                  
          // One line record array per type of polarizer is created
          aoaol.push_back(ArrayOfLineRecord()); // First is negative
          aoaol.push_back(ArrayOfLineRecord()); // Second is 0
          aoaol.push_back(ArrayOfLineRecord()); // Third is positive

          ArrayOfLineRecord& temp_abs_lines_sm = aoaol[aoaol.nelem()-3]; // sigma minus
          ArrayOfLineRecord& temp_abs_lines_pi = aoaol[aoaol.nelem()-2]; // pi
          ArrayOfLineRecord& temp_abs_lines_sp = aoaol[aoaol.nelem()-1]; // sigma plus

          temp_abs_lines_sm.resize(0);
          temp_abs_lines_sp.resize(0);
          temp_abs_lines_pi.resize(0);

          temp_abs_lines_sm.reserve(25000);
          temp_abs_lines_sp.reserve(25000);
          temp_abs_lines_pi.reserve(25000);

          // Else loop over all the lines in the species.
          for (Index ii = 0; ii< abs_lines_per_species[II].nelem(); ii++)
          {
              // local LineRecord
              LineRecord temp_LR = abs_lines_per_species[II][ii];
              const Numeric this_linestrength = temp_LR.I0();
              Rational hund;
                
              set_hund_case(hund, isotopologue_quantum, temp_LR);
              Numeric RS_sum     = 0; //Sum relative strength (which ought be close to one by the end)
              // Only look at lines which have no change in the main rotational number
              
              bool test=true;
              
              // Test that hund cases are properly defined
              if(hund==0)
              {
                  if(    temp_LR.QuantumNumbers().Upper()[QN_Omega].isUndefined() 
                      || temp_LR.QuantumNumbers().Lower()[QN_Omega].isUndefined())
                      test=false;
              }
              else if(hund==1)
              {
                  if(    temp_LR.QuantumNumbers().Upper()[QN_N].isUndefined() 
                      || temp_LR.QuantumNumbers().Lower()[QN_N].isUndefined())
                      test=false;
              }
              else //Case a is 0, case b is 1
              {
                  std::ostringstream os;
                  os << "There are undefined Hund cases for\n" << temp_LR;
                  throw std::runtime_error(os.str());
              }
              
              if(!test)
              {
                  std::ostringstream os;
                  os<<"Need to define the main quantum numbers for\n"<<
                  temp_LR;
                  throw std::runtime_error(os.str());
              }
              
              // Quantum numbers
              set_quantumnumbers(temp_LR, hund, isotopologue_quantum);
              const Rational J_up = temp_LR.QuantumNumbers().Upper()[QN_J], 
              J_lo = temp_LR.QuantumNumbers().Lower()[QN_J];
            
              test = J_up.isUndefined()||J_lo.isUndefined();
              if(test)
              {
                  std::ostringstream os;
                  os<<"Need to define J for\n"<<
                  temp_LR;
                  throw std::runtime_error(os.str());
              }
              
              test      = J_lo>J_up;
              const bool same = J_lo==J_up;
              const Rational  J    = (test?J_up:J_lo);
              
              for ( Rational M = -J; M<=J; M++ )
              {
                  /*
                   *                              Note that:
                   *                              sp := sigma plus,  which means DM =  1
                   *                              sm := sigma minus, which means DM = -1
                   *                              pi := planar,      which means DM =  0
                   */
                  if(test&&!same)
                  {
                      alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M, M-1  );
                      temp_abs_lines_sm.push_back(temp_LR);
                      
                      alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M, M    );
                      temp_abs_lines_pi.push_back(temp_LR);
                      
                      alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M, M+1  );
                      temp_abs_lines_sp.push_back(temp_LR);
                  }
                  else if(same)
                  {
                      alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M, M    );
                      temp_abs_lines_pi.push_back(temp_LR);
                      
                      if(M==-J)
                      {
                          alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M, M+1  );
                          temp_abs_lines_sp.push_back(temp_LR);
                      }
                      else if(M==J)
                      {
                          alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M, M-1  );
                          temp_abs_lines_sm.push_back(temp_LR);
                      }
                      else
                      {
                          alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M, M-1  );
                          temp_abs_lines_sm.push_back(temp_LR);                          
                          
                          alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M, M+1  );
                          temp_abs_lines_sp.push_back(temp_LR);
                      }
                  }
                  else
                  {
                      alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M+1, M  );
                      temp_abs_lines_sm.push_back(temp_LR);
                      
                      alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M  , M  );
                      temp_abs_lines_pi.push_back(temp_LR);
                      
                      alter_linerecord( temp_LR, RS_sum, this_linestrength, J_up, J_lo, M-1, M  );
                      temp_abs_lines_sp.push_back(temp_LR);
                  }
              }
              
              if (abs(RS_sum-1.)>margin) //Reasonable confidence?
              {
                  std::ostringstream os;
                  os << "The sum of relative strengths is not close to one. This is problematic and "
                  "you should look into why this happens.\nIt is currently " << RS_sum 
                  << " with J_lo: "<<J_lo<<", J_up: "<<J_up<<" for line: "<<
                  temp_LR <<"\n";
                  throw std::runtime_error(os.str());
              }
          }
      }
      
      const Index nzeeman = aoaol.nelem();
      z1_frequencyshift.resize(nzeeman);
      for(Index I=0;I<nzeeman;I++)
      {
          const Index nlines = aoaol[I].nelem();
          z1_frequencyshift[I].resize(nlines);
          for(Index J=0;J<nlines;J++)
          {
              Numeric GS;
              
              // Set necessary parameters from isotopologue_quantum
              set_GS(GS, isotopologue_quantum, aoaol[I][J]);
              
              // Store the frequency shift
              z1_frequencyshift[I][J] = frequency_change(aoaol[I][J], 1, GS);
          }
              
      }
      
}



void set_hund_case(Rational& hund, const SpeciesAuxData& isotopologue_quantum,const LineRecord& temp_LR)
{
  hund = (Index) isotopologue_quantum.getParam(temp_LR.Species(), 
                                               temp_LR.Isotopologue())[0].data[AuxIndex_Hund];
}

void set_GS(Numeric& GS, const SpeciesAuxData& isotopologue_quantum,const LineRecord& temp_LR)
{
    GS = isotopologue_quantum.getParam(temp_LR.Species(), 
                                       temp_LR.Isotopologue())[0].data[AuxIndex_GS];
}


void set_strength_partial_derivative_matrix(ArrayOfMatrix& A, ArrayOfMatrix& B, const ArrayOfRetrievalQuantity flag_partials, const Vector& f_grid)
{
    for(Index ii=0; ii<flag_partials.nelem(); ii++)
    {
        if(flag_partials[ii].MainTag() == "Zeeman" &&  flag_partials[ii].Subtag() == "Magnetic Strength" && flag_partials[ii].SubSubtag() == "From Propagation")
        {
            A[ii].resize(f_grid.nelem(),1);
            B[ii].resize(f_grid.nelem(),1);
            A[ii]=0;
            B[ii]=0;
            return;
        }
    }
}

Index part_mag_strength(const ArrayOfRetrievalQuantity& flag_partials)
{
    for(Index ii=0; ii<flag_partials.nelem(); ii++)
        if(flag_partials[ii].MainTag() == "Zeeman" &&  flag_partials[ii].Subtag() == "Magnetic Strength" && flag_partials[ii].SubSubtag() == "From Propagation")
            return ii;
    return -1;
}

Index part_mag_theta(const ArrayOfRetrievalQuantity& flag_partials)
{
    for(Index ii=0; ii<flag_partials.nelem(); ii++)
        if(flag_partials[ii].MainTag() == "Zeeman" &&  flag_partials[ii].Subtag() == "Magnetic Theta" && flag_partials[ii].SubSubtag() == "From Propagation")
            return ii;
    return -1;
}
