/**
   \file   continua.h

   Contains declarations of continuum functions.

   \author Stefan Buehler
   \date   2001-01-17
*/

#ifndef contiua_h
#define contiua_h

#include "matpackI.h"

const Numeric VMRCalcLimit = 1.000e-25;


void xsec_continuum_tag( MatrixView                    xsec,
			 const String&              name,
			 ConstVectorView              parameters,
			 ConstVectorView  	    f_mono,
			 ConstVectorView  	    p_abs,
			 ConstVectorView  	    t_abs,
			 ConstVectorView  	    n2_abs,
			 ConstVectorView  	    h2o_abs,
			 ConstVectorView              vmr );


void check_continuum_model(const String& name);


//////////////////////////////////////////////////////////////////////////// 
// water vapor absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void Rosenkranz_H2O_self_continuum( MatrixView           xsec,
				    Numeric	      C,
				    Numeric	      x,
				    ConstVectorView     f_mono,
				    ConstVectorView     p_abs,
				    ConstVectorView     t_abs,
				    ConstVectorView     vmr	 );

void Rosenkranz_H2O_foreign_continuum( MatrixView           xsec,
				       Numeric	         C,
				       Numeric	         x,
				       ConstVectorView     f_mono,
				       ConstVectorView     p_abs,
				       ConstVectorView     t_abs,
				       ConstVectorView     vmr	 );

void MPM93_H2O_continuum( MatrixView           xsec,
			  ConstVectorView     f_mono,
			  ConstVectorView     p_abs,
			  ConstVectorView     t_abs,
			  ConstVectorView     vmr	 );

//////////////////////////////////////////////////////////////////////////// 
// water droplet and ice particle absorption
//////////////////////////////////////////////////////////////////////////// 

void MPM93WaterDropletAbs( MatrixView           xsec,
			   ConstVectorView   f_mono,  // frequency vector
			   ConstVectorView    p_abs,  // pressure vector
			   ConstVectorView    t_abs,  // temperature vector
			   ConstVectorView      vmr); // suspended water droplet density vector

void MPM93IceCrystalAbs( MatrixView           xsec,
			 ConstVectorView   f_mono,    // frequency vector
			 ConstVectorView    p_abs,    // pressure vector
			 ConstVectorView    t_abs,    // temperature vector
			 ConstVectorView      vmr	 ); // suspended ice particle density vector, 
                                                    // valid range: 0-10.0e-3 kg/m3

//////////////////////////////////////////////////////////////////////////// 
// oxygen absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void MPM93O2AbsModel( MatrixView           xsec,
		      ConstVectorView     f_mono,
		      ConstVectorView     p_abs,
		      ConstVectorView     t_abs,
		      ConstVectorView     h2o_abs,
		      ConstVectorView     vmr );

void MPM93_O2_continuum( MatrixView           xsec,
			 ConstVectorView     f_mono,
			 ConstVectorView     p_abs,
			 ConstVectorView     t_abs,
			 ConstVectorView     h2o_abs,
			 ConstVectorView     vmr	 );

void Rosenkranz_O2_continuum( MatrixView           xsec,
			      ConstVectorView  	f_mono,
			      ConstVectorView  	p_abs,
			      ConstVectorView  	t_abs,
			      ConstVectorView     h2o_abs,
			      ConstVectorView     vmr	 );

//////////////////////////////////////////////////////////////////////////// 
// nitrogen absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void MPM93_N2_continuum( MatrixView           xsec,
			 ConstVectorView     f_mono,
			 ConstVectorView     p_abs,
			 ConstVectorView     t_abs,
			 ConstVectorView     h2o_abs,
			 ConstVectorView     vmr	 );

void Rosenkranz_N2_self_continuum( MatrixView           xsec,
				   ConstVectorView     f_mono,
				   ConstVectorView     p_abs,
				   ConstVectorView     t_abs,
				   ConstVectorView     vmr    );

void General_N2_self_continuum(    MatrixView           xsec,
                                   Numeric           C,
                                   Numeric           xf,
                                   Numeric           xt,
                                   Numeric           xp,
				   ConstVectorView     f_mono,
				   ConstVectorView     p_abs,
				   ConstVectorView     t_abs,
				   ConstVectorView     vmr    );

//////////////////////////////////////////////////////////////////////////// 
// carbon dioxide absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void Rosenkranz_CO2_self_continuum( MatrixView           xsec,
				    ConstVectorView     f_mono,
				    ConstVectorView     p_abs,
				    ConstVectorView     t_abs,
				    ConstVectorView     vmr	 );

void Rosenkranz_CO2_foreign_continuum( MatrixView           xsec,
				       ConstVectorView     f_mono,
				       ConstVectorView     p_abs,
				       ConstVectorView     t_abs,
				       ConstVectorView     n2_abs,
				       ConstVectorView     vmr	 );

//////////////////////////////////////////////////////////////////////////// 
// help functions
//////////////////////////////////////////////////////////////////////////// 

Numeric MPMLineShapeFunction( Numeric gamma, 
			      Numeric fl, 
			      Numeric f);

Numeric MPMLineShapeO2Function( Numeric gamma, 
				Numeric fl, 
				Numeric f,
                                Numeric delta);

Numeric WVSatPressureLiquidWater(Numeric t);

Numeric WVSatPressureIce(Numeric t);

#endif // contiua_h
