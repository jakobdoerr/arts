/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

/** \file
    Declarations required for the calculation of absorption coefficients.

    \author Stefan Buehler
*/

// This file contains 

#ifndef absorption_h
#define absorption_h

#include <iostream>
#include "vecmat.h"

/** The type that is used to store pointers to lineshape
    functions.  */
typedef void (*lsf_type)(VECTOR&,
			 VECTOR&,
			 Numeric,
			 Numeric,
			 Numeric,
			 VECTOR::subrange_type,
			 const INDEX);

/** Lineshape related information. There is one LineshapeRecord for
    each available lineshape function.

    \author Stefan Buehler
    \date   2000-08-21  */
class LineshapeRecord{
public:

  /** Default constructor. */
  LineshapeRecord(){};

  /** Initializing constructor, used to build the lookup table. */
  LineshapeRecord(const string& name,
		  const string& description,
		  lsf_type      function)
    : mname(name),
      mdescription(description),
      mfunction(function)
  { /* Nothing to do here. */ }
  /** Return the name of this lineshape. */
  const string&  Name()        const { return mname;        }   
  /** Return the description text. */
  const string&  Description() const { return mdescription; }
  /** Return pointer to lineshape function. */
  lsf_type Function() const { return mfunction; }
private:	
  string  mname;        ///< Name of the function (e.g., Lorentz).
  string  mdescription; ///< Short description.
  lsf_type mfunction;   ///< Pointer to lineshape function.

};

/** The type that is used to store pointers to lineshape
    normalization functions.  */
typedef void (*lsnf_type)(VECTOR&,
			  Numeric,
			  VECTOR::subrange_type,
			  const INDEX);

/** Lineshape related normalization function information. There is one
    LineshapeNormRecord for each available lineshape normalization
    function.

    \author Axel von Engeln
    \date   2000-11-30  */
class LineshapeNormRecord{
public:

  /** Default constructor. */
  LineshapeNormRecord(){};

  /** Initializing constructor, used to build the lookup table. */
  LineshapeNormRecord(const string& name,
		      const string& description,
		      lsnf_type      function)
    : mname(name),
      mdescription(description),
      mfunction(function)
  { /* Nothing to do here. */ }
  /** Return the name of this lineshape. */
  const string&  Name()        const { return mname;        }   
  /** Return the description text. */
  const string&  Description() const { return mdescription; }
  /** Return pointer to lineshape normalization function. */
  lsnf_type Function() const { return mfunction; }
private:	
  string  mname;        ///< Name of the function (e.g., linear).
  string  mdescription; ///< Short description.
  lsnf_type mfunction;  ///< Pointer to lineshape normalization function.
};

/** Lineshape related specification like which lineshape to use, the
normalizationfactor, and the cutoff.

    \author Axel von Engeln
    \date   2001-01-05  */
class LineshapeSpec{
public:

  /** Default constructor. */
  LineshapeSpec(){};

  /** Initializing constructor. */
  LineshapeSpec(const size_t&    ind_ls,
		const size_t&    ind_lsn,
		const Numeric&   cutoff)
    : mind_ls(ind_ls),
      mind_lsn(ind_lsn),
      mcutoff(cutoff)
  { /* Nothing to do here. */ }

  /** Return the index of this lineshape. */
  const size_t&  Ind_ls()        const { return mind_ls; }   
  /** Set it. */
  void SetInd_ls( size_t ind_ls ) { mind_ls = ind_ls; }

  /** Return the index of the normalization factor. */
  const size_t&  Ind_lsn()       const { return mind_lsn; }
  /** Set it. */
  void SetInd_lsn( size_t ind_lsn ) { mind_lsn = ind_lsn; }

  /** Return the cutoff frequency (in Hz). This is the distance from
      the line center outside of which the lineshape is defined to be
      zero. Negative means no cutoff.*/
  const Numeric& Cutoff() const { return mcutoff; }
  /** Set it. */
  void SetCutoff( Numeric cutoff ) { mcutoff = cutoff; }
private:	
  size_t  mind_ls;
  size_t  mind_lsn;
  Numeric mcutoff;
};

/** Holds a list of lineshape specifications: function, normalization, cutoff.
    \author Axel von Engeln */
typedef ARRAY<LineshapeSpec> ARRAYofLineshapeSpec;



/** Contains the lookup data for one isotope.
    \author Stefan Buehler */
class IsotopeRecord{
public:
  /** Default constructor. Needed by make_array. */
  IsotopeRecord() { /* Nothing to do here */ }
  /** Constructor that sets the values. */
  IsotopeRecord(const string&  	  name,
		const Numeric& 	  abundance,
		const Numeric& 	  mass,
		const int&     	  mytrantag,
		const int&     	  hitrantag,
		const ARRAY<int>& jpltags)
    : mname(name),
      mabundance(abundance),
      mmass(mass),
      mmytrantag(mytrantag),
      mhitrantag(hitrantag),
      mjpltags(jpltags.size())
  {
    // We need to use copy to initialize the ARRAY members. If we use
    // the assignment operator they end up all pointing to the same
    // data!
    copy(jpltags,mjpltags);

    // Some consistency checks whether the given data makes sense.
#ifndef NDEBUG
      {
	/* 1. All the tags must be positive or -1 */
	assert( (0<mmytrantag) || (-1==mmytrantag) );
	assert( (0<mhitrantag) || (-1==mhitrantag) );
	for ( size_t i=0; i<mjpltags.size(); ++i )
	  assert( (0<mjpltags[i]) || (-1==mjpltags[i]) );
      }
#endif // ifndef NDEBUG
  }

  /** Isotope name. */
  const string&       Name()         const { return mname;  }
  /** Normal abundance ( = isotopic ratio). (Absolute number.) */
  const Numeric&      Abundance()    const { return mabundance; }
  /** Mass of the isotope. (In unified atomic mass units u)
      If I understand this correctly this is the same as g/mol. */
  const Numeric&      Mass()         const { return mmass;    }
  /** MYTRAN2 tag numers for all isotopes. -1 means not included. */
  const int&          MytranTag()    const { return mmytrantag;    }
  /** HITRAN-96 tag numers for all isotopes. -1 means not included. */
  const int&          HitranTag()    const { return mhitrantag;    }
  /** JPL tag numbers for all isotopes. Empty array means not included. There
      can be more than one JPL tag for an isotopic species, because in
      JPL different vibrational states have different tags. */
  const ARRAY<int>&   JplTags()      const { return mjpltags;      }

  void SetPartitionFctCoeff( const ARRAY<Numeric>& qcoeff )
  {
    mqcoeff = qcoeff;
    mqcoeff_at_t_ref = -1.;
  }

  Numeric CalculatePartitionFctRatio( Numeric temperature ) const
  {
    Numeric qtemp = CalculatePartitionFctAtTemp( temperature );

    if ( qtemp < 0. ) 
      {
	ostringstream os;
	os << "Partition function of "
	   << "Isotope = " << mname
	   << "is unknown."
	   << endl;
	throw runtime_error(os.str());
      }
    return mqcoeff_at_t_ref / qtemp;
  }

  // calculate the partition function at the reference temperature
  void CalculatePartitionFctAtRefTemp( Numeric temperature )
  {
    //    if (mqcoeff_at_t_ref <= -1.0 )     
    mqcoeff_at_t_ref = CalculatePartitionFctAtTemp( temperature );
  }

private:

  // calculate the partition fct at a certain temperature
  // this is only the prototyping
  Numeric CalculatePartitionFctAtTemp( Numeric temperature ) const;

  string mname;
  Numeric mabundance;
  Numeric mmass;
  int mmytrantag;
  int mhitrantag;
  ARRAY<int> mjpltags;
  ARRAY<Numeric> mqcoeff;
  Numeric mqcoeff_at_t_ref;
};


/** Contains the lookup data for one species.

    \author Stefan Buehler  */
class SpeciesRecord{
public:

  /** Default constructor. */
  SpeciesRecord(){}  ;
  
  /** The constructor used in define_species_data. */
  SpeciesRecord(const char name[],
		const int degfr,
		const ARRAY<IsotopeRecord>& isotope)
    : mname(name),
      mdegfr(degfr),
      misotope(isotope.size())
  {
    // We need to use copy to initialize the ARRAY members. If we use
    // the assignment operator they end up all pointing to the same
    // data!
    copy(isotope,misotope);

#ifndef NDEBUG
      {
	/* Check that the isotopes are correctly sorted. */
	for ( size_t i=0; i<misotope.size()-1; ++i )
	  {
	    assert( misotope[i].Abundance() >= misotope[i+1].Abundance() );
	  }

	/* Check that the Mytran tags are correctly sorted. */
	for ( size_t i=0; i<misotope.size()-1; ++i )
	  {
	    if ( (0<misotope[i].MytranTag()) && (0<misotope[i+1].MytranTag()) )
	      {
		assert( misotope[i].MytranTag() < misotope[i+1].MytranTag() );
	    
		// Also check that the tags have the same base number:
		assert( misotope[i].MytranTag()/10 == misotope[i].MytranTag()/10 );
	      }
	  }

	/* Check that the Hitran tags are correctly sorted. */
	for ( size_t i=0; i<misotope.size()-1; ++i )
	  {
	    if ( (0<misotope[i].HitranTag()) && (0<misotope[i+1].HitranTag()) )
	      {
		assert( misotope[i].HitranTag() < misotope[i+1].HitranTag() );
	    
		// Also check that the tags have the same base number:
		assert( misotope[i].HitranTag()/10 == misotope[i+1].HitranTag()/10 );
	      }
	  }
      }
#endif // #ifndef NDEBUG
  }

  const string&               Name()     const { return mname;     }   
  int                         Degfr()    const { return mdegfr;    }
  const ARRAY<IsotopeRecord>& Isotope()  const { return misotope;  }
  ARRAY<IsotopeRecord>&       Isotope()        { return misotope;  }
  
private:
  /** Species name. */
  string mname;
  /** Degrees of freedom. */
  int mdegfr;
  /** Isotope data. */
  ARRAY<IsotopeRecord> misotope;
};


/** Spectral line catalog data. Here is a description of the ARTS
    catalogue format, largely taken from the Bredbeck book, except
    that some units are slightly changed:

    The line catalogue should not have any fixed column widths because the
    precision of the parameters should not be limited by the format.  The
    catalogue can then be stored as binary or ASCII. In the ASCII version
    the columns are separated by one or more blanks. The line format is
    then specified by only the order and the units of the columns. As
    the catalogue entry for each transition can be quite long, it can be
    broken across lines in the ASCII file. Each new transition is marked
    with a `@' character.

    The first column will contain the species and isotope, following
    the naming scheme described below.  Scientific notation is allowed,
    e.g. 501.12345e9.  The transitions of different isotopes should be
    kept in separate files (this is for the catalogue, the forward model
    should also be able to use a single line file). The suggested line
    format is:

    \verbatim
    Col  Variable                Label    Unit     Comment
    ------------------------------------------------------------------      
    0   `@'                      ENTRY       -     marks start of entry
    1   name                     NAME        -     e.g. O3-666
    2   center frequency            F       Hz     e.g. 501.12345e9 
    3   pressure shift of F       PSF    Hz/Pa    
    4   line intensity             I0   m^2/Hz 
    5   reference temp. for I0   T_I0        K
    6   lower state energy       ELOW     cm-1
    7   air broadened width      AGAM    Hz/Pa     values around 2 GHz/Pa
    8   self broadened width     SGAM    Hz/Pa
    9   AGAM temp. exponent      NAIR        -
    10   SGAM temp. exponent     NSELF       - 
    11   ref. temp. for AGAM, SGAM T_GAM     K
    12   number of aux. parameters N_AUX     -
    13   auxiliary parameter       AUX1      -
    14   ...
    15   error for F                DF      Hz
    16   error for AGAM          DAGAM       %
    17   error for SGAM          DSGAM       %
    18   error for NAIR          DNAIR       %
    19   error for NSELF        DNSELF       %
    20   error for PSF            DPSF       %
    21   quantum number code                       string or number ?
    22   lower state quanta                        string inside quotes
    23   upper state quanta                        string inside quotes
    24   information source of F
    25   information source of I0
    26   information source of line width variables
    27   information source of overlap constants
    \endverbatim

    One line could be:
    {\small
    \verbatim
    @ O3-666 110.83604e9 0 12.43453e-22 300 17.5973 1.52 2.03 0.73 0.73 296 ....   
    \endverbatim}

    The format used in the line file used by the FM
    can be a truncated version of the full line format.

    Some species need special parameters that are not needed by other
    species (for example overlap coefficients for O2). In the case of
    oxygen two parameters are sufficient to describe the overlap, but
    other species, e.g., methane, may need more coefficients. The
    default for \texttt{N\_AUX} is zero. In that case, no further
    \texttt{AUX} fields are present.

    The names of the private members and public access functions of
    this data structure follow the above table. The only difference is
    that underscores are omited and only the first letter of each name
    is capitalized. This is for consistency with the notation
    elsewhere in the program.

    \author Stefan Buehler */
class LineRecord {
public:

  /** Default constructor. Initialize to default values. The indices
      are initialized to large numbers, so that we at least get range
      errors when we try to used un-initialized data. */
  LineRecord()
    : mspecies (1000000),
      misotope (1000000),
      mf       (0.     ),
      mpsf     (0.     ),
      mi0      (0.     ),
      mti0     (0.     ),
      melow    (0.     ),
      magam    (0.     ),
      msgam    (0.     ),
      mnair    (0.     ),
      mnself   (0.     ),
      mtgam    (0.     ),
      maux     (       )
 { /* Nothing to do here. */ }

  /** Constructor that sets all data elements explicitly. If
      assertions are not disabled (i.e., if NDEBUG is not #defined),
      assert statements check that the species and isotope data
      exists. */
  LineRecord( size_t  	     	    species,
	      size_t  	     	    isotope,
	      Numeric 	     	    f,
	      Numeric 	     	    psf,
	      Numeric 	     	    i0,
	      Numeric 	     	    ti0,
	      Numeric 	     	    elow,
	      Numeric 	     	    agam,
	      Numeric 	     	    sgam,
	      Numeric 	     	    nair,
	      Numeric 	     	    nself,
	      Numeric 	     	    tgam,
	      const ARRAY<Numeric>& aux       )
    : mspecies (species	   ),
      misotope (isotope	   ),
      mf       (f      	   ),
      mpsf     (psf    	   ),
      mi0      (i0     	   ),
      mti0     (ti0    	   ),
      melow    (elow   	   ),
      magam    (agam   	   ),
      msgam    (sgam   	   ),
      mnair    (nair   	   ),
      mnself   (nself  	   ),
      mtgam    (tgam   	   ),  
      maux     (aux.size() )
  {
    // We need to use copy to initialize the ARRAY members. If we use
    // the assignment operator they end up all pointing to the same
    // data!
    copy(aux,maux);

    // Check if this species is legal, i.e., if species and isotope
    // data exists.
    extern ARRAY<SpeciesRecord> species_data;
    assert( mspecies < species_data.size() );
    assert( misotope < species_data[mspecies].Isotope().size() );
    // if ever this constructor is used, here is the calculation of
    // the partition fct at the reference temperature
    species_data[mspecies].Isotope()[misotope].CalculatePartitionFctAtRefTemp( mti0 ) ;
  }

  /** The index of the molecular species that this line belongs
      to. The species data can be accessed by species_data[Species()]. */
  size_t Species() const { return mspecies; }

  /** The index of the isotopic species that this line belongs
      to. The isotopic species data can be accessed by
      species_data[Species()].Isotope()[Isotope()].  */
  size_t Isotope() const { return misotope; }

  /** The full name of the species and isotope. E.g., `O3-666'. The
      name is found by looking up the information in species_data,
      using the species and isotope index. */
  string Name() const {
    // The species lookup data:
    extern const ARRAY<SpeciesRecord> species_data;
    const SpeciesRecord& sr = species_data[mspecies];
    return sr.Name() + "-" + sr.Isotope()[misotope].Name();
  }

  /** The matching SpeciesRecord from species_data. To get at the
      species data of a LineRecord lr, you can use:
      \begin{enumerate}
      \item species_data[lr.Species()]
      \item lr.SpeciesData()
      \end{enumerate}
      The only advantages of the latter are that the notation is
      slightly nicer and that you don't have to declare the external
      variable species_data. */
  const SpeciesRecord& SpeciesData() const {
    // The species lookup data:
    extern const ARRAY<SpeciesRecord> species_data;
    return species_data[mspecies];
  }

  /** The matching IsotopeRecord from species_data. The IsotopeRecord
      is a subset of the SpeciesRecord. To get at the isotope data of
      a LineRecord lr, you can use:
      \begin{enumerate}
      \item species_data[lr.Species()].Isotope()[lr.Isotope()]
      \item lr.SpeciesData().Isotope()[lr.Isotope()]
      \item lr.IsotopeData()
      \end{enumerate}
      The last option is clearly the shortest, and has the advantage
      that you don't have to declare the external variable
      species_data. */
  const IsotopeRecord& IsotopeData() const {
    // The species lookup data:
    extern const ARRAY<SpeciesRecord> species_data;
    return species_data[mspecies].Isotope()[misotope];
  }

  /** The line center frequency in <b> Hz</b>. */
  Numeric F() const     { return mf; }

  /** Set the line center frequency in <b> Hz</b>. */
  void setF( Numeric new_mf ) { mf = new_mf; }

  /** The pressure shift parameter in <b> Hz/Pa</b>. */
  Numeric Psf() const   { return mpsf; }

  /** The line intensity in <b> m^2*Hz</b> at the reference temperature \c Ti0. 

    The line intensity \f$I_0\f$ is defined by:
    
    \f[
    \alpha(\nu) = n \, x \, I_0(T) \, F(\nu)
    \f]

    where \f$\alpha\f$ is the absorption coefficient (in <b>
    m^-1</b>), \f$\nu\f$ is frequency, \f$n\f$ is the
    total number density, \f$x\f$ is the volume mixing ratio, and
    \f$F(\nu)\f$ is the lineshape function. */
  Numeric I0() const    { return mi0; }

  /** Reference temperature for I0 in <b> K</b>: */
  Numeric Ti0() const   { return mti0; }

  /** Lower state energy in <b> cm^-1</b>: */
  Numeric Elow() const  { return melow; }

  /** Air broadened width in <b> Hz/Pa</b>: */
  Numeric Agam() const  { return magam; }

  /** Self broadened width in <b> Hz/Pa</b>: */
  Numeric Sgam() const  { return msgam; }

  /** AGAM temperature exponent (dimensionless): */
  Numeric Nair() const  { return mnair; }

  /** SGAM temperature exponent (dimensionless): */
  Numeric Nself() const { return mnself; }

  /** Reference temperature for AGAM and SGAM in <b> K</b>: */
  Numeric Tgam() const  { return mtgam; }

  /** Number of auxiliary parameters. This function is actually
      redundant, since the number of auxiliary parameters can also be
      obtained directly with Aux.size(). I just added the function in
      order to have consistency of the interface with the catalgue
      format. */
  size_t Naux() const   { return maux.size(); }

  /** Auxiliary parameters. */
  const ARRAY<Numeric>& Aux() const { return maux; }

  /** Read one line from a stream associated with a HITRAN file. The HITRAN
    format is as follows (directly from the HITRAN documentation):

    \verbatim
    Each line consists of 100
    bytes of ASCII text data, followed by a line feed (ASCII 10) and
    carriage return (ASCII 13) character, for a total of 102 bytes per line.
    Each line can be read using the following READ and FORMAT statement pair
    (for a FORTRAN sequential access read):

          READ(3,800) MO,ISO,V,S,R,AGAM,SGAM,E,N,d,V1,V2,Q1,Q2,IERF,IERS,
         *  IERH,IREFF,IREFS,IREFH
    800   FORMAT(I2,I1,F12.6,1P2E10.3,0P2F5.4,F10.4,F4.2,F8.6,2I3,2A9,3I1,3I2)

    Each item is defined below, with its format shown in parenthesis.

      MO  (I2)  = molecule number
      ISO (I1)  = isotope number (1 = most abundant, 2 = second, etc)
      V (F12.6) = frequency of transition in wavenumbers (cm-1)
      S (E10.3) = intensity in cm-1/(molec * cm-2) at 296 Kelvin
      R (E10.3) = transition probability squared in Debyes**2
      AGAM (F5.4) = air-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
      SGAM (F5.4) = self-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
      E (F10.4) = lower state energy in wavenumbers (cm-1)
      N (F4.2) = coefficient of temperature dependence of air-broadened halfwidth
      d (F8.6) = shift of transition due to pressure (cm-1)
      V1 (I3) = upper state global quanta index
      V2 (I3) = lower state global quanta index
      Q1 (A9) = upper state local quanta
      Q2 (A9) = lower state local quanta
      IERF (I1) = accuracy index for frequency reference
      IERS (I1) = accuracy index for intensity reference
      IERH (I1) = accuracy index for halfwidth reference
      IREFF (I2) = lookup index for frequency
      IREFS (I2) = lookup index for intensity
      IREFH (I2) = lookup index for halfwidth

    The molecule numbers are encoded as shown in the table below:

      0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=   CO
      6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=  NH3
     12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=   HI
     18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=  HCN
     24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29= COF2
     30=  SF6   31=  H2S   32=HCOOH
    \endverbatim

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped.

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \author Stefan Buehler */
  bool ReadFromHitranStream(istream& is);



  /** Read one line from a stream associated with a MYTRAN2 file. The MYTRAN2
    format is as follows (directly taken from the abs_my.c documentation):

    \verbatim
    The MYTRAN format is as follows (FORTRAN notation):
    FORMAT(I2,I1,F13.4,1PE10.3,0P2F5.2,F10.4,2F4.2,F8.6,F6.4,2I3,2A9,4I1,3I2)
   
    Each item is defined below, with its FORMAT string shown in
    parenthesis.
   
       MO  (I2)      = molecule number
       ISO (I1)      = isotope number (1 = most abundant, 2 = second, etc)
    *  F (F13.4)     = frequency of transition in MHz
    *  errf (F8.4)   = error in f in MHz
       S (E10.3)     = intensity in cm-1/(molec * cm-2) at 296 K
    *  AGAM (F5.4)   = air-broadened halfwidth (HWHM) in MHz/Torr at Tref
    *  SGAM (F5.4)   = self-broadened halfwidth (HWHM) in MHz/Torr at Tref
       E (F10.4)     = lower state energy in wavenumbers (cm-1)
       N (F4.2)      = coefficient of temperature dependence of 
	 	       air-broadened halfwidth
    *  N_self (F4.2) = coefficient of temperature dependence of 
		       self-broadened halfwidth
    *  Tref (F7.2)   = reference temperature for AGAM and SGAM 
    *  d (F8.6)      = shift of transition due to pressure (MHz/Torr)
       V1 (I3) 	     = upper state global quanta index
       V2 (I3) 	     = lower state global quanta index
       Q1 (A9) 	     = upper state local quanta
       Q2 (A9) 	     = lower state local quanta
       IERS (I1)     = accuracy index for S
       IERH (I1)     = accuracy index for AGAM
    *  IERN (I1)     = accuracy index for N
       IREFF (I2)    = lookup index for F
       IREFS (I2)    = lookup index for S
       IREFH (I2)    = lookup index for AGAM
   
    The asterisks mark entries that are different from HITRAN.

    Note that AGAM and SGAM are for the temperature Tref, while S is
    still for 296 K!
   
    The molecule numbers are encoded as shown in the table below:
   
     0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=   CO
     6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=  NH3
    12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=   HI
    18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=  HCN
    24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29= COF2
    30=  SF6   31=  H2S   32=HCOOH   33= HO2    34=    O   35= CLONO2
    36=  NO+   37= Null   38= Null   39= Null   40=H2O_L   41= Null
    42= Null   43= OCLO   44= Null   45= Null   46=BRO     47= Null
    48= H2SO4  49=CL2O2

    All molecule numbers are from HITRAN, except for species with id's
    greater or equals 40, which are not included in HITRAN.
    (E.g.: For BrO, iso=1 is Br-79-O,iso=2 is  Br-81-O.)
    \endverbatim

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped.

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \date 31.10.00
    \author Axel von Engeln 
  */
  bool ReadFromMytran2Stream(istream& is);


  /** Read one line from a stream associated with a JPL file. The JPL
    format is as follows (directly taken from the jpl documentation):

    \verbatim 
    The catalog line files are composed of 80-character lines, with one
    line entry per spectral line.  The format of each line is:

    \label{lfmt}
    \begin{tabular}{@{}lccccccccr@{}}
    FREQ, & ERR, & LGINT, & DR, & ELO, & GUP, & TAG, & QNFMT, & QN${'}$, & QN${''}$\\ 
    (F13.4, & F8.4, & F8.4, & I2, & F10.4, & I3, & I7, & I4, & 6I2, & 6I2)\\
    \end{tabular}

    \begin{tabular}{lp{4.5in}} 
    FREQ: & Frequency of the line in MHz.\\ 
    ERR: & Estimated or experimental error of FREQ in MHz.\\ 
    LGINT: &Base 10 logarithm of the integrated intensity 
    in units of \linebreak nm$^2$$\cdot$MHz at 300 K. (See Section 3 for 
    conversions to other units.)\\ 
    DR: & Degrees of freedom in the rotational partition 
    function (0 for atoms, 2 for linear molecules, and 3 for nonlinear 
    molecules).\\ 
    ELO: &Lower state energy in cm$^{-1}$ relative to the lowest energy 
    spin--rotation level in ground vibronic state.\\ 
    GUP: & Upper state degeneracy.\\ 
    TAG: & Species tag or molecular identifier. 
    A negative value flags that the line frequency has 
    been measured in the laboratory.  The absolute value of TAG is then the 
    species tag and ERR is the reported experimental error.  The three most 
    significant digits of the species tag are coded as the mass number of the 
    species, as explained above.\\ 
    QNFMT: &Identifies the format of the quantum numbers 
    given in the field QN. These quantum number formats are given in Section 5 
    and are different from those in the first two editions of the catalog.\\ 
    QN${'}$: & Quantum numbers for the upper state coded 
    according to QNFMT.\\ 
    QN${''}$: & Quantum numbers for the lower state.\\
    \end{tabular} 
    \endverbatim

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped (unused in jpl).

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \date 01.11.00
    \author Axel von Engeln */
  bool ReadFromJplStream(istream& is);

  /** Read one line from a stream associated with an Arts file.

      Format: see arts distribution

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped.

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \date 15.12.00
    \author Oliver Lemke*/
  bool ReadFromArtsStream(istream& is);


private:
  // Molecular species index: 
  size_t mspecies;
  // Isotopic species index:
  size_t misotope;
  // The line center frequency in Hz:
  Numeric mf;
  // The pressure shift parameter in Hz/Pa:
  Numeric mpsf;
  // The line intensity in m^2/Hz:
  Numeric mi0;
  // Reference temperature for I0 in K:
  Numeric mti0;
  // Lower state energy in cm^-1:
  Numeric melow;
  // Air broadened width in Hz/Pa:
  Numeric magam;
  // Self broadened width in Hz/Pa:
  Numeric msgam;
  // AGAM temperature exponent (dimensionless):
  Numeric mnair;
  // SGAM temperature exponent (dimensionless):
  Numeric mnself;
  // Reference temperature for AGAM and SGAM in K:
  Numeric mtgam;
  // Array to hold auxiliary parameters:
  ARRAY<Numeric> maux;
};

// is needed to map jpl tags/arts identifier to the species/isotope data within arts
class SpecIsoMap{
public:
  SpecIsoMap():mspeciesindex(0), misotopeindex(0){}
  SpecIsoMap(const size_t& speciesindex,
		const size_t& isotopeindex)
    : mspeciesindex(speciesindex),
      misotopeindex(isotopeindex) 
  {}

  // Return the index to the species 
  const int& Speciesindex() const { return mspeciesindex; }
  // Return the index to the isotope
  const int& Isotopeindex() const { return misotopeindex; }

private:
  int mspeciesindex;
  int misotopeindex;
};



/** Holds a list of spectral line data.
    \author Stefan Buehler */
typedef ARRAY<LineRecord> ARRAYofLineRecord;

/** Holds a lists of spectral line data for each tag group.
    Dimensions: (tag_groups.size()) (number of lines for this tag)
    \author Stefan Buehler */
typedef ARRAY< ARRAY<LineRecord> > ARRAYofARRAYofLineRecord;



/** Output operator for LineRecord. The result should look like a
    catalogue line.

    \author Stefan Buehler */
ostream& operator << (ostream& os, const LineRecord& lr);



/** Define the species data map.

    \author Stefan Buehler  */
void define_species_map();



//------------------------------< Tag Group Stuff >------------------------------

/** A tag group can consist of the sum of several of these.

    \author Stefan Buehler */
class OneTag {
public:
  /** Default constructor. */
  OneTag() { /* Nothing to be done here. */ }

  /** Constructor from a tag definition string (Bredbeck
      notation). For examples see member function Name(). 

      \exception runtime_error The given string could not be mapped to
      a sensible tag description. */
  OneTag(string def); 

  /** Return the full name of this tag according to Bredbeck
      convention. Examples:
      \verbatim
      O3-*-*-*         : All O3 lines
      O3-666-*-*       : All O3-666 lines
      O3-*-500e9-501e9 : All O3 lines between 500 and 501 GHz.
      \endverbatim */
  string Name() const;
    
  /** Molecular species index. */
  size_t Species() const { return mspecies; }

  /** Isotopic species index.
      If this is equal to the number of isotopes (one more than
      allowed) it means all isotopes of this species. */ 
  size_t Isotope() const { return misotope; }

  /** The lower line center frequency in Hz.
      If this is <0 it means no lower limit. */
  Numeric Lf() const { return mlf; }

  /** The upper line center frequency in Hz:
      If this is <0 it means no upper limit. */
  Numeric Uf() const { return muf; }

private:
  // Molecular species index: 
  size_t mspecies;
  // Isotopic species index.
  // If this is equal to the number of isotopes (one more than
  // allowed) it means all isotopes of this species.
  size_t misotope;
  // The lower line center frequency in Hz.
  // If this is <0 it means no lower limit. 
  Numeric mlf;
  // The upper line center frequency in Hz:
  // If this is <0 it means no upper limit. 
  Numeric muf;
};


/** Output operator for OneTag. 

    \author Stefan Buehler */
ostream& operator << (ostream& os, const OneTag& ot);


/** Contains the available tag groups. Contrary to the Bredbeck
    definition, tag groups may only consist of tags belonging to the
    same species. The reason for this is that there is one VMR profile
    associated with each tag group.

    \author Stefan Buehler */
typedef  ARRAY< ARRAY<OneTag> > TagGroups;


void get_tagindex_for_strings( 
              ARRAYofsizet&   tags1_index, 
        const TagGroups&      tags1, 
        const ARRAYofstring&  tags2_strings );

void get_tag_group_index_for_tag_group( 
              size_t&         tags1_index, 
        const TagGroups&      tags1, 
        const ARRAY<OneTag>&  tags2 );

string get_tag_group_name( const ARRAY<OneTag>& tg );

// Doc header in absorption.cc
void write_lines_to_stream(ostream& os,
			   const ARRAYofLineRecord& lines);


/** Calculate line absorption cross sections for one tag group. All
    lines in the line list must belong to the same species. This must
    be ensured by lines_per_tgCreateFromLines, so it is only verified
    with assert. Also, the input vectors p_abs, and t_abs must all
    have the same dimension.

    This is mainly a copy of abs_species which is removed now, with
    the difference that the vmrs are removed from the absorption
    coefficient calculation. (the vmr is still used for the self
    broadening)

    Continua are not handled by this function, you have to call
    xsec_continuum_tag for those.

    \retval xsec   Cross section of one tag group.
    \param f_mono  Frequency grid.
    \param p_abs   Pressure grid.
    \param t_abs   Temperatures associated with p_abs.
    \param h2o_abs Total volume mixing ratio of water vapor.
    \param vmr     Volume mixing ratio of the calculated species.
    \param lines   The spectroscopic line list.
    \param ind_ls  Lineshape specifications.

    \author Stefan Buehler and Axel von Engeln
    \date   2001-01-11 */
void xsec_species( MATRIX&                 xsec,
		   const VECTOR&  	   f_mono,
		   const VECTOR&  	   p_abs,
		   const VECTOR&  	   t_abs,           
		   const VECTOR&  	   h2o_abs,           
		   const VECTOR&            vmr,
		   const ARRAYofLineRecord& lines,
		   const size_t             ind_ls,
		   const size_t             ind_lsn,
		   const Numeric            cutoff);

#endif // absorption_h
