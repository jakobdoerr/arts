/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
                      Patrick Eriksson <patrick@rss.chalmers.se>

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

/**
  \file  arts.h

  The global header file for ARTS. This file is included directly or
  indirectly by each and every ARTS source file. It must therefor not
  contain stuff that should not always be present.

  Note that you do not have to include this file explicitly in many
  cases, since it is included directly or indirectly by most ARTS
  header files.

  \author Stefan Buehler
  \date 16.05.1999 */

/** \mainpage
    
    This automatic documentation is still \e very experimental. We
    have just switched from doc++ to doxygen. You can use the HTML
    version to browse the source text. Just point and click, and
    eventually you will see the real implementation of functions and
    classes. 

    If you are looking for a more comprehensive text, check out the
    Arts User Guide that is also distributed along with the
    program. Section `Documentation' in Chapter `The art of developing
    ARTS' there also tells you how you should add documentation
    headers to your code if you are an ARTS developer.
 */

#ifndef arts_h
#define arts_h


//----------< First of all, include the configuration header >----------
// This header is generated by the configure script.
#if HAVE_CONFIG_H
#include "config.h"
#endif		
// FIXME: Put an else clause here that stops the compilation with an
// error message.

//----------< Standard library header files: >----------

// Decided that not all standard headers should be included
// everywhere. Now, each header file of ARTS must include the standad
// header files that it needs. Only the most basic ones are included
// here, notably our local sstream implementation. Also cassert.

// #include <iostream>		// Standard stream library
// #include <iomanip>
// #include <fstream>
// #include <string>		// Standard string library
// #include <vector>
// #include <map>
// #include <stdarg.h>
// #include <math.h>
// #include <cfloat>
// #include <typeinfo>
// #include <stdexcept>
// //#include <algorithm>
// #include <ctype.h>
// #include <climits>

// String stream library. This is included with the ARTS source code
// for now, since it is still missing in the current version of EGCS
// (egcs-2.91.66). Should be removed when stringstreams work as they
// should in the standard EGCS distribution.
#include "sstream.h"


//--------------------< Set floating point type >--------------------
/** The type to use for all floating point numbers. You should never
    use float or double explicitly, unless you have a very good
    reason. Always use this type instead.  */
typedef double Numeric;

//--------------------< Set integer type >--------------------
/** The type to use for all integer numbers and indices. You should never
    use int, long, or size_t explicitly, unless you have a very good
    reason. Always use this type instead.  */
typedef long Index;

//--------------------< Set string type >--------------------
/** The type to use for all strings. This is just to have consistent
    notation for all the atomic ARTS types. */ 
//typedef string String;

// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Index. */
typedef Array<Index> ArrayOfIndex;

/** An array of Numeric. */
typedef Array<Numeric> ArrayOfNumeric;




//--------------------< Set NDEBUG? >--------------------
/* Define this in order to turn off all debuggin features (TNT range
    checking, assertions, ...) */
#undef NDEBUG
//#define NDEBUG 1

// C Assert macro:
// Could be moved to config.h in the future. This must be included
// *after* the definition of NDEBUG, since NDEBUG changes the
// behaviour of assert (assert does nothing if NDEBUG is set).
#include <cassert>


//-----------< define a quick output for std::vector<> >----------
/* A quick hack output for std::vector<>. This is only for
    debugging purposes.
    \return As for all output operator functions, the output stream is 
            returned.
    \param  os Output stream
    \param  v  Vector to print                      
    \author SAB  */  
// template<class T>
// ostream& operator<<(ostream& os, const std::vector<T>& v)
// {
//   for (std::vector<T>::const_iterator i=v.begin(); i<v.end()-1; ++i)
//     os << *i << ", ";
//   os << *(v.end()-1);
//   return os;
// }



//---------------< Global variable declarations >---------------
// (Definitions of these are in FIXME: where?)


//---------------< Global function declarations: >---------------
// Documentations are with function definitions.
void define_wsv_group_names();	
void define_species_data();
void define_lineshape_data();
void define_lineshape_norm_data();

//
// Physical constants are now in constants.cc
//



#endif // arts_h



