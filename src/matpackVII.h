/* Copyright (C) 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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
   Implementation of Tensors of Rank 7.

   Dimensions are called: library, vitrine, shelf, book, page, row, column.
   or short:              l,       v,       s,     b,    p,    r,   c
  
   \author Stefan Buehler
   \date   2001-11-22
*/

#ifndef matpackVII_h
#define matpackVII_h

#include <iomanip>
#include "matpackVI.h"

/** The outermost iterator class for rank 7 tensors. This takes into
    account the defined strided. */
class Iterator7D {
public:
  // Constructors:
  Iterator7D();
  Iterator7D(const Iterator7D& o);
  Iterator7D(const Tensor6View& x, Index stride);

  // Operators:
  Iterator7D& operator++();
  bool operator!=(const Iterator7D& other) const;
  Tensor6View* const operator->();
  Tensor6View& operator*();
  
private:
  /** Current position. */
  Tensor6View msv;
  /** Stride. */
  Index mstride;
};

/** Const version of Iterator7D. */
class ConstIterator7D {
public:
  // Constructors:
  ConstIterator7D();
  ConstIterator7D(const ConstIterator7D& o);
  ConstIterator7D(const ConstTensor6View& x, Index stride);

  // Operators:
  ConstIterator7D& operator++();
  bool operator!=(const ConstIterator7D& other) const;
  const ConstTensor6View* operator->() const;
  const ConstTensor6View& operator*() const;

private:
  /** Current position. */
  ConstTensor6View msv;
  /** Stride. */
  Index mstride;
};


// Declare class Tensor7:
class Tensor7;


/** A constant view of a Tensor7.

This, together with the derived class Tensor7View, contains the
main implementation of a Tensor7. It defines the concepts of
Tensor7View. Plus additionally the recursive subrange operator,
which makes it possible to create a Tensor7View from a subrange of
a Tensor7View.

Dimensions are called: library, vitrine, shelf, book, page, row, column.
or short:              l,       v,       s,     b,    p,    r,   c

The class Tensor7 is just a special case of a Tensor7View
which also allocates storage. */
class ConstTensor7View {
public:
  // Member functions:
  Index nlibraries() const;
  Index nvitrines()  const;
  Index nshelves()   const;
  Index nbooks()     const;
  Index npages()     const;
  Index nrows()      const;
  Index ncols()      const;

  // Const index operators:

  // Result 7D (1 combination)
  // -------
  ConstTensor7View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 6D (7 combinations)
  // ------|
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // -----|-
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ----|--
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ---|---
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // --|----
  ConstTensor6View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|-----
  ConstTensor6View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |------
  ConstTensor6View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 5D (6+5+4+3+2+1 = 21 combinations)
  // -----||
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ----|-|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // ---|--|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // --|---|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // -|----|
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |-----|
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ----||-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // ---|-|-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // --|--|-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // -|---|-
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |----|-
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ---||--
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // --|-|--
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -|--|--
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |---|--
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // --||---
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|-|---
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |--|---
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -||----
  ConstTensor5View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-|----
  ConstTensor5View operator()( Index l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||-----
  ConstTensor5View operator()( Index l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 4D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ----|||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ---|-||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // --|--||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // -|---||
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |----||
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ---||-|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // --|-|-|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -|--|-|
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |---|-|
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // --||--|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|-|--|
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |--|--|
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -||---|
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |-|---|
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ||----|
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ---|||-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // --|-||-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -|--||-
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |---||-
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // --||-|-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|-|-|-
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |--|-|-
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -||--|-
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |-|--|-
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ||---|-
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // --|||--
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -|-||--
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |--||--
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -||-|--
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |-|-|--
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ||--|--
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -|||---
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-||---
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||-|---
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |||----
  ConstTensor4View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 3D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ||||---
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |||-|--
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ||-||--
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |-|||--
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -||||--
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |||--|-
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ||-|-|-
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |-||-|-
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|||-|-
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // ||--||-
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |-|-||-
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -||-||-
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |--|||-
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|-|||-
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // --||||-
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |||---|
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ||-|--|
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |-||--|
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|||--|
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // ||--|-|
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |-|-|-|
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -||-|-|
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |--||-|
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -|-||-|
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // --|||-|
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // ||---||
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |-|--||
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // -||--||
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |--|-||
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|-|-||
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // --||-||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |---|||
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -|--|||
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // --|-|||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ---||||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 2D (6+5+4+3+2+1 = 21 combinations)
  // |||||--
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // ||||-|-
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |||-||-
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // ||-|||-
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |-||||-
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|||||-
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // ||||--|
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |||-|-|
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // ||-||-|
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // |-|||-|
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -||||-|
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // |||--||
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ||-|-||
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |-||-||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|||-||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // ||--|||
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // |-|-|||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -||-|||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // |--||||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -|-||||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // --|||||
  ConstMatrixView  operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 1D (7 combinations)
  // ||||||-
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |||||-|
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // ||||-||
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |||-|||
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ||-||||
  ConstVectorView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // |-|||||
  ConstVectorView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -||||||
  ConstVectorView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result scalar (1 combination)
  // |||||||
  Numeric          operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;


  // Functions returning iterators:
  ConstIterator7D begin() const;
  ConstIterator7D end() const;
  
  // Friends:
  friend class Tensor7View;

  // Special constructor to make a Tensor7 view of a Tensor6.
  ConstTensor7View(const ConstTensor6View& a);

protected:
  // Constructors:
  ConstTensor7View();
  ConstTensor7View(Numeric *data,
		   const Range& l,
		   const Range& v, const Range& s, const Range& b,
		   const Range& p, const Range& r, const Range& c);
  ConstTensor7View(Numeric *data,
		   const Range& pl,
		   const Range& pv, const Range& ps, const Range& pb,
		   const Range& pp, const Range& pr, const Range& pc,
		   const Range& nl,
		   const Range& nv, const Range& ns, const Range& nb,
		   const Range& np, const Range& nr, const Range& nc);

  // Data members:
  // -------------
  /** The library range of mdata that is actually used. */
  Range mlr;
  /** The vitrine range of mdata that is actually used. */
  Range mvr;
  /** The shelf range of mdata that is actually used. */
  Range msr;
  /** The book range of mdata that is actually used. */
  Range mbr;
  /** The page range of mdata that is actually used. */
  Range mpr;
  /** The row range of mdata that is actually used. */
  Range mrr;
  /** The column range of mdata that is actually used. */
  Range mcr;
  /** Pointer to the plain C array that holds the data */
  Numeric *mdata;
};

/** The Tensor7View class

This contains the main implementation of a Tensor7. It defines
the concepts of Tensor7View. Plus additionally the recursive
subrange operator, which makes it possible to create a Tensor7View
from a subrange of a Tensor7View. 

The class Tensor7 is just a special case of a Tensor7View
which also allocates storage. */
class Tensor7View : public ConstTensor7View {
public:

  // Const index operators:

  // Result 7D (1 combination)
  // -------
  ConstTensor7View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 6D (7 combinations)
  // ------|
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // -----|-
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ----|--
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ---|---
  ConstTensor6View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // --|----
  ConstTensor6View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|-----
  ConstTensor6View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |------
  ConstTensor6View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 5D (6+5+4+3+2+1 = 21 combinations)
  // -----||
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ----|-|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // ---|--|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // --|---|
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // -|----|
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |-----|
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ----||-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // ---|-|-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // --|--|-
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // -|---|-
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |----|-
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ---||--
  ConstTensor5View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // --|-|--
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -|--|--
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |---|--
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // --||---
  ConstTensor5View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -|-|---
  ConstTensor5View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |--|---
  ConstTensor5View operator()( Index l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // -||----
  ConstTensor5View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-|----
  ConstTensor5View operator()( Index l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||-----
  ConstTensor5View operator()( Index l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 4D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ----|||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ---|-||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // --|--||
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // -|---||
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |----||
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ---||-|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // --|-|-|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -|--|-|
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |---|-|
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // --||--|
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|-|--|
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |--|--|
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -||---|
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // |-|---|
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ||----|
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ---|||-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // --|-||-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -|--||-
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |---||-
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // --||-|-
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|-|-|-
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |--|-|-
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -||--|-
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // |-|--|-
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ||---|-
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // --|||--
  ConstTensor4View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -|-||--
  ConstTensor4View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |--||--
  ConstTensor4View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -||-|--
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // |-|-|--
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ||--|--
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // -|||---
  ConstTensor4View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |-||---
  ConstTensor4View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // ||-|---
  ConstTensor4View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |||----
  ConstTensor4View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, const Range& c) const;

  // Result 3D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ||||---
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, const Range& c) const;
  // |||-|--
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, const Range& c) const;
  // ||-||--
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |-|||--
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // -||||--
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // |||--|-
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, const Range& c) const;
  // ||-|-|-
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |-||-|-
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // -|||-|-
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // ||--||-
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |-|-||-
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // -||-||-
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // |--|||-
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|-|||-
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // --||||-
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |||---|
  ConstTensor3View operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, const Range& r, Index        c) const;
  // ||-|--|
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |-||--|
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // -|||--|
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // ||--|-|
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |-|-|-|
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // -||-|-|
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // |--||-|
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -|-||-|
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // --|||-|
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // ||---||
  ConstTensor3View operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |-|--||
  ConstTensor3View operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // -||--||
  ConstTensor3View operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // |--|-||
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|-|-||
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // --||-||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |---|||
  ConstTensor3View operator()( Index        l,
			       const Range& v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -|--|||
  ConstTensor3View operator()( const Range& l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // --|-|||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ---||||
  ConstTensor3View operator()( const Range& l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 2D (6+5+4+3+2+1 = 21 combinations)
  // |||||--
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, const Range& c) const;
  // ||||-|-
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, const Range& c) const;
  // |||-||-
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, const Range& c) const;
  // ||-|||-
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |-||||-
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // -|||||-
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // ||||--|
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, const Range& r, Index        c) const;
  // |||-|-|
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, const Range& r, Index        c) const;
  // ||-||-|
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // |-|||-|
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // -||||-|
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // |||--||
  ConstMatrixView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       const Range& p, Index        r, Index        c) const;
  // ||-|-||
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |-||-||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // -|||-||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // ||--|||
  ConstMatrixView  operator()( Index        l,
			       Index        v, const Range& s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // |-|-|||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // -||-|||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // |--||||
  ConstMatrixView  operator()( Index        l,
			       const Range& v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -|-||||
  ConstMatrixView  operator()( const Range& l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // --|||||
  ConstMatrixView  operator()( const Range& l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result 1D (7 combinations)
  // ||||||-
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, const Range& c) const;
  // |||||-|
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, const Range& r, Index        c) const;
  // ||||-||
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       const Range& p, Index        r, Index        c) const;
  // |||-|||
  ConstVectorView  operator()( Index        l,
			       Index        v, Index        s, const Range& b,
			       Index        p, Index        r, Index        c) const;
  // ||-||||
  ConstVectorView  operator()( Index        l,
			       Index        v, const Range& s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // |-|||||
  ConstVectorView  operator()( Index        l,
			       const Range& v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;
  // -||||||
  ConstVectorView  operator()( const Range& l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;

  // Result scalar (1 combination)
  // |||||||
  Numeric          operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c) const;


  // Non-const index operators:

  // Result 7D (1 combination)
  // -------
  Tensor7View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 6D (7 combinations)
  // ------|
  Tensor6View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // -----|-
  Tensor6View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // ----|--
  Tensor6View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // ---|---
  Tensor6View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // --|----
  Tensor6View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);
  // -|-----
  Tensor6View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);
  // |------
  Tensor6View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 5D (6+5+4+3+2+1 = 21 combinations)
  // -----||
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // ----|-|
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // ---|--|
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // --|---|
  Tensor5View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // -|----|
  Tensor5View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // |-----|
  Tensor5View operator()( Index l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ----||-
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // ---|-|-
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // --|--|-
  Tensor5View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // -|---|-
  Tensor5View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // |----|-
  Tensor5View operator()( Index l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // ---||--
  Tensor5View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // --|-|--
  Tensor5View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // -|--|--
  Tensor5View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // |---|--
  Tensor5View operator()( Index l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // --||---
  Tensor5View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // -|-|---
  Tensor5View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // |--|---
  Tensor5View operator()( Index l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // -||----
  Tensor5View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);
  // |-|----
  Tensor5View operator()( Index l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);
  // ||-----
  Tensor5View operator()( Index l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 4D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ----|||
  Tensor4View operator()( const Range& l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // ---|-||
  Tensor4View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // --|--||
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // -|---||
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // |----||
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // ---||-|
  Tensor4View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // --|-|-|
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // -|--|-|
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // |---|-|
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // --||--|
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // -|-|--|
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // |--|--|
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // -||---|
  Tensor4View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // |-|---|
  Tensor4View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ||----|
  Tensor4View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ---|||-
  Tensor4View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // --|-||-
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // -|--||-
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // |---||-
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // --||-|-
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // -|-|-|-
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // |--|-|-
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // -||--|-
  Tensor4View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // |-|--|-
  Tensor4View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // ||---|-
  Tensor4View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // --|||--
  Tensor4View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // -|-||--
  Tensor4View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // |--||--
  Tensor4View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // -||-|--
  Tensor4View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // |-|-|--
  Tensor4View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // ||--|--
  Tensor4View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // -|||---
  Tensor4View operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // |-||---
  Tensor4View operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // ||-|---
  Tensor4View operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // |||----
  Tensor4View operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, const Range& c);

  // Result 3D (5+4+3+2+1 +4+3+2+1 +3+2+1 +2+1 +1 = 35 combinations)
  // ||||---
  Tensor3View operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, const Range& c);
  // |||-|--
  Tensor3View operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, const Range& c);
  // ||-||--
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // |-|||--
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // -||||--
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // |||--|-
  Tensor3View operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, const Range& c);
  // ||-|-|-
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // |-||-|-
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // -|||-|-
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // ||--||-
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // |-|-||-
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // -||-||-
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // |--|||-
  Tensor3View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // -|-|||-
  Tensor3View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // --||||-
  Tensor3View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // |||---|
  Tensor3View operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, const Range& r, Index        c);
  // ||-|--|
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // |-||--|
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // -|||--|
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // ||--|-|
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // |-|-|-|
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // -||-|-|
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // |--||-|
  Tensor3View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // -|-||-|
  Tensor3View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // --|||-|
  Tensor3View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // ||---||
  Tensor3View operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // |-|--||
  Tensor3View operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // -||--||
  Tensor3View operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // |--|-||
  Tensor3View operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // -|-|-||
  Tensor3View operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // --||-||
  Tensor3View operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // |---|||
  Tensor3View operator()( Index        l,
			  const Range& v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // -|--|||
  Tensor3View operator()( const Range& l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // --|-|||
  Tensor3View operator()( const Range& l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // ---||||
  Tensor3View operator()( const Range& l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);

  // Result 2D (6+5+4+3+2+1 = 21 combinations)
  // |||||--
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  Index        p, const Range& r, const Range& c);
  // ||||-|-
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  const Range& p, Index        r, const Range& c);
  // |||-||-
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  Index        p, Index        r, const Range& c);
  // ||-|||-
  MatrixView  operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  Index        p, Index        r, const Range& c);
  // |-||||-
  MatrixView  operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // -|||||-
  MatrixView  operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // ||||--|
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  const Range& p, const Range& r, Index        c);
  // |||-|-|
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  Index        p, const Range& r, Index        c);
  // ||-||-|
  MatrixView  operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  Index        p, const Range& r, Index        c);
  // |-|||-|
  MatrixView  operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // -||||-|
  MatrixView  operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // |||--||
  MatrixView  operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  const Range& p, Index        r, Index        c);
  // ||-|-||
  MatrixView  operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  const Range& p, Index        r, Index        c);
  // |-||-||
  MatrixView  operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // -|||-||
  MatrixView  operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // ||--|||
  MatrixView  operator()( Index        l,
			  Index        v, const Range& s, const Range& b,
			  Index        p, Index        r, Index        c);
  // |-|-|||
  MatrixView  operator()( Index        l,
			  const Range& v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // -||-|||
  MatrixView  operator()( const Range& l,
			  Index        v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // |--||||
  MatrixView  operator()( Index        l,
			  const Range& v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);
  // -|-||||
  MatrixView  operator()( const Range& l,
			  Index        v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);
  // --|||||
  MatrixView  operator()( const Range& l,
			  const Range& v, Index        s, Index        b,
			  Index        p, Index        r, Index        c);

  // Result 1D (7 combinations)
  // ||||||-
  VectorView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  Index        p, Index        r, const Range& c);
  // |||||-|
  VectorView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  Index        p, const Range& r, Index        c);
  // ||||-||
  VectorView  operator()( Index        l,
			  Index        v, Index        s, Index        b,
			  const Range& p, Index        r, Index        c);
  // |||-|||
  VectorView  operator()( Index        l,
			  Index        v, Index        s, const Range& b,
			  Index        p, Index        r, Index        c);
  // ||-||||
  VectorView  operator()( Index        l,
			  Index        v, const Range& s, Index        b,
			  Index        p, Index        r, Index        c);
  // |-|||||
  VectorView  operator()( Index        l,
			  const Range& v, Index        s, Index        b,
			  Index        p, Index        r, Index        c);
  // -||||||
  VectorView  operator()( const Range& l,
			  Index        v, Index        s, Index        b,
			  Index        p, Index        r, Index        c);

  // Result scalar (1 combination)
  // |||||||
  Numeric&         operator()( Index        l,
			       Index        v, Index        s, Index        b,
			       Index        p, Index        r, Index        c);


  // Functions returning const iterators:
  ConstIterator7D begin() const;
  ConstIterator7D end() const;
  // Functions returning iterators:
  Iterator7D begin();
  Iterator7D end();
  
  // Assignment operators:
  Tensor7View& operator=(const ConstTensor7View& v);
  Tensor7View& operator=(const Tensor7View& v);
  Tensor7View& operator=(const Tensor7& v);
  Tensor7View& operator=(Numeric x);

  // Other operators:
  Tensor7View& operator*=(Numeric x);
  Tensor7View& operator/=(Numeric x);
  Tensor7View& operator+=(Numeric x);
  Tensor7View& operator-=(Numeric x);

  Tensor7View& operator*=(const ConstTensor7View& x);
  Tensor7View& operator/=(const ConstTensor7View& x);
  Tensor7View& operator+=(const ConstTensor7View& x);
  Tensor7View& operator-=(const ConstTensor7View& x);

  // Friends:

  // Special constructor to make a Tensor7 view of a Tensor6.
  Tensor7View(const Tensor6View& a);

protected:
  // Constructors:
  Tensor7View();
  Tensor7View(Numeric *data,
	      const Range& l,
	      const Range& v, const Range& s, const Range& b,
	      const Range& p, const Range& r, const Range& c);
  Tensor7View(Numeric *data,
	      const Range& pl,
	      const Range& pv, const Range& ps, const Range& pb,
	      const Range& pp, const Range& pr, const Range& pc,
	      const Range& nl,
	      const Range& nv, const Range& ns, const Range& nb,
	      const Range& np, const Range& nr, const Range& nc);
};

/** The Tensor7 class. This is a Tensor7View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor7View. Additionally defined here
    are: 

    1. Constructors and destructor.
    2. Assignment operators.
    3. Resize function. */
class Tensor7 : public Tensor7View {
public:
  // Constructors:
  Tensor7();
  Tensor7(Index        l,
	  Index        v, Index        s, Index        b,
	  Index        p, Index        r, Index        c);
  Tensor7(Index        l,
	  Index        v, Index        s, Index        b,
	  Index        p, Index        r, Index        c,
	  Numeric fill);
  Tensor7(const ConstTensor7View& v);
  Tensor7(const Tensor7& v);

  // Assignment operators:
  Tensor7& operator=(const Tensor7& x);
  Tensor7& operator=(Numeric x);

  // Resize function:
  void resize(Index        l,
	      Index        v, Index        s, Index        b,
	      Index        p, Index        r, Index        c);

  // Destructor:
  ~Tensor7();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator7D origin,
                 const ConstIterator7D& end,
                 Iterator7D target);

inline void copy(Numeric x,
                 Iterator7D target,
                 const Iterator7D& end);

void transform( Tensor7View y,
                double (&my_func)(double),
                ConstTensor7View x );

Numeric max(const ConstTensor7View& x);

Numeric min(const ConstTensor7View& x);

std::ostream& operator<<(std::ostream& os, const ConstTensor7View& v);


// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Tensor7. */
typedef Array<Tensor7> ArrayOfTensor7;




#endif    // matpackVII_h
