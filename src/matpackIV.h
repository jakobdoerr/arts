/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>
                      Wolfram-Andre Haas <wolhaas@hermes.fho-emden.de>

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
  Implementation of Tensors of Rank 4.

  Based on Tensor3 by Stefan Buehler.

  The four dimensions are called: book, page, row, column.

  \author Wolfram-Andre Haas
  \date   2002-03-01
 */

#ifndef matpackIV_h
#define matpackIV_h

#include <iomanip>
#include "matpackI.h"
#include "matpackIII.h"

/** The outermost iterator class for rank 4 tensors. This takes into
    account the defined strided. */
class Iterator4D {
public:
  // Constructors:
  Iterator4D();
  Iterator4D(const Iterator4D& o);
  Iterator4D(const Tensor3View& x, Index stride);

  // Operators:
  Iterator4D& operator++();
  bool operator!=(const Iterator4D& other) const;
  Tensor3View* const operator->();
  Tensor3View& operator*();

private:
  /** Current position. */
  Tensor3View msv;
  /** Stride. */
  Index mstride;
};

/** Const version of Iterator4D. */
class ConstIterator4D {
public:
  // Constructors:
  ConstIterator4D();
  ConstIterator4D(const ConstIterator4D& o);
  ConstIterator4D(const ConstTensor3View& x, Index stride);

  // Operators:
  ConstIterator4D& operator++();
  bool operator!=(const ConstIterator4D& other) const;
  const ConstTensor3View* operator->() const;
  const ConstTensor3View& operator*()  const;

private:
  /** Current position. */
  ConstTensor3View msv;
  /** Stride. */
  Index mstride;
};


// Declare class Tensor4:
class Tensor4;


/** A constant view of a Tensor4.

    This, together with the derived class Tensor4View, contains the
    main implementation of a Tensor4. It defines the concepts of
    Tensor4View. Plus additionally the recursive subrange operator,
    which makes it possible to create a Tensor4View from a subrange of
    a Tensor4View.

    The four dimensions of the tensor are called: book, page, row, column.

    The class Tensor4 is just a special case of a Tensor4View
    which also allocates storage. */
class ConstTensor4View {
public:
  // Member functions:
  Index nbooks() const;
  Index npages() const;
  Index nrows()  const;
  Index ncols()  const;

  // Const index operators:
  ConstTensor4View operator()( const Range& b, const Range& p, const Range& r, const Range& c ) const;

  ConstTensor3View operator()( const Range& b, const Range& p, const Range& r, Index c        ) const;
  ConstTensor3View operator()( const Range& b, const Range& p, Index r,        const Range& c ) const;
  ConstTensor3View operator()( const Range& b, Index p,        const Range &r, const Range& c ) const;
  ConstTensor3View operator()( Index b,        const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& b, const Range& p, Index r,        Index c        ) const;
  ConstMatrixView  operator()( const Range& b, Index p,        const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& b, Index p,        Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index b,        const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index b,        const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( Index b,        Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( const Range& b, Index p,        Index r,        Index c        ) const;
  ConstVectorView  operator()( Index b,        const Range& p, Index r,        Index c        ) const;
  ConstVectorView  operator()( Index b,        Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( Index b,        Index p,        Index r,        const Range& c ) const;

  Numeric          operator()( Index b,        Index p,        Index r,        Index c        ) const;

  // Functions returning iterators:
  ConstIterator4D begin() const;
  ConstIterator4D end()   const;

  // Friends:
  friend class Tensor4View;
  friend class ConstIterator5D;
  friend class ConstTensor5View;
  friend class ConstTensor6View;
  friend class ConstTensor7View;

protected:
  // Constructors:
  ConstTensor4View();
  ConstTensor4View(Numeric *data,
                   const Range& b, const Range& p, const Range& r, const Range& c);
  ConstTensor4View(Numeric *data,
                   const Range& pb, const Range& pp, const Range& pr, const Range& pc,
                   const Range& nb, const Range& np, const Range& nr, const Range& nc);

  // Data members:
  // -------------
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

/** The Tensor4View class

    This contains the main implementation of a Tensor4. It defines
    the concepts of Tensor4View. Plus additionally the recursive
    subrange operator, which makes it possible to create a Tensor4View
    from a subrange of a Tensor4View.

    The class Tensor4 is just a special case of a Tensor4View
    which also allocates storage. */
class Tensor4View : public ConstTensor4View {
public:

  // Const index operators:
  ConstTensor4View operator()( const Range& b, const Range& p, const Range& r, const Range& c ) const;

  ConstTensor3View operator()( const Range& b, const Range& p, const Range& r, Index c        ) const;
  ConstTensor3View operator()( const Range& b, const Range& p, Index r,        const Range& c ) const;
  ConstTensor3View operator()( const Range& b, Index p,        const Range &r, const Range& c ) const;
  ConstTensor3View operator()( Index b,        const Range& p, const Range& r, const Range& c ) const;

  ConstMatrixView  operator()( const Range& b, const Range& p, Index r,        Index c        ) const;
  ConstMatrixView  operator()( const Range& b, Index p,        const Range& r, Index c        ) const;
  ConstMatrixView  operator()( const Range& b, Index p,        Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index b,        const Range& p, Index r,        const Range& c ) const;
  ConstMatrixView  operator()( Index b,        const Range& p, const Range& r, Index c        ) const;
  ConstMatrixView  operator()( Index b,        Index p,        const Range& r, const Range& c ) const;

  ConstVectorView  operator()( const Range& b, Index p,        Index r,        Index c        ) const;
  ConstVectorView  operator()( Index b,        const Range& p, Index r,        Index c        ) const;
  ConstVectorView  operator()( Index b,        Index p,        const Range& r, Index c        ) const;
  ConstVectorView  operator()( Index b,        Index p,        Index r,        const Range& c ) const;

  Numeric          operator()( Index b,        Index p,        Index r,        Index c        ) const;

  // Non-const index operators:

  Tensor4View operator()( const Range& b, const Range& p, const Range& r, const Range& c );

  Tensor3View operator()( const Range& b, const Range& p, const Range& r, Index c        );
  Tensor3View operator()( const Range& b, const Range& p, Index r,        const Range& c );
  Tensor3View operator()( const Range& b, Index p,        const Range &r, const Range& c );
  Tensor3View operator()( Index b,        const Range& p, const Range& r, const Range& c );

  MatrixView  operator()( const Range& b, const Range& p, Index r,        Index c        );
  MatrixView  operator()( const Range& b, Index p,        const Range& r, Index c        );
  MatrixView  operator()( const Range& b, Index p,        Index r,        const Range& c );
  MatrixView  operator()( Index b,        const Range& p, Index r,        const Range& c );
  MatrixView  operator()( Index b,        const Range& p, const Range& r, Index c        );
  MatrixView  operator()( Index b,        Index p,        const Range& r, const Range& c );

  VectorView  operator()( const Range& b, Index p,        Index r,        Index c        );
  VectorView  operator()( Index b,        const Range& p, Index r,        Index c        );
  VectorView  operator()( Index b,        Index p,        const Range& r, Index c        );
  VectorView  operator()( Index b,        Index p,        Index r,        const Range& c );

  Numeric&    operator()( Index b,        Index p,        Index r,        Index c        );

  // Functions returning const iterators:
  ConstIterator4D begin() const;
  ConstIterator4D end()   const;
  // Functions returning iterators:
  Iterator4D begin();
  Iterator4D end();

  // Assignment operators:
  Tensor4View& operator=(const ConstTensor4View& v);
  Tensor4View& operator=(const Tensor4View& v);
  Tensor4View& operator=(const Tensor4& v);
  Tensor4View& operator=(Numeric x);

  // Other operators:
  Tensor4View& operator*=(Numeric x);
  Tensor4View& operator/=(Numeric x);
  Tensor4View& operator+=(Numeric x);
  Tensor4View& operator-=(Numeric x);

  Tensor4View& operator*=(const ConstTensor4View& x);
  Tensor4View& operator/=(const ConstTensor4View& x);
  Tensor4View& operator+=(const ConstTensor4View& x);
  Tensor4View& operator-=(const ConstTensor4View& x);

  // Friends:
  // friend class VectorView;
  // friend ConstTensor4View transpose(ConstTensor4View m);
  // friend Tensor4View transpose(Tensor4View m);
  friend class Iterator5D;
  friend class Tensor5View;
  friend class Tensor6View;
  friend class Tensor7View;

protected:
  // Constructors:
  Tensor4View();
  Tensor4View(Numeric *data,
              const Range& b, const Range& p, const Range& r, const Range& c);
  Tensor4View(Numeric *data,
              const Range& pb, const Range& pp, const Range& pr, const Range& pc,
              const Range& nb, const Range& np, const Range& nr, const Range& nc);
};

/** The Tensor4 class. This is a Tensor4View that also allocates storage
    automatically, and deallocates it when it is destroyed. We take
    all the functionality from Tensor4View. Additionally defined here
    are:

    1. Constructors and destructor.
    2. Assignment operators.
    3. Resize function. */
class Tensor4 : public Tensor4View {
public:
  // Constructors:
  Tensor4();
  Tensor4(Index b, Index p, Index r, Index c);
  Tensor4(Index b, Index p, Index r, Index c, Numeric fill);
  Tensor4(const ConstTensor4View& v);
  Tensor4(const Tensor4& v);

  // Assignment operators:
  Tensor4& operator=(const Tensor4& x);
  Tensor4& operator=(Numeric x);

  // Resize function:
  void resize(Index b, Index p, Index r, Index c);

  // Destructor:
  ~Tensor4();
};


// Function declarations:
// ----------------------

inline void copy(ConstIterator4D origin,
                 const ConstIterator4D& end,
                 Iterator4D target);

inline void copy(Numeric x,
                 Iterator4D target,
                 const Iterator4D& end);


// Declare the existance of class Array:
template<class base>
class Array;

/** An array of Tensor4. */
typedef Array<Tensor4> ArrayOfTensor4;


// Functions for Iterator4D
// ------------------------

/** Default constructor. */
inline Iterator4D::Iterator4D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline Iterator4D::Iterator4D(const Iterator4D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline Iterator4D::Iterator4D(const Tensor3View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here.
}

/** Prefix increment operator. */
inline Iterator4D& Iterator4D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool Iterator4D::operator!=(const Iterator4D& other) const
{
  if ( msv.mdata +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart
       !=
       other.msv.mdata +
       other.msv.mpr.mstart +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 3D iterators. */
inline Tensor3View* const Iterator4D::operator->()
{
  return &msv;
}

/** Dereferencing. */
inline Tensor3View& Iterator4D::operator*()
{
  return msv;
}

// Functions for ConstIterator4D
// -----------------------------

/** Default constructor. */
inline ConstIterator4D::ConstIterator4D()
{
  // Nothing to do here.
}

/** Copy constructor. */
inline ConstIterator4D::ConstIterator4D(const ConstIterator4D& o) :
  msv(o.msv), mstride(o.mstride)
{
  // Nothing to do here.
}

/** Explicit constructor. */
inline ConstIterator4D::ConstIterator4D(const ConstTensor3View& x, Index stride) :
  msv(x), mstride(stride)
{
  // Nothing to do here.
}

/** Prefix increment operator. */
inline ConstIterator4D& ConstIterator4D::operator++()
{
  msv.mdata += mstride;
  return *this;
}

/** Not equal operator, needed for algorithms like copy. */
inline bool ConstIterator4D::operator!=(const ConstIterator4D& other) const
{
  if ( msv.mdata +
       msv.mpr.mstart +
       msv.mrr.mstart +
       msv.mcr.mstart
       !=
       other.msv.mdata +
       other.msv.mpr.mstart +
       other.msv.mrr.mstart +
       other.msv.mcr.mstart )
    return true;
  else
    return false;
}

/** The -> operator is needed, so that we can write i->begin() to get
    the 3D iterators. */
inline const ConstTensor3View* ConstIterator4D::operator->() const
{
  return &msv;
}

/** Dereferencing. */
inline const ConstTensor3View& ConstIterator4D::operator*() const
{
  return msv;
}



// Functions for ConstTensor4View:
// ------------------------------

/** Returns the number of books. */
inline Index ConstTensor4View::nbooks() const
{
  return mbr.mextent;
}

/** Returns the number of pages. */
inline Index ConstTensor4View::npages() const
{
  return mpr.mextent;
}

/** Returns the number of rows. */
inline Index ConstTensor4View::nrows() const
{
  return mrr.mextent;
}

/** Returns the number of columns. */
inline Index ConstTensor4View::ncols() const
{
  return mcr.mextent;
}

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor4. This allows
    correct recursive behavior.  */
inline ConstTensor4View ConstTensor4View::operator()(const Range& b,
                                                     const Range& p,
                                                     const Range& r,
                                                     const Range& c) const
{
  return ConstTensor4View(mdata,
                          mbr, mpr, mrr, mcr,
                          b,   p,   r,   c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) */
inline ConstTensor3View ConstTensor4View::operator()(const Range& b,
                                                     const Range& p,
                                                     const Range& r,
                                                     Index c       ) const
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return ConstTensor3View(mdata +
                          mcr.mstart + c * mcr.mstride,
                          mbr, mpr, mrr,
                          b,   p,   r);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) */
inline ConstTensor3View ConstTensor4View::operator()(const Range& b,
                                                     const Range& p,
                                                     Index r,
                                                     const Range& c) const
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r < mrr.mextent );

  return ConstTensor3View(mdata +
                          mrr.mstart + r * mrr.mstride,
                          mbr, mpr, mcr,
                          b,   p,   c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) */
inline ConstTensor3View ConstTensor4View::operator()(const Range& b,
                                                     Index p,
                                                     const Range& r,
                                                     const Range& c) const
{
  // Check that p is valid:
  assert( 0 <= p );
  assert( p < mpr.mextent );

  return ConstTensor3View(mdata +
                          mpr.mstart + p * mpr.mstride,
                          mbr, mrr, mcr,
                          b,   r,   c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) */
inline ConstTensor3View ConstTensor4View::operator()(Index b,
                                                     const Range& p,
                                                     const Range& r,
                                                     const Range& c) const
{
  // Check that b is valid:
  assert( 0 <= b );
  assert( b < mbr.mextent );

  return ConstTensor3View(mdata +
                          mbr.mstart + b * mbr.mstride,
                          mpr, mrr, mcr,
                          p,   r,   c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
inline ConstMatrixView ConstTensor4View::operator()(const Range& b,
                                                    const Range& p,
                                                    Index r,
                                                    Index c       ) const
{
  // Check that r and c are valid:
  assert( 0 <= r );
  assert( 0 <= c );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata +
                         mrr.mstart + r * mrr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mbr, mpr,
                         b,   p);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
inline ConstMatrixView ConstTensor4View::operator()(const Range& b,
                                                    Index p,
                                                    const Range& r,
                                                    Index c       ) const
{
  // Check that p and c are valid:
  assert( 0 <= p );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata +
                         mpr.mstart + p * mpr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mbr, mrr,
                         b,   r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
inline ConstMatrixView ConstTensor4View::operator()(const Range& b,
                                                    Index p,
                                                    Index r,
                                                    const Range& c) const
{
  // Check that p and r are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return ConstMatrixView(mdata +
                         mpr.mstart + p * mpr.mstride +
                         mrr.mstart + r * mrr.mstride,
                         mbr, mcr,
                         b,   c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
inline ConstMatrixView ConstTensor4View::operator()(Index b,
                                                    const Range& p,
                                                    Index r,
                                                    const Range& c) const
{
  // Check that b and r are valid:
  assert( 0 <= b );
  assert( 0 <= r );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );

  return ConstMatrixView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mrr.mstart + r * mrr.mstride,
                         mpr, mcr,
                         p,   c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
inline ConstMatrixView ConstTensor4View::operator()(Index b,
                                                    const Range& p,
                                                    const Range& r,
                                                    Index c       ) const
{
  // Check that b and c are valid:
  assert( 0 <= b );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( c < mcr.mextent );

  return ConstMatrixView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mpr, mrr,
                         p,   r);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) */
inline ConstMatrixView ConstTensor4View::operator()(Index b,
                                                    Index p,
                                                    const Range& r,
                                                    const Range& c) const
{
  // Check that b and p are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );

  return ConstMatrixView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mpr.mstart + p * mpr.mstride,
                         mrr, mcr,
                         r,   c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) */
inline ConstVectorView ConstTensor4View::operator()(const Range& b,
                                                    Index p,
                                                    Index r,
                                                    Index c       ) const
{
  // Check that p, r and c are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata +
                         mpr.mstart + p * mpr.mstride +
                         mrr.mstart + r * mrr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mbr,
                         b);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) */
inline ConstVectorView ConstTensor4View::operator()(Index b,
                                                    const Range& p,
                                                    Index r,
                                                    Index c       ) const
{
  // Check that b, r and c are valid:
  assert( 0 <= b );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mrr.mstart + r * mrr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mpr,
                         p);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) */
inline ConstVectorView ConstTensor4View::operator()(Index b,
                                                    Index p,
                                                    const Range& r,
                                                    Index c       ) const
{
  // Check that b, p and c are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return ConstVectorView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mcr.mstart + c * mcr.mstride,
                         mrr,
                         r);
}

/** Const index operator returning an object of type
    ConstVectorView. Reducing the dimension by three.) */
inline ConstVectorView ConstTensor4View::operator()(Index b,
                                                    Index p,
                                                    Index r,
                                                    const Range& c) const
{
  // Check that b, p and r are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return ConstVectorView(mdata +
                         mbr.mstart + b * mbr.mstride +
                         mpr.mstart + p * mpr.mstride +
                         mrr.mstart + r * mrr.mstride,
                         mcr,
                         c);
}

/** Plain const index operator. */
inline Numeric ConstTensor4View::operator()(Index b,
                                            Index p,
                                            Index r,
                                            Index c) const
{
  // Check if indices are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return *( mdata +
            mbr.mstart + b * mbr.mstride +
            mpr.mstart + p * mpr.mstride +
            mrr.mstart + r * mrr.mstride +
            mcr.mstart + c * mcr.mstride );
}

/** Return const iterator to first book. */
inline ConstIterator4D ConstTensor4View::begin() const
{
  return ConstIterator4D( ConstTensor3View(mdata + mbr.mstart,
                                           mpr, mrr, mcr),
                          mbr.mstride );
}

/** Return const iterator behind last book. */
inline ConstIterator4D ConstTensor4View::end() const
{
  return ConstIterator4D( ConstTensor3View(mdata + mbr.mstart +
                                           (mbr.mextent) * mbr.mstride,
                                           mpr, mrr, mcr),
                          mbr.mstride );
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for derived classes. */
inline ConstTensor4View::ConstTensor4View() :
  mbr(0,0,1),
  mpr(0,0,1),
  mrr(0,0,1),
  mcr(0,0,1),
  mdata(NULL)
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor4 to initialize
    its own Tensor4View part. The page range pr must have a stride to
    account for the length of one page. The book range br must have a
    stride to account for the length of one book. */
inline ConstTensor4View::ConstTensor4View(Numeric *data,
                                          const Range& br,
                                          const Range& pr,
                                          const Range& rr,
                                          const Range& cr) :
  mbr(br),
  mpr(pr),
  mrr(rr),
  mcr(cr),
  mdata(data)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct sub-tensors from
    sub-tensors. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range. */
inline ConstTensor4View::ConstTensor4View(Numeric *data,
                                          const Range& pb,
                                          const Range& pp,
                                          const Range& pr,
                                          const Range& pc,
                                          const Range& nb,
                                          const Range& np,
                                          const Range& nr,
                                          const Range& nc) :
  mbr(pb,nb),
  mpr(pp,np),
  mrr(pr,nr),
  mcr(pc,nc),
  mdata(data)
{
  // Nothing to do here.
}

/** Output operator. This demonstrates how iterators can be used to
    traverse the tensor. We use the standard output operator for
    Tensor to print each book in turn. */
inline std::ostream& operator<<(std::ostream& os, const ConstTensor4View& v)
{
  // Page iterators:
  ConstIterator4D ib = v.begin();
  const ConstIterator4D end_book = v.end();

  if ( ib != end_book ) {
    os << *ib;
    ++ib;
  }

  for ( ; ib != end_book; ++ib ) {
    os << "\n\n";
    os << *ib;
  }

  return os;
}


// Functions for Tensor4View:
// -------------------------

/** Const index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor4. This allows
    correct recursive behavior. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline ConstTensor4View Tensor4View::operator()(const Range& b,
                                                const Range& p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor4View::operator()(const Range& b,
                                                const Range& p,
                                                const Range& r,
                                                Index c       ) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor4View::operator()(const Range& b,
                                                const Range& p,
                                                Index r,
                                                const Range& c) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor4View::operator()(const Range& b,
                                                Index p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstTensor3View. (Reducing the dimension by one.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstTensor3View Tensor4View::operator()(Index b,
                                                const Range& p,
                                                const Range& r,
                                                const Range& c) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor4View::operator()(const Range& b,
                                               const Range& p,
                                               Index r,
                                               Index c       ) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor4View::operator()(const Range& b,
                                               Index p,
                                               const Range& r,
                                               Index c       ) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor4View::operator()(const Range& b,
                                               Index p,
                                               Index r,
                                               const Range& c) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor4View::operator()(Index b,
                                               const Range& p,
                                               Index r,
                                               const Range& c) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor4View::operator()(Index b,
                                               const Range& p,
                                               const Range& r,
                                               Index c       ) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstMatrixView. (Reducing the dimension by two.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstMatrixView Tensor4View::operator()(Index b,
                                               Index p,
                                               const Range& r,
                                               const Range& c) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor4View::operator()(const Range& b,
                                               Index p,
                                               Index r,
                                               Index c       ) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor4View::operator()(Index b,
                                               const Range& p,
                                               Index r,
                                               Index c       ) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor4View::operator()(Index b,
                                               Index p,
                                               const Range& r,
                                               Index c       ) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Const index operator returning an object of type
    ConstVectorView. (Reducing the dimension by three.) Has to be
    redefined here, since it is hiden by the non-const operator of the
    derived class. */
inline ConstVectorView Tensor4View::operator()(Index b,
                                               Index p,
                                               Index r,
                                               const Range& c) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Plain const index operator. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class. */
inline Numeric Tensor4View::operator()(Index b,
                                       Index p,
                                       Index r,
                                       Index c) const
{
  return ConstTensor4View::operator()(b,p,r,c);
}

/** Index operator for subrange. We have to also account for the
    case, that *this is already a subrange of a Tensor4. This allows
    correct recursive behavior.  */
inline Tensor4View Tensor4View::operator()(const Range& b,
                                           const Range& p,
                                           const Range& r,
                                           const Range& c)
{
  return Tensor4View(mdata,
                     mbr, mpr, mrr, mcr,
                     b,   p,   r,   c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by one.) */
inline Tensor3View Tensor4View::operator()(const Range& b,
                                           const Range& p,
                                           const Range& r,
                                           Index c)
{
  // Check that c is valid:
  assert( 0 <= c );
  assert( c < mcr.mextent );

  return Tensor3View(mdata +
                     mcr.mstart + c * mcr.mstride,
                     mbr, mpr, mrr,
                     b,   p,   r);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by one.) */
inline Tensor3View Tensor4View::operator()(const Range& b,
                                           const Range& p,
                                           Index r,
                                           const Range& c)
{
  // Check that r is valid:
  assert( 0 <= r );
  assert( r < mrr.mextent );

  return Tensor3View(mdata +
                     mrr.mstart + r * mrr.mstride,
                     mbr, mpr, mcr,
                     b,   p,   c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by one.) */
inline Tensor3View Tensor4View::operator()(const Range& b,
                                           Index p,
                                           const Range& r,
                                           const Range& c)
{
  // Check that p is valid:
  assert( 0 <= p );
  assert( p < mpr.mextent );

  return Tensor3View(mdata +
                     mpr.mstart + p * mpr.mstride,
                     mbr, mrr, mcr,
                     b,   r,   c);
}

/** Index operator returning an object of type
    Tensor3View. (Reducing the dimension by one.) */
inline Tensor3View Tensor4View::operator()(Index b,
                                           const Range& p,
                                           const Range& r,
                                           const Range& c)
{
  // Check that b is valid:
  assert( 0 <= b );
  assert( b < mbr.mextent );

  return Tensor3View(mdata +
                     mbr.mstart + b * mbr.mstride,
                     mpr, mrr, mcr,
                     p,   r,   c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
inline MatrixView Tensor4View::operator()(const Range& b,
                                          const Range& p,
                                          Index r,
                                          Index c)
{
  // Check that r and c are valid:
  assert( 0 <= r );
  assert( 0 <= c );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return MatrixView(mdata +
                    mrr.mstart + r * mrr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mbr, mpr,
                    b,   p);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
inline MatrixView Tensor4View::operator()(const Range& b,
                                          Index p,
                                          const Range& r,
                                          Index c)
{
  // Check that p and c are valid:
  assert( 0 <= p );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return MatrixView(mdata +
                    mpr.mstart + p * mpr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mbr, mrr,
                    b,   r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
inline MatrixView Tensor4View::operator()(const Range& b,
                                          Index p,
                                          Index r,
                                          const Range& c)
{
  // Check that p and r are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return MatrixView(mdata +
                    mpr.mstart + p * mpr.mstride +
                    mrr.mstart + r * mrr.mstride,
                    mbr, mcr,
                    b,   c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
inline MatrixView Tensor4View::operator()(Index b,
                               const Range& p,
                               Index r,
                               const Range& c)
{
  // Check that b and r are valid:
  assert( 0 <= b );
  assert( 0 <= r );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );

  return MatrixView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mrr.mstart + r * mrr.mstride,
                    mpr, mcr,
                    p,   c);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
inline MatrixView Tensor4View::operator()(Index b,
                                          const Range& p,
                                          const Range& r,
                                          Index c)
{
  // Check that b and c are valid:
  assert( 0 <= b );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( c < mcr.mextent );

  return MatrixView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mpr, mrr,
                    p,   r);
}

/** Index operator returning an object of type
    MatrixView. (Reducing the dimension by two.) */
inline MatrixView Tensor4View::operator()(Index b,
                                          Index p,
                                          const Range& r,
                                          const Range& c)
{
  // Check that b and p are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );

  return MatrixView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mpr.mstart + p * mpr.mstride,
                    mrr, mcr,
                    r,   c);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by three.) */
inline VectorView Tensor4View::operator()(const Range& b,
                                          Index p,
                                          Index r,
                                          Index c)
{
  // Check that p, r and c are valid:
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return VectorView(mdata +
                    mpr.mstart + p * mpr.mstride +
                    mrr.mstart + r * mrr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mbr,
                    b);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by three.) */
inline VectorView Tensor4View::operator()(Index b,
                                          const Range& p,
                                          Index r,
                                          Index c)
{
  // Check that b, r and c are valid:
  assert( 0 <= b );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return VectorView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mrr.mstart + r * mrr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mpr,
                    p);
}

/** Index operator returning an object of type
    VectorView. (Reducing the dimension by three.) */
inline VectorView Tensor4View::operator()(Index b,
                                          Index p,
                                          const Range& r,
                                          Index c)
{
  // Check that b, p and c are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( c < mcr.mextent );

  return VectorView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mcr.mstart + c * mcr.mstride,
                    mrr,
                    r);
}

/** Index operator returning an object of type
    VectorView. Reducing the dimension by three.) */
inline VectorView Tensor4View::operator()(Index b,
                                          Index p,
                                          Index r,
                                          const Range& c)
{
  // Check that b, p and r are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );

  return VectorView(mdata +
                    mbr.mstart + b * mbr.mstride +
                    mpr.mstart + p * mpr.mstride +
                    mrr.mstart + r * mrr.mstride,
                    mcr,
                    c);
}

/** Plain non-const index operator. */
inline Numeric& Tensor4View::operator()(Index b,
                                        Index p,
                                        Index r,
                                        Index c)
{
  // Check if indices are valid:
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );
  assert( b < mbr.mextent );
  assert( p < mpr.mextent );
  assert( r < mrr.mextent );
  assert( c < mcr.mextent );

  return *( mdata +
            mbr.mstart + b * mbr.mstride +
            mpr.mstart + p * mpr.mstride +
            mrr.mstart + r * mrr.mstride +
            mcr.mstart + c * mcr.mstride );
}

/** Return const iterator to first book. Has to be redefined here, since it is
    hiden by the non-const operator of the derived class.*/
inline ConstIterator4D Tensor4View::begin() const
{
  return ConstTensor4View::begin();
}

/** Return const iterator behind last book. */
inline ConstIterator4D Tensor4View::end() const
{
  return ConstTensor4View::end();
}

/** Return iterator to first book. */
inline Iterator4D Tensor4View::begin()
{
  return Iterator4D( Tensor3View(mdata + mbr.mstart,
                                 mpr, mrr, mcr),
                     mbr.mstride );
}

/** Return iterator behind last book. */
inline Iterator4D Tensor4View::end()
{
  return Iterator4D( Tensor3View(mdata + mbr.mstart +
                                 (mbr.mextent) * mbr.mstride,
                                 mpr, mrr, mcr),
                     mbr.mstride );
}

/** Assignment operator. This copies the data from another Tensor4View
    to this Tensor4View. Dimensions must agree! Resizing would destroy
    the selection that we might have done in this Tensor4View by
    setting its range. */
inline Tensor4View& Tensor4View::operator=(const ConstTensor4View& m)
{
  // Check that sizes are compatible:
  assert( mbr.mextent == m.mbr.mextent );
  assert( mpr.mextent == m.mpr.mextent );
  assert( mrr.mextent == m.mrr.mextent );
  assert( mcr.mextent == m.mcr.mextent );

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from Tensor4View to Tensor4View. This is a tricky
    one. The problem is that since Tensor4View is derived from
    ConstTensor4View, a default = operator is generated by the
    compiler, which does not do what we want. So we need this one to
    override the default. */
inline Tensor4View& Tensor4View::operator=(const Tensor4View& m)
{
  // Check that sizes are compatible:
  assert( mbr.mextent == m.mbr.mextent );
  assert( mpr.mextent == m.mpr.mextent );
  assert( mrr.mextent == m.mrr.mextent );
  assert( mcr.mextent == m.mcr.mextent );

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment from a Tensor4. This must exist to overide the
    automatically generated assignment operators, which don't copy the
    contents! */
inline Tensor4View& Tensor4View::operator=(const Tensor4& m)
{
  // Check that sizes are compatible:
  assert( mbr.mextent == m.mbr.mextent );
  assert( mpr.mextent == m.mpr.mextent );
  assert( mrr.mextent == m.mrr.mextent );
  assert( mcr.mextent == m.mcr.mextent );

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assigning a scalar to a Tensor4View will set all elements to this
    value. */
inline Tensor4View& Tensor4View::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

// Some little helper functions:
//------------------------------

/** Multiplication by scalar. */
inline Tensor4View& Tensor4View::operator*=(Numeric x)
{
  const Iterator4D eb = end();
  for ( Iterator4D b = begin(); b != eb ; ++b )
  {
    *b *= x;
  }
  return *this;
}

/** Division by scalar. */
inline Tensor4View& Tensor4View::operator/=(Numeric x)
{
  const Iterator4D eb = end();
  for ( Iterator4D b = begin(); b != eb ; ++b )
  {
    *b /= x;
  }
  return *this;
}

/** Addition of scalar. */
inline Tensor4View& Tensor4View::operator+=(Numeric x)
{
  const Iterator4D eb = end();
  for ( Iterator4D b = begin(); b != eb ; ++b )
  {
    *b += x;
  }
  return *this;
}

/** Subtraction of scalar. */
inline Tensor4View& Tensor4View::operator-=(Numeric x)
{
  const Iterator4D eb = end();
  for ( Iterator4D b = begin(); b != eb ; ++b )
  {
    *b -= x;
  }
  return *this;
}

/** Element-vise multiplication by another Tensor4. */
inline Tensor4View& Tensor4View::operator*=(const ConstTensor4View& x)
{
  assert( nbooks() == x.nbooks() );
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()  );
  assert( ncols()  == x.ncols()  );
  ConstIterator4D  xb = x.begin();
  Iterator4D        b = begin();
  const Iterator4D eb = end();
  for ( ; b != eb ; ++b, ++xb )
    {
      *b *= *xb;
    }
  return *this;
}

/** Element-vise division by another Tensor4. */
inline Tensor4View& Tensor4View::operator/=(const ConstTensor4View& x)
{
  assert( nbooks() == x.nbooks() );
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()  );
  assert( ncols()  == x.ncols()  );
  ConstIterator4D  xb = x.begin();
  Iterator4D        b = begin();
  const Iterator4D eb = end();
  for ( ; b != eb ; ++b, ++xb )
    {
      *b /= *xb;
    }
  return *this;
}

/** Element-vise addition of another Tensor4. */
inline Tensor4View& Tensor4View::operator+=(const ConstTensor4View& x)
{
  assert( nbooks() == x.nbooks() );
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()  );
  assert( ncols()  == x.ncols()  );
  ConstIterator4D  xb = x.begin();
  Iterator4D        b = begin();
  const Iterator4D eb = end();
  for ( ; b != eb ; ++b, ++xb )
    {
      *b += *xb;
    }
  return *this;
}

/** Element-vise subtraction of another Tensor4. */
inline Tensor4View& Tensor4View::operator-=(const ConstTensor4View& x)
{
  assert( nbooks() == x.nbooks() );
  assert( npages() == x.npages() );
  assert( nrows()  == x.nrows()  );
  assert( ncols()  == x.ncols()  );
  ConstIterator4D  xb = x.begin();
  Iterator4D        b = begin();
  const Iterator4D eb = end();
  for ( ; b != eb ; ++b, ++xb )
    {
      *b -= *xb;
    }
  return *this;
}

/** Default constructor. This is necessary, so that we can have a
    default constructor for the derived class Tensor4. */
inline Tensor4View::Tensor4View() :
  ConstTensor4View()
{
  // Nothing to do here.
}

/** Explicit constructor. This one is used by Tensor4 to initialize its
    own Tensor4View part. The row range rr must have a
    stride to account for the length of one row. */
inline Tensor4View::Tensor4View(Numeric *data,
                                const Range& br,
                                const Range& pr,
                                const Range& rr,
                                const Range& cr) :
  ConstTensor4View(data, br, pr, rr, cr)
{
  // Nothing to do here.
}

/** Recursive constructor. This is used to construct SubMatrices from
    SubMatrices. That means that the new ranges have to be interpreted
    relative to the original ranges.

    The new ranges may contain -1 for the extent which acts as a
    joker. However, the used Range constructor converts this to an
    explicit range, consistent with the original Range.

    \param *data The actual data.
    \param p Previous range.
    \param n New Range.  */
inline Tensor4View::Tensor4View(Numeric *data,
                                const Range& pb,
                                const Range& pp,
                                const Range& pr,
                                const Range& pc,
                                const Range& nb,
                                const Range& np,
                                const Range& nr,
                                const Range& nc) :
  ConstTensor4View(data, pb, pp, pr, pc, nb, np, nr, nc)
{
  // Nothing to do here.
}

/** Copy data between begin and end to target. Target must be a valid
    area of memory. Note that the strides in the iterators can be
    different, so that we can copy data between different
    kinds of subtensors. */
inline void copy(ConstIterator4D origin,
                 const ConstIterator4D& end,
                 Iterator4D target)
{
  for ( ; origin != end ; ++origin, ++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy( origin->begin(), origin->end(), target->begin() );
    }
}

/** Copy a scalar to all elements. */
inline void copy(Numeric x,
                 Iterator4D target,
                 const Iterator4D& end)
{
  for ( ; target != end ; ++target )
    {
      // We use the copy function for the next smaller rank of tensor
      // recursively:
      copy( x, target->begin(), target->end() );
    }
}


// Functions for Tensor4:
// ---------------------

/** Default constructor. */
inline Tensor4::Tensor4() :
  Tensor4View::Tensor4View()
{
  // Nothing to do here. However, note that the default constructor
  // for Tensor4View has been called in the initializer list. That is
  // crucial, otherwise internal range objects will not be properly
  // initialized.
}

/** Constructor setting size. This constructor has to set the strides
    in the book, page and row ranges correctly! */
inline Tensor4::Tensor4(Index b, Index p, Index r, Index c) :
  Tensor4View( new Numeric[b*p*r*c],
               Range( 0, b, p*r*c ),
               Range( 0, p, r*c ),
               Range( 0, r, c ),
               Range( 0, c) )
{
  // Nothing to do here.
}

/** Constructor setting size and filling with constant value. */
inline Tensor4::Tensor4(Index b, Index p, Index r, Index c, Numeric fill) :
  Tensor4View( new Numeric[b*p*r*c],
               Range( 0, b, p*r*c ),
               Range( 0, p, r*c ),
               Range( 0, r, c ),
               Range( 0, c) )
{
  // Here we can access the raw memory directly, for slightly
  // increased efficiency:
  for ( Numeric *x = mdata; x < mdata + b*p*r*c; ++x )
    *x = fill;
}

/** Copy constructor from Tensor4View. This automatically sets the size
    and copies the data. */
inline Tensor4::Tensor4(const ConstTensor4View& m) :
  Tensor4View( new Numeric[m.nbooks()*m.npages()*m.nrows()*m.ncols()],
               Range( 0, m.nbooks(), m.npages()*m.nrows()*m.ncols() ),
               Range( 0, m.npages(), m.nrows()*m.ncols() ),
               Range( 0, m.nrows(), m.ncols() ),
               Range( 0, m.ncols() ) )
{
  copy( m.begin(), m.end(), begin() );
}

/** Copy constructor from Tensor4. This automatically sets the size
    and copies the data. */
inline Tensor4::Tensor4(const Tensor4& m) :
  Tensor4View( new Numeric[m.nbooks()*m.npages()*m.nrows()*m.ncols()],
               Range( 0, m.nbooks(), m.npages()*m.nrows()*m.ncols() ),
               Range( 0, m.npages(), m.nrows()*m.ncols() ),
               Range( 0, m.nrows(), m.ncols() ),
               Range( 0, m.ncols() ) )
{
  // There is a catch here: If m is an empty tensor, then it will have
  // dimensions of size 0. But these are used to initialize the stride
  // for higher dimensions! Thus, this method has to be consistent
  // with the behaviour of Range::Range. For now, Range::Range allows
  // also stride 0.
  copy( m.begin(), m.end(), begin());
}

/** Assignment operator from another tensor. It is important that this
    operator exists. Otherwise the = operator seems to copy references
    instead of content in some cases.

    The Behavior of this one is a bit special: If the size of the
    target tensor is 0 then it will be automatically resized to match
    (this is needed to have the correct initialization for constructed
    classes that use the assignment operator to initialize their
    data).
*/
inline Tensor4& Tensor4::operator=(const Tensor4& m)
{
  //  cout << "Tensor4 copy: m = " << m.nrows() << " " << m.ncols() << "\n";
  //  cout << "              n = " << nrows() << " " << ncols() << "\n";

  // None of the extents can be zero for a valid tensor, so we just
  // have to check one.
  if ( 0 == mcr.mextent )
    {
      // Adjust if previously empty.
      resize( m.mbr.mextent, m.mpr.mextent, m.mrr.mextent, m.mcr.mextent );
    }
  else
    {
      // Check that sizes are compatible:
      assert( mbr.mextent == m.mbr.mextent );
      assert( mpr.mextent == m.mpr.mextent );
      assert( mrr.mextent == m.mrr.mextent );
      assert( mcr.mextent == m.mcr.mextent );
    }

  copy( m.begin(), m.end(), begin() );
  return *this;
}

/** Assignment operator from scalar. Assignment operators are not
    inherited. */
inline Tensor4& Tensor4::operator=(Numeric x)
{
  copy( x, begin(), end() );
  return *this;
}

/** Resize function. If the size is already correct this function does
    nothing. All data is lost after resizing! The new tensor is not
    initialized, so it will contain random values.*/
inline void Tensor4::resize(Index b, Index p, Index r, Index c)
{
  assert( 0 <= b );
  assert( 0 <= p );
  assert( 0 <= r );
  assert( 0 <= c );

  if ( mbr.mextent != b ||
       mpr.mextent != p ||
       mrr.mextent != r ||
       mcr.mextent != c )
    {
      delete mdata;
      mdata = new Numeric[b*p*r*c];

      mbr.mstart = 0;
      mbr.mextent = b;
      mbr.mstride = p*r*c;

      mpr.mstart = 0;
      mpr.mextent = p;
      mpr.mstride = r*c;

      mrr.mstart = 0;
      mrr.mextent = r;
      mrr.mstride = c;

      mcr.mstart = 0;
      mcr.mextent = c;
      mcr.mstride = 1;
    }
}

/** Destructor for Tensor4. This is important, since Tensor4 uses new to
    allocate storage. */
inline Tensor4::~Tensor4()
{
//   cout << "Destroying a Tensor4:\n"
//        << *this << "\n........................................\n";
  delete mdata;
}


/** A generic transform function for tensors, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for tensors! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

    transform(y,sin,x) computes y = sin(x)

    The two views may be the same one, in which case the
    conversion happens in place.

    \retval   y   the results of the function acting on each element of x
    \param    my_func a function (e.g., sqrt)
    \param    x   a tensor */
inline void transform( Tensor4View y,
                       double (&my_func)(double),
                       ConstTensor4View x )
{
  // Check dimensions:
  assert( y.nbooks() == x.nbooks() );
  assert( y.npages() == x.npages() );
  assert( y.nrows()  == x.nrows()  );
  assert( y.ncols()  == x.ncols()  );

  const ConstIterator4D xe = x.end();
  ConstIterator4D       xi = x.begin();
  Iterator4D            yi = y.begin();
  for ( ; xi != xe; ++xi, ++yi )
    {
      // Use the transform function of lower dimensional tensors
      // recursively:
      transform( *yi, my_func, *xi );
    }
}

/** Max function, tensor version. */
inline Numeric max(const ConstTensor4View& x)
{
  const ConstIterator4D xe = x.end();
  ConstIterator4D       xi = x.begin();

  // Initial value for max:
  Numeric themax = max(*xi);
  ++xi;

  for ( ; xi != xe ; ++xi )
    {
      // Use the max function of lower dimensional tensors
      // recursively:
      Numeric maxi = max(*xi);
      if ( maxi > themax )
        themax = maxi;
    }

  return themax;
}

/** Min function, tensor version. */
inline Numeric min(const ConstTensor4View& x)
{
  const ConstIterator4D xe = x.end();
  ConstIterator4D       xi = x.begin();

  // Initial value for min:
  Numeric themin = min(*xi);
  ++xi;

  for ( ; xi != xe ; ++xi )
    {
      // Use the min function of lower dimensional tensors
      // recursively:
      Numeric mini = min(*xi);
      if ( mini < themin )
        themin = mini;
    }

  return themin;
}

#endif    // matpackIV_h
