//                              -*- Mode: C++ -*- 
// BasePair.h
// Copyright Â© 2006 Laboratoire de Biologie Informatique et ThÃ©orique
//                  UniversitÃ© de MontrÃ©al.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Wed Jul 12 16:27:16 2006
// $Revision$
// $Id$
// 


#ifndef _annotate_BasePair_h_
#define _annotate_BasePair_h_

#include <algorithm>
#include <utility>

#include "mccore/GraphModel.h"
#include "mccore/ResId.h"

using namespace mccore;
using namespace std;



namespace annotate
{

  class BasePair : public pair< GraphModel::label, GraphModel::label >
  {

  public:

    ResId fResId;

    ResId rResId;
    
    // LIFECYCLE ------------------------------------------------------------

    BasePair (GraphModel::label l, const ResId &fResId, GraphModel::label r, const ResId &rResId)
      : pair< GraphModel::label, GraphModel::label > (l, r),
	fResId (fResId),
	rResId (rResId)
    { }

    ~BasePair () { }
    
    // OPERATORS ------------------------------------------------------------

    /**
     * Assigns the right BasePair to this.
     * @param right the BasePair to copy into this.
     * @return a reference to this.
     */
    BasePair& operator= (const BasePair &right)
    {
      if (this != &right)
	{
	  first = right.first;
	  second = right.second;
	  fResId = right.fResId;
	  rResId = right.rResId;
	}
      return *this;
    }

    /**
     * Compares this BasePair to right.
     * @param right the BasePair to compare with this.
     * @return true if both BasePair are equal.
     */
    bool operator== (const BasePair &right) const
    {
      return (&right == this
	      || (first == right.first && second == right.second));
    }

    /**
     * Compares this BasePair to right.
     * @param right the BasePair to compare with this.
     * @return false if both BasePair are equal.
     */
    bool operator!= (const BasePair &right) const
    {
      return ! operator== (right);
    }

    /**
     */
    bool operator< (const BasePair &right) const
    {
      return (&right != this
	      && (fResId < right.fResId
		  || (fResId == right.fResId
		      && rResId < right.rResId)));
    }
    
    // ACCESS ---------------------------------------------------------------

    // METHODS --------------------------------------------------------------

    void reverse ()
    {
      std::swap (first, second);
      std::swap (fResId, rResId);
    }

    // I/O  -----------------------------------------------------------------
    
  };
  
}

#endif
