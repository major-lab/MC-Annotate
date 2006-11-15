//                              -*- Mode: C++ -*- 
// BaseLink.h
// Copyright Â© 2006 Laboratoire de Biologie Informatique et ThÃ©orique
//                  UniversitÃ© de MontrÃ©al.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Wed Aug  2 18:34:59 2006
// $Revision$
// $Id$
// 


#ifndef _annotate_BaseLink_h_
#define _annotate_BaseLink_h_

#include <utility>

#include "mccore/GraphModel.h"
#include "mccore/ResId.h"

using namespace mccore;
using namespace std;



namespace annotate
{

  class BaseLink : public pair< GraphModel::label, GraphModel::label >
  {

  public:

    ResId fResId;

    ResId rResId;
    
    // LIFECYCLE ------------------------------------------------------------

    BaseLink (GraphModel::label l, const ResId &fResId, GraphModel::label r, const ResId &rResId)
      : pair< GraphModel::label, GraphModel::label > (l, r),
	fResId (fResId),
	rResId (rResId)
    { }

    ~BaseLink () { }
    
    // OPERATORS ------------------------------------------------------------

    /**
     * Assigns the right BaseLink to this.
     * @param right the BaseLink to copy into this.
     * @return a reference to this.
     */
    BaseLink& operator= (const BaseLink &right)
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
     * Compares this BaseLink to right.
     * @param right the BaseLink to compare with this.
     * @return true if both BaseLink are equal.
     */
    bool operator== (const BaseLink &right) const
    {
      return (&right == this
	      || (first == right.first && second == right.second));
    }

    /**
     * Compares this BaseLink to right.
     * @param right the BaseLink to compare with this.
     * @return false if both BaseLink are equal.
     */
    bool operator!= (const BaseLink &right) const
    {
      return ! operator== (right);
    }

    /**
     */
    bool operator< (const BaseLink &right) const
    {
      return (&right != this
	      && (fResId < right.fResId
		  || (fResId == right.fResId
		      && rResId < right.rResId)));
    }
    
    // ACCESS ---------------------------------------------------------------

    // METHODS --------------------------------------------------------------

    // I/O  -----------------------------------------------------------------
    
  };
  
}

#endif
