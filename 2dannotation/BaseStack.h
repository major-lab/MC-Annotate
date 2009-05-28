//                              -*- Mode: C++ -*- 
// BaseStack.h
// Copyright Â© 2006 Laboratoire de Biologie Informatique et ThÃ©orique
//                  UniversitÃ© de MontrÃ©al.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Wed Aug  2 18:11:28 2006
// $Revision: 59 $
// $Id: BaseStack.h 59 2006-11-15 21:25:50Z larosem $
// 


#ifndef _annotate_BaseStack_h_
#define _annotate_BaseStack_h_

#include <utility>

#include "mccore/GraphModel.h"
#include "mccore/ResId.h"

using namespace mccore;
using namespace std;



namespace annotate
{

  class BaseStack : public pair< GraphModel::label, GraphModel::label >
  {
    
  public:
    
    ResId fResId;

    ResId rResId;

    // LIFECYCLE ------------------------------------------------------------

    BaseStack (GraphModel::label l, const ResId &fResId, GraphModel::label r, const ResId &rResId)
      : pair< GraphModel::label, GraphModel::label > (l, r),
	fResId (fResId),
	rResId (rResId)
    { }

    ~BaseStack () { }
    
    // OPERATORS ------------------------------------------------------------

    /**
     * Assigns the right BaseStack to this.
     * @param right the BaseStack to copy into this.
     * @return a reference to this.
     */
    BaseStack& operator= (const BaseStack &right)
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
     * Compares this BaseStack to right.
     * @param right the BaseStack to compare with this.
     * @return true if both BaseStack are equal.
     */
    bool operator== (const BaseStack &right) const
    {
      return (&right == this
	      || (first == right.first && second == right.second));
    }

    /**
     * Compares this BaseStack to right.
     * @param right the BaseStack to compare with this.
     * @return false if both BaseStack are equal.
     */
    bool operator!= (const BaseStack &right) const
    {
      return ! operator== (right);
    }

    /**
     */
    bool operator< (const BaseStack &right) const
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
