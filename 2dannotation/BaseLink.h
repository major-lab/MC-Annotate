//                              -*- Mode: C++ -*-
// BaseLink.h
// Copyright Â© 2006 Laboratoire de Biologie Informatique et ThÃ©orique
//                  UniversitÃ© de MontrÃ©al.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Wed Aug  2 18:34:59 2006
// $Revision: 59 $
// $Id: BaseLink.h 59 2006-11-15 21:25:50Z larosem $
//


#ifndef _annotate_BaseLink_h_
#define _annotate_BaseLink_h_

#include "BaseInteraction.h"

namespace annotate
{

  class BaseLink : public BaseInteraction
  {

  public:

    // LIFECYCLE ------------------------------------------------------------

    BaseLink (	mccore::GraphModel::label l,
    			const mccore::ResId &fResId,
    			mccore::GraphModel::label r,
    			const mccore::ResId &rResId)
      : BaseInteraction(l, fResId, r, rResId)
    {
    	meType = BaseInteraction::eLINK;
    }

    BaseLink ( const BaseLink& aLink )
    : BaseInteraction(aLink.first, aLink.fResId, aLink.second, aLink.rResId)
    {
    	meType = BaseInteraction::eLINK;
    }

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
      return BaseInteraction::operator == (right);
    }

    /**
     * Compares this BaseLink to right.
     * @param right the BaseLink to compare with this.
     * @return false if both BaseLink are equal.
     */
    bool operator!= (const BaseLink &right) const
    {
		return BaseInteraction::operator != (right);
    }

    /**
     */
    bool operator< (const BaseLink &right) const
    {
		return BaseInteraction::operator<(right);
    }

    // ACCESS ---------------------------------------------------------------

    // METHODS --------------------------------------------------------------

    // I/O  -----------------------------------------------------------------

  };

}

#endif
