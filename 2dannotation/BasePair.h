//                              -*- Mode: C++ -*-
// BasePair.h
// Copyright Â© 2006 Laboratoire de Biologie Informatique et ThÃ©orique
//                  UniversitÃ© de MontrÃ©al.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Wed Jul 12 16:27:16 2006
// $Revision: 59 $
// $Id: BasePair.h 59 2006-11-15 21:25:50Z larosem $
//


#ifndef _annotate_BasePair_h_
#define _annotate_BasePair_h_

#include <algorithm>
#include "BaseInteraction.h"

namespace annotate
{

  class BasePair : public BaseInteraction
  {

  public:
	  typedef std::vector< std::pair< const mccore::PropertyType*, const mccore::PropertyType* > > face_vector;
	  enum enOrientation
	  {
		  eCis,
		  eTrans,
		  eUnknown
	  };

    // LIFECYCLE ------------------------------------------------------------

    BasePair (	mccore::GraphModel::label l,
    			const mccore::ResId &fResId,
    			mccore::GraphModel::label r,
    			const mccore::ResId &rResId,
    			const face_vector& aFaces)
      : BaseInteraction(l, fResId, r, rResId),
      mFaces(aFaces),
      meOrientation(eUnknown)
    {
    	meType = BaseInteraction::ePAIR;
    }

    BasePair (const BasePair& aPair)
	: BaseInteraction(aPair.first, aPair.fResId, aPair.second, aPair.rResId),
		mFaces(aPair.mFaces),
		meOrientation(aPair.meOrientation)
    {
    	meType = BaseInteraction::ePAIR;
    }

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
    	  mFaces = right.mFaces;
    	  meOrientation = right.meOrientation;
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
      return BaseInteraction::operator == (right);
    }

    /**
     * Compares this BasePair to right.
     * @param right the BasePair to compare with this.
     * @return false if both BasePair are equal.
     */
    bool operator!= (const BasePair &right) const
    {
		return BaseInteraction::operator != (right);
    }

    /**
     */
    bool operator< (const BasePair &right) const
    {
      return BaseInteraction::operator < (right);
    }

    // ACCESS ---------------------------------------------------------------
    const face_vector& faces() const {return mFaces;}
    const enOrientation& orientation() const {return meOrientation;}
    void orientation(const enOrientation& aOrient) {meOrientation = aOrient;}

    // METHODS --------------------------------------------------------------
	bool areContiguous(const BasePair &aBasePair) const
	{
		bool bContiguous = fResId.areContiguous(aBasePair.fResId);
		bContiguous = bContiguous && rResId.areContiguous(aBasePair.rResId);
		return bContiguous;
	}

    void reverse ()
    {
      std::swap (first, second);
      std::swap (fResId, rResId);
      mFaces = reverseFaces(mFaces);
    }

    static face_vector reverseFaces(const face_vector& aFaces);



    // I/O  -----------------------------------------------------------------

	private:
		std::vector< std::pair< const mccore::PropertyType*, const mccore::PropertyType* > > mFaces;
		enOrientation meOrientation;
	};
}

#endif
