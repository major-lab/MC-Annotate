#ifndef _annotate_BaseInteraction_h_
#define _annotate_BaseInteraction_h_

#include <utility>

#include "mccore/GraphModel.h"
#include "mccore/ResId.h"

namespace annotate
{

	class BaseInteraction : public pair< mccore::GraphModel::label, mccore::GraphModel::label >
	{
	public:
		enum enType
		{
			ePAIR,
			eLINK,
			eSTACK,
			eUNKNOWN,
		};

	    mccore::ResId fResId;
	    mccore::ResId rResId;

    	// LIFECYCLE ------------------------------------------------------------

    	BaseInteraction (
    		mccore::GraphModel::label l,
    		const mccore::ResId &fResId,
    		mccore::GraphModel::label r,
    		const mccore::ResId &rResId)
		: pair< mccore::GraphModel::label, mccore::GraphModel::label > (l, r),
		fResId (fResId),
		rResId (rResId),
		meType (eUNKNOWN)
		{ }

    	~BaseInteraction () { }

    	// ACCESS --------------------------------------------------------------
    	enType& type() {return meType;}
    	const enType& type() const {return meType;}

    	// OPERATORS ------------------------------------------------------------

	    /**
	     */
	    virtual bool operator< (const BaseInteraction &right) const
	    {
	    	mccore::ResId thisMin = std::min(fResId, rResId);
	    	mccore::ResId rightMin = std::min(right.fResId, right.rResId);

	    	mccore::ResId thisMax = std::max(fResId, rResId);
	    	mccore::ResId rightMax = std::max(right.fResId, right.rResId);

			return (&right != this
			&& (thisMin < rightMin
				|| (thisMin == rightMin && (thisMax < rightMax
					|| (thisMax == rightMax && meType < right.meType)))));
	    }

	    bool sameResidues(const BaseInteraction &right) const
	    {
	    	return (&right == this
			      || (first == right.first && second == right.second));
	    }

	    protected:

			/**
		     * Compares this BaseInteraction to right.
		     * @param right the BaseInteraction to compare with this.
		     * @return true if both BaseInteraction are equal.
		     */
		    virtual bool operator== (const BaseInteraction &right) const
		    {
		      return sameResidues(right);
		    }

		    /**
		     * Compares this BaseInteraction to right.
		     * @param right the BaseInteraction to compare with this.
		     * @return false if both BaseInteraction are equal.
		     */
		    virtual bool operator!= (const BaseInteraction &right) const
		    {
		      return ! operator== (right);
		    }

		    enType meType;
    };
}

#endif /*_annotate_BaseInteraction_h_*/
