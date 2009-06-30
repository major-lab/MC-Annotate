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
		rResId (rResId)
		{ }

    	~BaseInteraction () { }
    	
    	// OPERATORS ------------------------------------------------------------
	
	    /**
	     */
	    virtual bool operator< (const BaseInteraction &right) const
	    {
	      return (&right != this
		      && (fResId < right.fResId
			  || (fResId == right.fResId
			      && rResId < right.rResId)));
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
    };
}

#endif /*_annotate_BaseInteraction_h_*/
