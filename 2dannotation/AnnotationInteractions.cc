#include "AnnotationInteractions.h"
#include "AnnotateModel.h"
#include "mccore/Pdbstream.h"
#include "mccore/GraphModel.h"

#include <sstream>
#include <iterator>

namespace annotate
{
	static const unsigned int PAIRING_MARK = 1;
	
	std::string AnnotationInteractions::mstrAnnotationName = "Interactions";
		
	AnnotationInteractions::AnnotationInteractions() : mpModel(NULL)
	{}
	
	AnnotationInteractions::~AnnotationInteractions()
	{
		clear();
	}
	
	void AnnotationInteractions::clear()
	{
		mPairs.clear();
    	mStacks.clear();
		mLinks.clear();
    	mMarks.clear();
    	
    	std::multiset< BaseInteraction*, less_ptr<BaseInteraction> >::const_iterator it;
    	for(it = mInteractions.begin(); it != mInteractions.end(); ++ it)
    	{
    		*it;
    	}
    	mInteractions.clear();
	}
	
	void AnnotationInteractions::update(const mccore::GraphModel& aModel)
	{
		mpModel = &aModel;
		
		clear();
		mMarks.resize (aModel.size (), 0);
		mccore::GraphModel::edge_const_iterator eit;
		for (eit = aModel.edge_begin (); aModel.edge_end () != eit; ++eit)
		{
			const mccore::Residue *ref;
			const mccore::Residue *res;

			if ((ref = (*eit)->getRef ())->getResId () < (res = (*eit)->getRes ())->getResId ())
			{
				const mccore::ResId refId = ref->getResId();
				const mccore::ResId resId = res->getResId();
				mccore::GraphModel::label refLabel = aModel.getVertexLabel (const_cast< mccore::Residue* > (ref));
				mccore::GraphModel::label resLabel = aModel.getVertexLabel (const_cast< mccore::Residue* > (res));
				
				if ((*eit)->isPairing ())
				{
					mMarks[refLabel] |= PAIRING_MARK;
					mMarks[resLabel] |= PAIRING_MARK;
					BasePair* pInteraction = new BasePair(refLabel, refId, resLabel, resId);
					mPairs.push_back(*pInteraction);
					mInteractions.insert(pInteraction);
				}
				if ((*eit)->is (PropertyType::pAdjacent))
				{
					BaseLink* pInteraction = new BaseLink(refLabel, refId, resLabel, resId);
					mLinks.push_back(*pInteraction);
					mInteractions.insert(pInteraction);
				}
				if ((*eit)->isStacking ())
				{
					BaseStack* pInteraction = new BaseStack(refLabel, refId, resLabel, resId);
					mStacks.push_back(*pInteraction);
					mInteractions.insert(pInteraction);
				}
			}
		}
		std::sort (mPairs.begin (), mPairs.end ());
    	std::sort (mStacks.begin (), mStacks.end ());
    	std::sort (mLinks.begin (), mLinks.end ());
	}
	
	void AnnotationInteractions::update(const AnnotateModel& aModel)
	{
		const mccore::GraphModel* pModel = &aModel;
		update(*pModel);
	}
	
	std::string AnnotationInteractions::output() const
	{
		std::ostringstream oss;
		outputStacks (oss);
    	oss << "Base-pairs ------------------------------------------------------" << endl;
    	outputPairs (oss);
		return oss.str();
	}
	
	void AnnotationInteractions::outputStacks (std::ostringstream& oss) const
	{
		if(NULL != mpModel)
		{
			std::vector< BaseStack > nonAdjacentStacks;
			std::vector< BaseStack >::const_iterator bsit;

			oss << "Adjacent stackings ----------------------------------------------" << endl;

			for (bsit = mStacks.begin (); mStacks.end () != bsit; ++bsit)
			{
				const mccore::Relation *rel = mpModel->internalGetEdge (bsit->first, bsit->second);
				if (rel->is (PropertyType::pAdjacent))
				{
					const std::set< const mccore::PropertyType* > &labels = rel->getLabels ();
	
					oss << bsit->fResId << "-" << bsit->rResId << " : ";
		    		std::copy (labels.begin (), labels.end (), ostream_iterator< const mccore::PropertyType* > (oss, " "));
					oss << endl;
				}
				else
				{
					nonAdjacentStacks.push_back (*bsit);
				}
			}
    
		    oss << "Non-Adjacent stackings ------------------------------------------" << endl;
	    
			for (bsit = nonAdjacentStacks.begin (); nonAdjacentStacks.end () != bsit; ++bsit)
			{
				const std::set< const mccore::PropertyType* > &labels = mpModel->internalGetEdge (bsit->first, bsit->second)->getLabels ();
			
				oss << bsit->fResId << "-" << bsit->rResId << " : ";
				std::copy (labels.begin (), labels.end (), ostream_iterator< const mccore::PropertyType* > (oss, " "));
				oss << endl;
			}
	
	    	oss << "Number of stackings = " << mStacks.size () << endl
		    	<< "Number of adjacent stackings = " << mStacks.size () - nonAdjacentStacks.size () << endl
			    << "Number of non adjacent stackings = " << nonAdjacentStacks.size () << endl;
		}
	}
  
	void AnnotationInteractions::outputPair (std::ostringstream& oss, const BasePair& aBasePair) const
	{
		if(NULL != mpModel)
		{
			const mccore::Relation &rel = *mpModel->internalGetEdge (aBasePair.first, aBasePair.second);
			const std::set< const mccore::PropertyType* > &labels = rel.getLabels ();
			const std::vector< std::pair< const mccore::PropertyType*, const mccore::PropertyType* > > &faces = rel.getPairedFaces ();
			std::vector< pair< const mccore::PropertyType*, const mccore::PropertyType* > >::const_iterator pfit;
	
			oss << aBasePair.fResId << '-' << aBasePair.rResId << " : ";
			oss << mccore::Pdbstream::stringifyResidueType (rel.getRef ()->getType())
				<< "-"
				<< mccore::Pdbstream::stringifyResidueType (rel.getRes ()->getType ())
				<< " ";
			for (pfit = faces.begin (); faces.end () != pfit; ++pfit)
			{
				oss << *pfit->first << "/" << *pfit->second << ' ';
			}
			copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (oss, " "));
			oss << endl;
		}
	}
  

	void AnnotationInteractions::outputPairs (std::ostringstream& oss) const
	{
    	std::vector< BasePair >::const_iterator bpit;

		for (bpit = mPairs.begin (); mPairs.end () != bpit; ++bpit)
		{
			outputPair(oss, *bpit);
		}
	}
	
	std::list<const BaseInteraction*> AnnotationInteractions::getInteractions(
			const mccore::ResId ref, 
			const mccore::ResId res) const
	{
		std::list<const BaseInteraction*> interactions;
		std::multiset< BaseInteraction*, less_ptr<BaseInteraction> >::const_iterator it;
		BaseInteraction key(0, std::min(ref, res), 0, std::max(ref, res)); // Dummy interaction for search
		// TODO : Return a list of interactions
		std::pair<interaction_multiset_const_it, interaction_multiset_const_it> range;
		range = mInteractions.equal_range(&key);
		for(it = range.first; it != range.second; ++it)
		{
			interactions.push_back(*it);
		}
		return interactions;
	}
	
	std::list<const BaseInteraction*> AnnotationInteractions::getInteractions(
			const std::set<mccore::ResId>& aResIds) const
	{
		std::list<const BaseInteraction*> interactions;
		std::multiset< BaseInteraction*, less_ptr<BaseInteraction> >::const_iterator it;
		it = mInteractions.begin();
		for(it = mInteractions.begin(); it != mInteractions.end(); ++ it)
		{
			if(	(aResIds.find((*it)->fResId) != aResIds.end()) 
				&& (aResIds.find((*it)->rResId) != aResIds.end()))
			{
				interactions.push_back(*it);
			}
		}
		return interactions;
	}
}