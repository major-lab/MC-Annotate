#include "AnnotationInteractions.h"
#include "AnnotateModel.h"
#include "mccore/Pdbstream.h"
#include "mccore/GraphModel.h"

#include <cassert>
#include <iterator>
#include <sstream>


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
		mccore::GraphModel::edge_const_iterator eit;
		for (eit = aModel.edge_begin (); aModel.edge_end () != eit; ++eit)
		{
			const mccore::Residue *ref;
			const mccore::Residue *res;

			// Insure we're entering the relation only once
			if ((ref = (*eit)->getRef ())->getResId () < (res = (*eit)->getRes ())->getResId ())
			{
				const mccore::ResId refId = ref->getResId();
				const mccore::ResId resId = res->getResId();
				mccore::GraphModel::label refLabel = aModel.getVertexLabel (const_cast< mccore::Residue* > (ref));
				mccore::GraphModel::label resLabel = aModel.getVertexLabel (const_cast< mccore::Residue* > (res));

				if ((*eit)->isPairing ())
				{
					BasePair* pInteraction = newBasePair(aModel, *(*eit));
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
					BaseStack* pInteraction = newBaseStack(aModel, *(*eit));
					mStacks.push_back(*pInteraction);
					mInteractions.insert(pInteraction);
				}
			}
		}
		std::sort (mPairs.begin (), mPairs.end ());
    	std::sort (mStacks.begin (), mStacks.end ());
    	std::sort (mLinks.begin (), mLinks.end ());

    	mWWPairs = getWWBasePairs(aModel);
	}

	void AnnotationInteractions::update(AnnotateModel& aModel)
	{
		const mccore::GraphModel* pModel = &aModel;
		update(*pModel);
	}

	BaseStack* AnnotationInteractions::newBaseStack(
		const mccore::GraphModel& aModel,
		const mccore::Relation& aRelation) const
	{
		assert(aRelation.isStacking());
		const mccore::Residue *ref = aRelation.getRef ();
		const mccore::Residue *res = aRelation.getRes ();

		// Insure the relation is in the right order
		assert(ref->getResId () < res->getResId ());

		const mccore::ResId refId = ref->getResId();
		const mccore::ResId resId = res->getResId();
		mccore::GraphModel::label refLabel = aModel.getVertexLabel (const_cast< mccore::Residue* > (ref));
		mccore::GraphModel::label resLabel = aModel.getVertexLabel (const_cast< mccore::Residue* > (res));
		BaseStack::enAdjacency eAdjacency = BaseStack::eUnconnected;
		if(aRelation.is(PropertyType::pAdjacent5p))
		{
			eAdjacency = BaseStack::eAdjacent5p;
		}else if(aRelation.is(PropertyType::pAdjacent3p))
		{
			eAdjacency = BaseStack::eAdjacent3p;
		}
		BaseStack::enStackingType eStacking = BaseStack::eUpward;
		if(aRelation.is(PropertyType::pUpward))
		{
			eStacking = BaseStack::eUpward;
		}else if(aRelation.is(PropertyType::pDownward))
		{
			eStacking = BaseStack::eDownward;
		}else if(aRelation.is(PropertyType::pInward))
		{
			eStacking = BaseStack::eInward;
		}else if(aRelation.is(PropertyType::pOutward))
		{
			eStacking = BaseStack::eOutward;
		}
		BaseStack* pInteraction = new BaseStack(refLabel, refId, resLabel, resId, eAdjacency);
		return pInteraction;
	}

	BasePair* AnnotationInteractions::newBasePair(
			const mccore::GraphModel& aModel,
			const mccore::Relation& aRelation) const
	{
		assert(aRelation.isPairing());
		const mccore::Residue *ref = aRelation.getRef ();
		const mccore::Residue *res = aRelation.getRes ();

		// Insure the relation is in the right order
		assert(ref->getResId () < res->getResId ());

		const std::vector< std::pair< const mccore::PropertyType*, const mccore::PropertyType* > > &faces = aRelation.getPairedFaces ();

		const mccore::ResId refId = ref->getResId();
		const mccore::ResId resId = res->getResId();
		mccore::GraphModel::label refLabel = aModel.getVertexLabel (const_cast< mccore::Residue* > (ref));
		mccore::GraphModel::label resLabel = aModel.getVertexLabel (const_cast< mccore::Residue* > (res));
		BasePair* pInteraction = new BasePair(refLabel, refId, resLabel, resId, faces);

		if(aRelation.is(mccore::PropertyType::pCis))
		{
			pInteraction->orientation(BasePair::eCis);
		}else if(aRelation.is(mccore::PropertyType::pTrans))
		{
			pInteraction->orientation(BasePair::eTrans);
		}else
		{
			pInteraction->orientation(BasePair::eUnknown);
		}
		return pInteraction;
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
			std::pair<std::list< BaseStack >, std::list<BaseStack> > stacks;
			stacks = getSplitStacksAdjacencies();

			oss << "Adjacent stackings ----------------------------------------------" << endl;
			outputStacks(oss, stacks.first);

		    oss << "Non-Adjacent stackings ------------------------------------------" << endl;
		    outputStacks(oss, stacks.second);
	    	oss << "Number of stackings = " << mStacks.size () << endl
		    	<< "Number of adjacent stackings = " << stacks.first.size () - stacks.second.size () << endl
			    << "Number of non adjacent stackings = " << stacks.second.size () << endl;
		}
	}

	void AnnotationInteractions::outputStacks(
		std::ostringstream& oss,
		const std::list< BaseStack >& aStacks) const
	{
		std::list< BaseStack >::const_iterator it;
		for (it = aStacks.begin (); aStacks.end () != it; ++ it)
		{
			const mccore::Relation *rel = mpModel->internalGetEdge (it->first, it->second);
			const std::set< const mccore::PropertyType* > &labels = rel->getLabels ();
			oss << it->fResId << "-" << it->rResId << " : ";
			std::copy (labels.begin (), labels.end (), ostream_iterator< const mccore::PropertyType* > (oss, " "));
			oss << endl;
		}
	}

	std::pair<std::list< BaseStack >, std::list<BaseStack> >
	AnnotationInteractions::getSplitStacksAdjacencies() const
	{
		std::pair<std::list< BaseStack >, std::list<BaseStack> > stacks;
		std::vector<BaseStack>::const_iterator it;
		for (it = mStacks.begin (); mStacks.end () != it; ++ it)
		{
			const mccore::Relation *rel = mpModel->internalGetEdge (it->first, it->second);
			if (rel->is (PropertyType::pAdjacent))
			{
				stacks.first.push_back(*it);
			}
			else
			{
				stacks.second.push_back(*it);
			}
		}
		return stacks;
	}

	void AnnotationInteractions::outputPair (std::ostringstream& oss, const BasePair& aBasePair) const
	{
		if(NULL != mpModel)
		{
			const mccore::Relation &rel = *mpModel->internalGetEdge (aBasePair.first, aBasePair.second);
			const std::set< const mccore::PropertyType* > &labels = rel.getLabels ();
			std::vector< pair< const mccore::PropertyType*, const mccore::PropertyType* > >::const_iterator pfit;

			oss << aBasePair.fResId << '-' << aBasePair.rResId << " : ";
			oss << mccore::Pdbstream::stringifyResidueType (rel.getRef ()->getType())
				<< "-"
				<< mccore::Pdbstream::stringifyResidueType (rel.getRes ()->getType ())
				<< " ";
			for (pfit = aBasePair.faces().begin (); aBasePair.faces().end () != pfit; ++pfit)
			{
				oss << *pfit->first << "/" << *pfit->second << ' ';
			}
			copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (oss, " "));
		}
	}


	void AnnotationInteractions::outputPairs (std::ostringstream& oss) const
	{
    	std::vector< BasePair >::const_iterator bpit;

		for (bpit = mPairs.begin (); mPairs.end () != bpit; ++bpit)
		{
			outputPair(oss, *bpit);
			oss << std::endl;
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

	std::set< BasePair > AnnotationInteractions::getWWBasePairs(
		const mccore::GraphModel& aModel) const
	{
		std::set< BasePair> oWWBasePairs;
		std::set< BasePair> excludedPairs;

		std::vector<const BasePair*> resToPair;
		resToPair.resize(aModel.size());
		std::vector< BasePair >::const_iterator it;

		for(it = mPairs.begin(); mPairs.end() != it; ++it)
		{
			const mccore::Relation &rel = *aModel.internalGetEdge(it->first, it->second);

			// Filter on nucleotides
			if( isWatsonCrick(rel, false))
			{
				BasePair oWWPair = *it;
				if(oWWPair.rResId < oWWPair.fResId)
				{
					oWWPair.reverse();
				}
				oWWBasePairs.insert(oWWPair);

				// Check if residues share pairs
				if(NULL == resToPair[it->first] && NULL == resToPair[it->second])
				{
					resToPair[it->first] = &(*it);
					resToPair[it->second] = &(*it);
				}
				else
				{
					bool bRelationStrict = checkFacesStrict(rel);
					if(NULL != resToPair[it->first])
					{
						const BasePair* oldPair = resToPair[it->first];
						const mccore::Relation &oldRel = *aModel.internalGetEdge (oldPair->first, oldPair->second);
						if(!checkFacesStrict(oldRel) && bRelationStrict)
						{
							excludedPairs.insert(*oldPair);
							resToPair[it->first] = &(*it);
						}
						else
						{
							excludedPairs.insert(*it);
						}
					}

					if(NULL != resToPair[it->second])
					{
						const BasePair* oldPair = resToPair[it->second];
						const mccore::Relation &oldRel =
							*aModel.internalGetEdge (oldPair->first, oldPair->second);
						if(!checkFacesStrict(oldRel) && bRelationStrict)
						{
							excludedPairs.insert(*oldPair);
							resToPair[it->second] = &(*it);
						}
						else
						{
							excludedPairs.insert(*it);
						}
					}
				}
			}
		}
		return SetDifference<BasePair>(oWWBasePairs, excludedPairs);
	}

	bool AnnotationInteractions::checkFacesStrict(const mccore::Relation &aRelation) const
	{
		bool bFaces = false;
		const std::vector< pair< const PropertyType*, const PropertyType* > > &faces = aRelation.getPairedFaces ();
		std::vector< pair< const PropertyType*, const PropertyType* > >::const_iterator pfit;
		for (pfit = faces.begin (); faces.end () != pfit && !bFaces; ++pfit)
		{
			const PropertyType* pProp1 = pfit->first;
			const PropertyType* pProp2 = pfit->second;
			if(pProp1->toString() == "Ww" && pProp2->toString() == "Ww")
			{
				bFaces = true;
			}
		}
		return bFaces;
	}

	bool AnnotationInteractions::isWatsonCrick(
		const mccore::Relation &aRelation,
		bool bStrict) const
	{
		bool bIsWatsonCrick = checkNucleotides(aRelation)
			&& checkOrientation(aRelation);
		if(bIsWatsonCrick)
		{
			if(bStrict)
			{
				bIsWatsonCrick = checkFacesStrict(aRelation);
			}
			else
			{
				bIsWatsonCrick = checkFacesRelax(aRelation);
			}
		}

		return bIsWatsonCrick;
	}

	bool AnnotationInteractions::checkNucleotides(
		const mccore::Relation &aRelation) const
	{
		bool bNucleotides = false;
		const mccore::ResidueType* pRefType = aRelation.getRef()->getType();
		const mccore::ResidueType* pResType = aRelation.getRes()->getType();
		if((pRefType->isG() && (pResType->isC() || pResType->isU()))
			|| (pRefType->isA() && pResType->isU())
			|| (pRefType->isC() && pResType->isG())
			|| (pRefType->isU() && !pResType->isU()))
		{
			bNucleotides = true;
		}
		return bNucleotides;
	}

	bool AnnotationInteractions::checkFacesRelax(
		const mccore::Relation &aRelation) const
	{
		bool bFaces = false;
		const std::vector< pair< const PropertyType*, const PropertyType* > > &faces = aRelation.getPairedFaces ();
		std::vector< pair< const PropertyType*, const PropertyType* > >::const_iterator pfit;
		for (pfit = faces.begin (); faces.end () != pfit && !bFaces; ++pfit)
		{
			const PropertyType* pProp1 = pfit->first;
			const PropertyType* pProp2 = pfit->second;
			if(pProp1->isW() && pProp2->isW())
			{
				bFaces = true;
			}
		}
		return bFaces;
	}

	bool AnnotationInteractions::checkOrientation(
			const mccore::Relation &aRelation) const
	{
		// Filter on orientation
		bool bCis = false;
		bool bAntiparallel = false;
		const std::set< const PropertyType* > &labels = aRelation.getLabels ();
		std::set< const PropertyType* >::const_iterator it;
		for(it = labels.begin(); it != labels.end() && !(bCis && bAntiparallel); ++it)
		{
			if((*it)->toString() == "cis")
			{
				bCis = true;
			}
			else if((*it)->toString() == "antiparallel")
			{
				bAntiparallel = true;
			}
		}
		return bCis && bAntiparallel;
	}

	bool AnnotationInteractions::areContiguous(
			const mccore::ResId ref,
			const mccore::ResId res) const
	{
		bool bContiguous = false;
		std::vector<BaseLink>::const_iterator it;
		for(it = mLinks.begin(); it != mLinks.end() && !bContiguous; ++ it)
		{
			if(std::min(ref, res) == std::min(it->fResId, it->rResId)
				&& std::max(ref, res) == std::max(it->fResId, it->rResId))
			{
				bContiguous = true;
			}
		}
		return bContiguous;
	}
}
