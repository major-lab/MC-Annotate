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

    	mWWPairs = getWWBasePairs(aModel);
	}

	void AnnotationInteractions::update(AnnotateModel& aModel)
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
		std::list<const BaseInteraction*> interactions =
			getInteractions(ref, res);

		std::list<const BaseInteraction*>::const_iterator it;
		for(it = interactions.begin();
			it != interactions.end() && !bContiguous;
			++it)
		{
			const BaseLink* pLink = dynamic_cast<const BaseLink*>(*it);
			bContiguous = (NULL != pLink);
		}
		return bContiguous;
	}
}
