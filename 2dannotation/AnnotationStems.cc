#include "AlgorithmExtra.h"
#include "AnnotationStems.h"
#include "AnnotateModel.h"
#include "AnnotationInteractions.h"

#include <cassert>
#include <sstream>

namespace annotate
{
	// Static members
	std::string AnnotationStems::mstrAnnotationName = "Stems";

	// Methods
	AnnotationStems::AnnotationStems()
	{
		addRequirement<AnnotationInteractions>();
	}

	AnnotationStems::~AnnotationStems()
	{
		clear();
	}

	void AnnotationStems::update(AnnotateModel& aModel)
	{
		// Remove previous annotation if any
		clear();

		// Create all the potential stems
		std::vector<Stem> stems;
		getPotentialStems(aModel, stems);

		// Remove the stems of 1 pair only
		removeInvalids(stems);

		// Filter overlapping stems
		removeOverlaps(stems);

		// Keep the remaining stems for the annotation
		std::vector<Stem>::iterator itStem = stems.begin();
		for(itStem = stems.begin(); itStem != stems.end(); ++itStem)
		{
			// Name the stem properly
			std::string name = stemName(aModel, *itStem);
			itStem->name(name);
			mStems.push_back(*itStem);
		}
	}

	std::string AnnotationStems::stemName(
		const AnnotateModel& aModel,
		const Stem& aStem)
	{
		std::ostringstream oss;
    	mccore::ResId r1 = aStem.basePairs().front().fResId;
    	mccore::ResId r2 = aStem.basePairs().back().fResId;
    	mccore::ResId r3= aStem.basePairs().front().rResId;
    	mccore::ResId r4 = aStem.basePairs().back().rResId;
    	oss << r1 << "-" << r2 << ", " << r3 << "-" << r4;

		return oss.str();
	}

	std::string AnnotationStems::output() const
	{
		std::ostringstream oss;
    	vector< Stem >::const_iterator it;

    	for (it = mStems.begin (); mStems.end () != it; ++it)
	    {
	    	oss << it->name() << std::endl;
    	}
		return oss.str();
	}

	const std::vector< Stem >& AnnotationStems::getStems() const
	{
		return mStems;
	}

	void AnnotationStems::clear()
	{
		mStems.clear();
	}

	bool AnnotationStems::checkFacesStrict(const mccore::Relation &aRelation) const
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

	bool AnnotationStems::isBestPartner(
		const AnnotateModel& aModel,
		const BasePair& aPair,
		const std::set<BasePair>& aPairs) const
	{
		bool bBestPartner = true;
		const mccore::ResId resId = aPair.fResId;
		std::set<BasePair>::const_iterator it = aPairs.begin();
		// Find the first occurence
		while(it != aPairs.end() && (*it).fResId < resId)
		{
			++it;
		}

		// Check for a better partner
		while(it != aPairs.end() && (*it).fResId == resId)
		{
			if(aPair != *it)
			{
				const mccore::Relation &rel = *aModel.internalGetEdge (it->first, it->second);
				if(checkFacesStrict(rel))
				{
					bBestPartner = false;
				}
			}
			++ it;
		}
		return bBestPartner;
	}

	void AnnotationStems::getPotentialStems(
		const AnnotateModel& aModel,
		std::vector<Stem>& aStems) const
	{
		std::vector<Stem> potentialStems;
		std::set< BasePair > potentials;

		const AnnotationInteractions* pInteractions = NULL;
		pInteractions = aModel.getAnnotation<AnnotationInteractions>();
		assert(NULL != pInteractions);
		potentials = pInteractions->pairsWW();

		// Filter out multi-chain pairs
		potentials = filterOutMultiChainsPairs(potentials);

		// Create all the potential stems
		std::set< BasePair >::const_iterator it = potentials.begin();
		for( ; it != potentials.end(); ++ it)
		{
			bool bContinues = false;
			std::vector<Stem>::iterator itStem = potentialStems.begin();
			while(itStem != potentialStems.end())
			{
				bool bCont = itStem->continues(*it);

				if(bCont)
				{
					const mccore::Relation &rel = *aModel.internalGetEdge (it->first, it->second);
					if(!checkFacesStrict(rel))
					{
						bCont = isBestPartner(aModel, *it, potentials);
					}
				}

				if(bCont)
				{
					itStem->push_back(*it);
					bContinues = true;
					++itStem;
				}
				else
				{
					if(!itStem->basePairs().back().areContiguous(*it))
					{
						aStems.push_back(*itStem);
						itStem = potentialStems.erase(itStem);
					}
					else
					{
						++itStem;
					}
				}
			}

			if(!bContinues)
			{
				Stem currentStem;
				currentStem.push_back(*it);
				potentialStems.push_back(currentStem);
			}
		}

		for(std::vector<Stem>::const_iterator itS =	potentialStems.begin();
			itS != potentialStems.end();
			++itS)
		{
			aStems.push_back(*itS);
		}
	}

	std::set< BasePair > AnnotationStems::filterOutMultiChainsPairs(
		std::set< BasePair >& aPairs) const
	{
		std::set< BasePair > filteredSet;
		for(std::set<BasePair>::const_iterator it = aPairs.begin();
			aPairs.end() != it;
			++ it)
		{
			if(it->fResId.getChainId() == it->rResId.getChainId())
			{
				filteredSet.insert(*it);
			}
		}
		return filteredSet;
	}

	void AnnotationStems::removeInvalids(std::vector<Stem>& aStems) const
	{
		std::vector<Stem>::iterator itStem = aStems.begin();
		while(itStem != aStems.end())
		{
			if(2 > itStem->size())
			{
				itStem = aStems.erase(itStem);
			}
			else
			{
				++itStem;
			}
		}
	}

	void AnnotationStems::removeOverlaps(std::vector<Stem>& aStems) const
	{
		std::vector<Stem>::iterator itStem = aStems.begin();
		while(itStem != aStems.end())
		{
			if(2 > itStem->size())
			{
				itStem = aStems.erase(itStem);
			}
			else
			{
				++itStem;
			}
		}
	}
}
