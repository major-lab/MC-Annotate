#include "AlgorithmExtra.h"
#include "AnnotationStems.h"
#include "AnnotateModel.h"
#include "AnnotationInteractions.h"
#include <sstream>

namespace annotate
{
	AnnotationStems::AnnotationStems()
	{
		addRequirement(AnnotationInteractions().provides());
	}
	
	AnnotationStems::~AnnotationStems()
	{
		clear();
	}
	
	const std::string AnnotationStems::provides() const
	{
		std::string strAnnotationName = "Stems";
		return strAnnotationName;
	}
		
	void AnnotationStems::update(const AnnotateModel& aModel)
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
			mStems.push_back(*itStem);
		}
	}	
	
	std::string AnnotationStems::output() const
	{
		std::ostringstream oss;
		int i=0;
    	vector< Stem >::const_iterator it;	
    
    	for (it = mStems.begin (); mStems.end () != it; ++it)
	    {
	    	ResId r1, r2, r3, r4;
	    	r1 = (*it).basePairs().front().fResId;
	    	r2 = (*it).basePairs().back().fResId;
	    	r3= (*it).basePairs().front().rResId;
	    	r4 = (*it).basePairs().back().rResId;
	    	oss << "Stem " << i << " : " << r1 << "-" << r2;
	    	oss << ", " << r3 << "-" << r4 << std::endl; 
	    	i ++;
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
	
	bool AnnotationStems::isWatsonCrick(
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
	
	std::set< BasePair > AnnotationStems::getWWBasePairs(
		const AnnotateModel& aModel) const
	{
		std::set< BasePair> oWWBasePairs;
		std::set< BasePair> excludedPairs;
		const AnnotationInteractions* pAInteractions = NULL;
		pAInteractions = aModel.getAnnotation<AnnotationInteractions>(AnnotationInteractions().provides());
		
		if(NULL != pAInteractions)
		{
			std::vector<const BasePair*> resToPair;
			resToPair.resize(aModel.size());
		  	std::vector< BasePair >::const_iterator bpit = pAInteractions->getPairs().begin();
	
			for(;pAInteractions->getPairs().end () != bpit; ++bpit)
			{
				const mccore::Relation &rel = *aModel.internalGetEdge (bpit->first, bpit->second);
				
				// Filter on nucleotides
				if( isWatsonCrick(rel, false))
				{
					BasePair oWWPair = *bpit;
		  			if(oWWPair.rResId < oWWPair.fResId)
					{
						oWWPair.reverse();	  				
					}
					oWWBasePairs.insert(oWWPair);
					
					// Check if residues share pairs
					if(NULL == resToPair[bpit->first] && NULL == resToPair[bpit->second])
					{
						resToPair[bpit->first] = &(*bpit);
						resToPair[bpit->second] = &(*bpit);
					}
					else 
					{
						bool bRelationStrict = checkFacesStrict(rel);
						if(NULL != resToPair[bpit->first])
						{
							const BasePair* oldPair = resToPair[bpit->first];
							const mccore::Relation &oldRel = *aModel.internalGetEdge (oldPair->first, oldPair->second);
							if(!checkFacesStrict(oldRel) && bRelationStrict)
							{
								excludedPairs.insert(*oldPair);
								resToPair[bpit->first] = &(*bpit);
							} 
							else
							{
								excludedPairs.insert(*bpit);
							}
						}
						
						if(NULL != resToPair[bpit->second])
						{
							const BasePair* oldPair = resToPair[bpit->second];
							const mccore::Relation &oldRel = 
								*aModel.internalGetEdge (oldPair->first, oldPair->second);
							if(!checkFacesStrict(oldRel) && bRelationStrict)
							{
								excludedPairs.insert(*oldPair);
								resToPair[bpit->second] = &(*bpit);
							}
							else
							{
								excludedPairs.insert(*bpit);
							}
						}
					}
				}
			}		
	  	}
  		return SetDifference<BasePair>(oWWBasePairs, excludedPairs);
	}
	
	bool AnnotationStems::checkNucleotides(
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
	
	bool AnnotationStems::checkFacesRelax(const mccore::Relation &aRelation) const
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
	
	bool AnnotationStems::checkOrientation(const mccore::Relation &aRelation) const
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
					gOut (0) << it->fResId << "-" << it->rResId << " is considered a best partner" << std::endl;
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
		potentials = getWWBasePairs(aModel);
		
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