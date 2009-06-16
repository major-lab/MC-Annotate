#include "AnnotationStems.h"
#include "AnnotateModel.h"
#include <sstream>

namespace annotate
{
	AnnotationStems::AnnotationStems()
	{
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
	
	std::set< BasePair > AnnotationStems::getWWBasePairs(
		const AnnotateModel& aModel) const
	{
	  	set< BasePair> oWWBasePairs;
	  	vector< BasePair >::const_iterator bpit = aModel.getBasePairs().begin();

		for(;aModel.getBasePairs().end () != bpit; ++bpit)
		{
			const mccore::Relation &rel = *aModel.internalGetEdge (bpit->first, bpit->second);
			
			// Filter on nucleotides
			if( checkNucleotides(rel) 
				&& checkFaces(rel) 
				&& checkOrientation(rel))
			{
				// Add the pair if it passed all tests
				BasePair oWWPair = *bpit;
				if(oWWPair.rResId < oWWPair.fResId)
				{
					oWWPair.reverse();	  				
				}
				oWWBasePairs.insert(oWWPair);
			}			
	  	}
  		return oWWBasePairs;
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
	
	bool AnnotationStems::checkFaces(const mccore::Relation &aRelation) const
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
	
	bool AnnotationStems::checkOrientation(const mccore::Relation &aRelation) const
	{
		// Filter on orientation
		bool bCis = false;
		const std::set< const PropertyType* > &labels = aRelation.getLabels ();
		std::set< const PropertyType* >::const_iterator it;
		for(it = labels.begin(); it != labels.end() && !bCis; ++it)
		{	  		
			if((*it)->toString() == "cis")
			{
				bCis = true;
			}
		}
		return bCis;
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
				if(itStem->continues(*it))
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