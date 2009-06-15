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
	
	std::set< BasePair > AnnotationStems::getWWBasePairs(const AnnotateModel& aModel) const
	{
	  	set< BasePair> oWWBasePairs;
	  	vector< BasePair >::const_iterator bpit;

		for(bpit = aModel.getBasePairs().begin (); 
			aModel.getBasePairs().end () != bpit; 
			++bpit)
		{
			const Relation &rel = *aModel.internalGetEdge (bpit->first, bpit->second);
			const vector< pair< const PropertyType*, const PropertyType* > > &faces = rel.getPairedFaces ();
			vector< pair< const PropertyType*, const PropertyType* > >::const_iterator pfit;

			for (pfit = faces.begin (); faces.end () != pfit; ++pfit)
		  	{
		  		const PropertyType* pProp1 = pfit->first;
	  			const PropertyType* pProp2 = pfit->second;
	  			std::string face1 = pProp1->toString();
		  		std::string face2 = pProp2->toString();
		  		if(	(0 < face1.size() && face1[0] == 'W') 
					&& (0 < face2.size() && face2[0] == 'W'))
		  		{
		  			BasePair oWWPair = *bpit;
	  				if(oWWPair.rResId < oWWPair.fResId)
	  				{
	  					oWWPair.reverse();	  				
		  			}
		  			oWWBasePairs.insert(oWWPair);
	  			}
			}
	  	}
  		return oWWBasePairs;
	}
	
	void AnnotationStems::getPotentialStems(
		const AnnotateModel& aModel, 
		std::vector<Stem>& aStems) const
	{
		std::set< BasePair > potentials;
		potentials = getWWBasePairs(aModel);
		
		// Create all the potential stems
		std::set< BasePair >::const_iterator it = potentials.begin();
		for( ; it != potentials.end(); ++ it)
		{
			bool bContinues = false;
			std::vector<Stem>::iterator itStem;
			for(itStem = aStems.begin(); itStem != aStems.end(); ++ itStem)
			{
				if(itStem->continues(*it))
				{
					itStem->push_back(*it);
					bContinues = true;
				}
			}
			
			if(!bContinues)
			{
				Stem currentStem;
				currentStem.push_back(*it);
				aStems.push_back(currentStem);
			}
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