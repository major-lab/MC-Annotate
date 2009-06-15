#include "AnnotationCycles.h"
#include "AnnotateModel.h"
#include "mccore/Molecule.h"
#include <sstream>

namespace annotate
{
	AnnotationCycles::AnnotationCycles()
	{
	}
	
	AnnotationCycles::~AnnotationCycles()
	{
		clear();
	}
	
	void AnnotationCycles::clear()
	{
		mCycles.clear();
	}
	
	const std::string AnnotationCycles::provides() const
	{
		std::string strAnnotationName = "Cycles";
		return strAnnotationName;
	}
		
	void AnnotationCycles::update(const AnnotateModel& aModel)
	{
		clear();
		
		if (! aModel.getCyclesMolecule().empty ())
		{
			Molecule::const_iterator itMol;
			ostringstream oss;
		  
			for(itMol = aModel.getCyclesMolecule().begin(); 
				itMol != aModel.getCyclesMolecule().end(); 
				++ itMol)
			{
				mCycles.push_back(*itMol);
			}
		}
	}
	
	std::string AnnotationCycles::output() const
	{
		ostringstream oss;
		int i = 0;
		std::vector<Cycle>::const_iterator itCycle;
		for (itCycle = mCycles.begin(); mCycles.end() != itCycle; ++itCycle)
	    {
	    	AbstractModel::const_iterator it;
	    	oss << "Cycle " << i << " : ";
	    	oss << "[";
			for (it = itCycle->begin (); itCycle->end () != it; ++it)
			{
				oss << " " << it->getResId () << *it->getType ();
			}
			oss << " ] " << itCycle->size () << endl;
			++ i;
	    }
	  	
		return oss.str();		
	}
	
	const std::vector< Cycle >& AnnotationCycles::getCycles() const
	{
		return mCycles;
	}
}