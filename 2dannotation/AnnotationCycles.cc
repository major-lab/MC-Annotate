#include "AnnotationCycles.h"
#include "AnnotateModel.h"
#include "mccore/Molecule.h"
#include <sstream>

namespace annotate
{
	std::string AnnotationCycles::mstrAnnotationName = "Cycles";
	
	AnnotationCycles::AnnotationCycles(unsigned int auiMaxCycleSize)
	{
		muiMaxCycleSize = auiMaxCycleSize;
	}
	
	AnnotationCycles::~AnnotationCycles()
	{
		clear();
	}
	
	void AnnotationCycles::clear()
	{
		mCycles.clear();
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
				if(0 == muiMaxCycleSize || itMol->size() <= muiMaxCycleSize)
				{
					Cycle cycle(*itMol);
					mCycles.push_back(cycle);
				}				
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
			for (it = itCycle->getModel().begin (); itCycle->getModel().end () != it; ++it)
			{
				oss << " " << it->getResId () << *it->getType ();
			}
			oss << " ] " << itCycle->getModel().size () << endl;
			++ i;
	    }
	  	
		return oss.str();		
	}
	
	const std::vector< Cycle >& AnnotationCycles::getCycles() const
	{
		return mCycles;
	}
}