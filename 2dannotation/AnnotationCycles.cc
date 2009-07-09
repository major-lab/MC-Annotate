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
					Cycle cycle(*itMol, aModel.relationMask());
					std::string name = cycleName(aModel, cycle);
					cycle.name(name);
					mCycles.push_back(cycle);
				}				
			}
		}
	}
	
	std::string AnnotationCycles::cycleName(
		const AnnotateModel& aModel, 
		const Cycle& aCycle) const
	{
		std::ostringstream oss;
		oss << "[" << aModel.name() << "] ";
		
		std::list<mccore::ResId>::const_iterator it;
		for(it = aCycle.getResidues().begin();
			it != aCycle.getResidues().end(); 
			++ it)
		{
			oss << *it << " ";
		}		
		return oss.str();
	}
	
	std::string AnnotationCycles::output() const
	{
		ostringstream oss;
		std::list<Cycle>::const_iterator itCycle;
		for (itCycle = mCycles.begin(); mCycles.end() != itCycle; ++itCycle)
	    {
	    	oss << itCycle->name() << std::endl;
	    }
	  	
		return oss.str();		
	}
	
	const std::list< Cycle >& AnnotationCycles::getCycles() const
	{
		return mCycles;
	}
}