#include "TertiaryStructure.h"

#include "AnnotateModel.h"

namespace annotate
{
	TertiaryStructure::TertiaryStructure()
	{
	}
	
	TertiaryStructure::~TertiaryStructure()
	{
	}
	
	void TertiaryStructure::update(const AnnotateModel& aModel)
	{
		mModel.clear();
		std::set<mccore::ResId>::const_iterator itRes;
		for(itRes = mResidues.begin(); itRes != mResidues.end(); ++ itRes)
		{
			GraphModel::const_iterator itResidue = aModel.find(*itRes);
			mModel.insert(*itResidue);
		}
	}
	
	void TertiaryStructure::addCycle(const Cycle& aCycle) 
	{
		mCycles.push_back(aCycle);
		mResidues.insert(aCycle.getResidues().begin(), aCycle.getResidues().end());
	}
}