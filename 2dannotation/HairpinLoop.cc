#include "HairpinLoop.h"

namespace annotate
{
	HairpinLoop::HairpinLoop() 
	{ 
	}
	
	HairpinLoop::HairpinLoop(const std::vector<const Residue *>& aResidues)
	{
		mResidues = aResidues;
	}
	
	HairpinLoop::~HairpinLoop() 
	{ 
		mResidues.clear(); 
	}
	
	const std::vector< const Residue* >& HairpinLoop::getResidues() const
	{
		return mResidues;
	}
}