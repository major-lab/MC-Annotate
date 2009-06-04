#include "Loop.h"

namespace annotate
{
	Loop::Loop() 
	{ 
	}
	
	Loop::Loop(const std::vector< Linker >& aLinkers)
	{
		mLinkers = aLinkers;
	}
	
	Loop::~Loop() 
	{
		clear();
	}
	
	const std::vector< Linker >& Loop::getLinkers() const
	{
		return mLinkers;
	}
	
	void Loop::clear()
	{
		mLinkers.clear();
	}
}