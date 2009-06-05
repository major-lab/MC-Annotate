#include "Annotation.h"

namespace annotate
{
	const std::set<std::string >& Annotation::requires() const
	{
		return mRequirements;
	}
	
	void Annotation::addRequirement(const std::string& aRequirement)
	{
		mRequirements.insert(aRequirement);
	}
}