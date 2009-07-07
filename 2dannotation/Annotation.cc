#include "Annotation.h"

namespace annotate
{
	const std::set<std::string >& Annotation::requires() const
	{
		return mRequirements;
	}
}