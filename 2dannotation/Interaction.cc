#include "Interaction.h"

namespace annotate {

Interaction::Interaction(const std::string& aRes1, const std::string& aRes2)
: std::pair<std::string, std::string>(aRes1, aRes2)
{
}

	// OPERATORS ---------------------------------------------------------------
bool Interaction::operator <(const Interaction& aRight) const
{
	return ((std::min(first, second) < std::min(aRight.first, aRight.second))
	|| (std::min(first, second) == std::min(aRight.first, aRight.second)
		&& (std::max(first, second) < std::max(aRight.first, aRight.second))));
}

};
