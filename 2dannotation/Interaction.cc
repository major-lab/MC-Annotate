#include "Interaction.h"

namespace annotate {

Interaction::Interaction(const mccore::ResId& aRes1, const mccore::ResId& aRes2)
: std::pair<mccore::ResId, mccore::ResId>(aRes1, aRes2)
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
