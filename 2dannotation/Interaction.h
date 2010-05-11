#ifndef _Interaction_H_
#define _Interaction_H_

#include <string>
#include <utility>

#include <mccore/ResId.h>

namespace annotate {

class Interaction : public std::pair<mccore::ResId, mccore::ResId>
{
public:
	// LIFECYLCE ---------------------------------------------------------------
	Interaction(const mccore::ResId& aRes1, const mccore::ResId& aRes2);

	// OPERATORS ---------------------------------------------------------------
	bool operator <(const Interaction& aRight) const;
};

}; // namespace annotate

#endif /*_Interaction_H_*/
