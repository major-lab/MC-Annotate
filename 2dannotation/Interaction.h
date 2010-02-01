#ifndef _Interaction_H_
#define _Interaction_H_

#include <string>
#include <utility>

namespace annotate {

class Interaction : public std::pair<std::string, std::string>
{
public:
	// LIFECYLCE ---------------------------------------------------------------
	Interaction(const std::string& aRes1, const std::string& aRes2);

	// OPERATORS ---------------------------------------------------------------
	bool operator <(const Interaction& aRight) const;
};

}; // namespace annotate

#endif /*_Interaction_H_*/
