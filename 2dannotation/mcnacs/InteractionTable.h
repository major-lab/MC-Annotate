/*
 * InteractionTable.h
 *
 *  Created on: Jan 18, 2010
 *      Author: blanchmf
 */

#ifndef INTERACTIONTABLE_H_
#define INTERACTIONTABLE_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace annotate {
	class CycleInfo;
};

class InteractionTable
{
public:
	typedef std::pair<annotate::CycleInfo, annotate::CycleInfo> interacting_pair;
	typedef std::set<interacting_pair> interacting_set;

	// LIFECYCLE ---------------------------------------------------------------
	InteractionTable();
	~InteractionTable();

	void addInteraction(const annotate::CycleInfo& aRow, const annotate::CycleInfo& aColumn);

	const interacting_set& operator()(
		const std::string& astrRow,
		const std::string& astrColumn) const;

	// METHODS -----------------------------------------------------------------
	std::string toString() const;
	std::list<std::string> getProfileList() const;
private:
	std::map<std::string, unsigned int> mIndices;
	std::vector<std::vector<interacting_set> > mTable;
};


#endif /* INTERACTIONTABLE_H_ */
