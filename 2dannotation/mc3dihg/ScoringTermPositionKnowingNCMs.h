/*
 * ScoringTermPositionKnowingNCMs.h
 *
 *  Created on: Jun 2, 2010
 *      Author: blanchmf
 */

#ifndef _scoringtermpositionknowingncms_H_
#define _scoringtermpositionknowingncms_H_

#include "ScoringTerm.h"

#include <map>
#include <string>
#include <vector>

class ScoringTermPositionKnowingNCMs : public ScoringTerm
{
public:

	// LIFECYCLE ---------------------------------------------------------------
	ScoringTermPositionKnowingNCMs();


	virtual void read(const std::string& astrFile);
	virtual float score(
		const enProfile& aeProfile1,
		const enProfile& aeProfile2,
		unsigned int auiPos1,
		unsigned int auiPos2,
		unsigned int auiTypeId,
		const enNucleotide& aeNuc1,
		const enNucleotide& aeNuc2) const;
private:
	typedef std::pair<unsigned int, unsigned int> interaction_coordinate;
	std::map<std::string, unsigned int> mProfileMap;
	std::vector<std::vector<std::map<interaction_coordinate, float> > > mFrequencies;
};

#endif /* _scoringtermpositionknowingncms_H_ */
