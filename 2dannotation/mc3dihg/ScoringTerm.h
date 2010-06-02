/*
 * scoringterm.h
 *
 *  Created on: Jun 2, 2010
 *      Author: blanchmf
 */

#ifndef _scoringterm_H_
#define _scoringterm_H_

#include <string>

class ScoringTerm
{
public:
	enum enNucleotide
	{
		eNucleotideA,
		eNucleotideC,
		eNucleotideG,
		eNucleotideU
	};

	enum enProfile
	{
		eProfile3 = 0,
		eProfile4 = 1,
		eProfile5 = 2,
		eProfile6 = 3,
		eProfile2_2 = 4,
		eProfile2_3 = 5,
		eProfile2_4 = 6,
		eProfile2_5 = 7,
		eProfile2_6 = 8,
		eProfile3_3 = 9,
		eProfile3_4 = 10,
		eProfile3_5 = 11,
		eProfile4_4 = 12
	};

	virtual void read(const std::string& aFile) = 0;
	virtual float score(
		const enProfile& aeProfile1,
		const enProfile& aeProfile2,
		unsigned int auiPos1,
		unsigned int auiPos2,
		unsigned int auiTypeId,
		const enNucleotide& aeNuc1,
		const enNucleotide& aeNuc2) const = 0;
};

#endif /* _scoringterm_H_ */
