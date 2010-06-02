/*
 * ScoringTermTypeKnowingNCMs.h
 *
 *  Created on: Jun 2, 2010
 *      Author: blanchmf
 */

#ifndef _ScoringTermTypeKnowingNCMs_H_
#define _ScoringTermTypeKnowingNCMs_H_

#include "ScoringTerm.h"

#include "InteractionInfo.h"

#include <map>
#include <string>
#include <vector>

class ScoringTermTypeKnowingNCMs : public ScoringTerm
{
public:

	// LIFECYCLE ---------------------------------------------------------------
	ScoringTermTypeKnowingNCMs();

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
	std::map<std::string, enProfile> mProfileMap;

	typedef annotate::InteractionInfo::enFace face;
	typedef std::pair<face, face> face_pair;
	typedef std::pair<face_pair, annotate::BasePair::enOrientation> interaction_type;
	std::map<interaction_type, unsigned int> mInteractionTypeMap;
	std::vector<std::vector<std::vector<float> > > mFrequencies; // NCM x NCM x type

	unsigned int interactionTypeStringToId(const std::string& astrType) const;
};

#endif /* _ScoringTermTypeKnowingNCMs_H_ */
