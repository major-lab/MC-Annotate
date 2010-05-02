/*
 * t3tablebuilder.h
 *
 *  Created on: May 1, 2010
 *      Author: blanchmf
 */

#ifndef _t3tablebuilder_H_
#define _t3tablebuilder_H_

#include "CycleInfo.h"
#include "RNAGroupFile.h"
#include "RNAGroupModel.h"
#include <map>
#include <set>
#include <string>
#include <vector>

class T3TableBuilder
{
public:
	// LIFECYCLE
	T3TableBuilder(const std::string& astrGroupFile, const std::string& aInteractions);
	void computeTable();

protected:
	typedef std::pair<unsigned int, unsigned int> interaction_coordinate;
	typedef annotate::InteractionInfo::enFace face;
	typedef std::pair<face, face> face_pair;
	typedef std::pair<face_pair, annotate::BasePair::enOrientation> interaction_type;
	typedef annotate::tuple3<annotate::InteractionInfo, annotate::CycleInfo, annotate::CycleInfo> interaction_cycle_pair;

	std::vector<std::vector<std::vector<std::vector<std::vector<float> > > > > mFrequencies;
	std::map<std::string, unsigned int> mProfileMap;
	std::map<interaction_type, unsigned int> mInteractionTypeMap;
	std::map<unsigned int, annotate::RNAGroupModel> mGroupModelMap;
	annotate::RNAGroupFile mGroupFile;
	std::set<interaction_cycle_pair> mInteractions;

	unsigned int profileId(const std::string& astrProfile) const;
	virtual void computeStatistics() = 0;
	unsigned int profileId(const annotate::CycleInfo& aCycle) const;
	unsigned int interactionTypeId(const annotate::InteractionInfo& aInteraction) const;
	unsigned int flipInteractionTypeId(const annotate::InteractionInfo& aInteraction) const;
	std::string profileString(unsigned int auiProfile) const;
	std::string interactionTypeString(unsigned int auiTypeId) const;

	interaction_coordinate getInteractionCoordinates(const interaction_cycle_pair& aInteraction) const;
	unsigned int getNucleotidePosition(const annotate::CycleInfo& aCycle, const mccore::ResId& aResId) const;
	unsigned int profileSize(unsigned int auiProfileId) const;
private:
	void readInteractions(const std::string& astrInteractions);
	void readPDBGroups(const std::string& astrGroups);

	annotate::CycleInfo getCycleInfo(
		const std::string& astrPDB, unsigned int uiModel,
		const std::string& astrProfile,
		const std::string& astrResIds,
		const std::string& astrSequence) const;

	annotate::CycleInfo::residue_strand getResIds(
		const std::string& aResidues) const;

	annotate::InteractionInfo::enFace getFace(
		const std::string& astrGeneralFaces) const;
	annotate::BasePair::enOrientation getOrientation(
		const std::string& astrOrientation) const;
};


#endif /* _t3tablebuilder_H_ */
