/*
 * MC3DInteractionHypothesisGenerator.h
 *
 *  Created on: May 3, 2010
 *      Author: blanchmf
 */

#ifndef _mc3dihg_H_
#define _mc3dihg_H_

#include <string>
#include <map>
#include <vector>

#include "InteractionInfo.h"
#include "BasePair.h"
#include "CycleInfo.h"

#include "ScoringTermTypeKnowingNCMs.h"
#include "ScoringTermPositionKnowingNCMs.h"

class InteractionHypothesis
{
public:
	InteractionHypothesis(
		const annotate::CycleInfo& aCycle1,
		const annotate::CycleInfo& aCycle2,
		unsigned int auiType,
		unsigned int auiPos1,
		unsigned int auiPos2,
		float afProb1,
		float afProb2) :
	mpCycle1(&aCycle1),
	mpCycle2(&aCycle2),
	muiType(auiType),
	muiPos1(auiPos1),
	muiPos2(auiPos2),
	mProbabilities(afProb1, afProb2)
	{
	}

	bool operator <(const InteractionHypothesis& aRight) const
	{
		bool bLess = false;
		if(*mpCycle1 < *aRight.mpCycle1)
		{
			bLess = true;
		}else if(*mpCycle1 == *aRight.mpCycle1)
		{
			if(*mpCycle2 < *aRight.mpCycle2)
			{
				bLess = true;
			}else if(*mpCycle2 == *aRight.mpCycle2)
			{
				if(muiType < aRight.muiType)
				{
					bLess = true;
				}else if(muiType == aRight.muiType)
				{
					if(muiPos1 < aRight.muiPos1)
					{
						bLess = true;
					}else if(muiPos1 == aRight.muiPos1)
					{
						if(muiPos2 < aRight.muiPos2)
						{
							bLess = true;
						}else if(muiPos2 == aRight.muiPos2)
						{
							if(muiPos2 < aRight.muiPos2)
							{
								bLess = true;
							}
						}
					}
				}
			}
		}
		return bLess;
	}

	const annotate::CycleInfo* mpCycle1;
	const annotate::CycleInfo* mpCycle2;
	unsigned int muiType;
	unsigned int muiPos1;
	unsigned int muiPos2;
	std::pair<float, float> mProbabilities;
};

class MC3DInteractionHypothesisGenerator
{
public:
	MC3DInteractionHypothesisGenerator(int argc, char * argv []);

	// METHODS ---------------------------------------------
	void version () const;
	void usage () const;
	void help () const;
private:
	typedef std::pair<annotate::CycleInfo, float> cycle_prob_pair;
	typedef std::vector<cycle_prob_pair > indexed_cycles;

	typedef annotate::InteractionInfo::enFace face;
	typedef std::pair<face, face> face_pair;
	typedef std::pair<face_pair, annotate::BasePair::enOrientation> interaction_type;
	float mfAppVersion;
	std::string mstrAppName;

	std::map<std::string, std::map<std::string, float> > mScores;

	std::string mstrFirstCycles;
	std::string mstrSecondCycles;
	std::string mstrTerm1File;
	std::string mstrTerm2File;
	std::string mstrTerm3File;
	std::string mstrTerm4File;
	indexed_cycles mIndexedCycleFile1;
	indexed_cycles mIndexedCycleFile2;
	std::vector<std::vector<float> > mTerm1; // NCM x NCM
	ScoringTermPositionKnowingNCMs mTerm2;
	std::vector<std::vector<std::vector<std::vector<std::vector<float> > > > > mTerm3; // NCM x NCM x type x pos1 x pos2
	std::vector<std::vector<std::vector<float> > > mTerm4; // Nuc1 x Nuc2 x type

	std::map<std::string, unsigned int> mProfileMap;
	std::map<interaction_type, unsigned int> mInteractionTypeMap;
	std::map<std::string, unsigned int> mNucleotideStringIdMap;
	std::map<unsigned int, std::string> mNucleotideIdStringMap;
	std::multimap<float, InteractionHypothesis> mScoredHypothesis;

	void read_options (int argc, char* argv[]);

	void initializeMaps();
	void read_database();
	void readTerm1();
	void readTerm2();
	void readTerm3();
	void readTerm4();

	indexed_cycles readIndexedCycleFile(const std::string& astrFileName) const;

	std::string flipSequence(unsigned int auiStrand1, unsigned int auiStrand2, const std::string& astrSequence) const;

	unsigned int profileSize(unsigned int auiProfileId) const;

	unsigned int interactionTypeStringToId(const std::string& astrType) const;
	std::string interactionTypeIdToString(unsigned int auiTypeId) const;

	std::string nucleotideString(unsigned int auiNucId) const;
	unsigned int nucleotideId(const std::string& astrNucleotide) const;

	void computeScores();
	void displayHypothesis() const;
	annotate::CycleInfo::residue_strand getResIds(const std::string& aResidues) const;
	annotate::CycleInfo getCycleInfo(
		const std::string& astrIdentifier,
		const std::string& astrProfile,
		const std::string& astrResIds,
		const std::string& astrSequence) const;
};
#endif // _mc3dihg_H_
