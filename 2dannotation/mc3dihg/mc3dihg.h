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

class InteractionHypothesis
{
public:
	InteractionHypothesis(
		const std::string& astrCycle1,
		const std::string& astrCycle2,
		unsigned int auiType,
		unsigned int auiPos1,
		unsigned int auiPos2) :
	mstrCycle1(astrCycle1),
	mstrCycle2(astrCycle2),
	muiType(auiType),
	muiPos1(auiPos1),
	muiPos2(auiPos2)
	{
	}

	bool operator <(const InteractionHypothesis& aRight) const
	{
		bool bLess = false;
		if(mstrCycle1 < aRight.mstrCycle1)
		{
			bLess = true;
		}else if(mstrCycle1 == aRight.mstrCycle1)
		{
			if(mstrCycle2 < aRight.mstrCycle2)
			{
				bLess = true;
			}else if(mstrCycle2 == aRight.mstrCycle2)
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

	std::string mstrCycle1;
	std::string mstrCycle2;
	unsigned int muiType;
	unsigned int muiPos1;
	unsigned int muiPos2;
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
	typedef std::pair<std::string, std::string> cycle_key;
	typedef std::pair<cycle_key, cycle_key> junction_key;
	typedef std::map<std::string, cycle_key> indexed_cycles;

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
	std::vector<std::vector<std::vector<float> > > mTerm2; // NCM x NCM x type
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

	cycle_key cycleKey(const std::string& astrKey) const;
	std::string cycleKeyToString(const cycle_key& aKey) const;
	cycle_key forceSymmetric(const cycle_key& aKey) const;
	std::string flipSequence(unsigned int auiStrand1, unsigned int auiStrand2, const std::string& astrSequence) const;

	unsigned int profileSize(unsigned int auiProfileId) const;

	unsigned int interactionTypeStringToId(const std::string& astrType) const;
	std::string interactionTypeIdToString(unsigned int auiTypeId) const;

	std::string nucleotideString(unsigned int auiNucId) const;
	unsigned int nucleotideId(const std::string& astrNucleotide) const;

	std::string getNucleotide(const cycle_key& aCycle, unsigned int auiPosition) const;

	void computeScores();
	void displayHypothesis() const;
};
#endif // _mc3dihg_H_
