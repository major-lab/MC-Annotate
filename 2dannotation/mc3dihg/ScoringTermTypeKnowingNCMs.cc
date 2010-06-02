/*
 * ScoringTermTypeKnowingNCMs.cc
 *
 *  Created on: Jun 2, 2010
 *      Author: blanchmf
 */

#include "ScoringTermTypeKnowingNCMs.h"

#include "StringUtil.h"

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>

ScoringTermTypeKnowingNCMs::ScoringTermTypeKnowingNCMs()
{
	// Initialize the profile map
	mProfileMap["3"] = eProfile3;
	mProfileMap["4"] = eProfile4;
	mProfileMap["5"] = eProfile5;
	mProfileMap["6"] = eProfile6;
	mProfileMap["2_2"] = eProfile2_2;
	mProfileMap["2_3"] = eProfile2_3;
	mProfileMap["2_4"] = eProfile2_4;
	mProfileMap["2_5"] = eProfile2_5;
	mProfileMap["2_6"] = eProfile2_6;
	mProfileMap["3_3"] = eProfile3_3;
	mProfileMap["3_4"] = eProfile3_4;
	mProfileMap["3_5"] = eProfile3_5;
	mProfileMap["4_4"] = eProfile4_4;

	// Initialize the interaction type map
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eWatson, annotate::InteractionInfo::eWatson), annotate::BasePair::eCis), 		0));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eWatson, annotate::InteractionInfo::eWatson), annotate::BasePair::eTrans), 	1));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eWatson, annotate::InteractionInfo::eHoogsteen), annotate::BasePair::eCis), 	2));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eWatson, annotate::InteractionInfo::eHoogsteen), annotate::BasePair::eTrans), 	3));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eWatson, annotate::InteractionInfo::eSugar), annotate::BasePair::eCis), 		4));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eWatson, annotate::InteractionInfo::eSugar), annotate::BasePair::eTrans), 		5));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eWatson, annotate::InteractionInfo::eRibose), annotate::BasePair::eUnknown), 	6));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eWatson, annotate::InteractionInfo::ePhosphate), annotate::BasePair::eUnknown), 7));

	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eHoogsteen, annotate::InteractionInfo::eWatson), annotate::BasePair::eCis), 		8));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eHoogsteen, annotate::InteractionInfo::eWatson), annotate::BasePair::eTrans), 		9));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eHoogsteen, annotate::InteractionInfo::eHoogsteen), annotate::BasePair::eCis), 	10));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eHoogsteen, annotate::InteractionInfo::eHoogsteen), annotate::BasePair::eTrans), 	11));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eHoogsteen, annotate::InteractionInfo::eSugar), annotate::BasePair::eCis), 		12));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eHoogsteen, annotate::InteractionInfo::eSugar), annotate::BasePair::eTrans), 		13));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eHoogsteen, annotate::InteractionInfo::eRibose), annotate::BasePair::eUnknown), 	14));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eHoogsteen, annotate::InteractionInfo::ePhosphate), annotate::BasePair::eUnknown), 	15));

	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eSugar, annotate::InteractionInfo::eWatson), annotate::BasePair::eCis), 		16));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eSugar, annotate::InteractionInfo::eWatson), annotate::BasePair::eTrans), 		17));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eSugar, annotate::InteractionInfo::eHoogsteen), annotate::BasePair::eCis), 	18));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eSugar, annotate::InteractionInfo::eHoogsteen), annotate::BasePair::eTrans), 	19));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eSugar, annotate::InteractionInfo::eSugar), annotate::BasePair::eCis), 		20));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eSugar, annotate::InteractionInfo::eSugar), annotate::BasePair::eTrans), 		21));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eSugar, annotate::InteractionInfo::eRibose), annotate::BasePair::eUnknown), 	22));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eSugar, annotate::InteractionInfo::ePhosphate), annotate::BasePair::eUnknown), 	23));

	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eRibose, annotate::InteractionInfo::eWatson), annotate::BasePair::eUnknown), 		24));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eRibose, annotate::InteractionInfo::eHoogsteen), annotate::BasePair::eUnknown), 	25));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eRibose, annotate::InteractionInfo::eSugar), annotate::BasePair::eUnknown), 		26));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eRibose, annotate::InteractionInfo::eRibose), annotate::BasePair::eUnknown), 		27));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::eRibose, annotate::InteractionInfo::ePhosphate), annotate::BasePair::eUnknown), 	28));

	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::ePhosphate, annotate::InteractionInfo::eWatson), annotate::BasePair::eUnknown), 	29));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::ePhosphate, annotate::InteractionInfo::eHoogsteen), annotate::BasePair::eUnknown), 	30));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::ePhosphate, annotate::InteractionInfo::eSugar), annotate::BasePair::eUnknown), 		31));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::ePhosphate, annotate::InteractionInfo::eRibose), annotate::BasePair::eUnknown), 	32));
	mInteractionTypeMap.insert(std::pair<interaction_type, unsigned int>(interaction_type(face_pair(annotate::InteractionInfo::ePhosphate, annotate::InteractionInfo::ePhosphate), annotate::BasePair::eUnknown), 	33));
}

void ScoringTermTypeKnowingNCMs::read(const std::string& astrFile)
{
	unsigned int uiNbProfiles = mProfileMap.size();
	mFrequencies.resize(uiNbProfiles);
	for(unsigned int i = 0; i < uiNbProfiles; ++ i)
	{
		// By default assume everything is probable
		mFrequencies[i].resize(uiNbProfiles);
		for(unsigned int j = 0; j < uiNbProfiles; ++ j)
		{
			mFrequencies[i][j].resize(mInteractionTypeMap.size(), 1.0f);
		}
	}

	if(0 < astrFile.size())
	{
		std::ifstream infile;
		infile.open(astrFile.c_str(), std::ios_base::in);
		if(infile.good())
		{
			std::string strLine;
			while(std::getline(infile, strLine).good())
			{
				// Remove comments
				strLine = annotate::cutStringAfter(strLine, "//");
				if(0 < strLine.size())
				{
					// Remove whitespaces
					annotate::cleanString(strLine, ' ');

					// Read the data
					std::vector<std::string> fields = annotate::splitStringFields(strLine, ":");
					std::string strProfile1 = fields[0];
					std::string strProfile2 = fields[1];
					std::string strInterType = fields[2];
					float fFreq = atof(fields[3].c_str());

					unsigned int uiTypeId = interactionTypeStringToId(strInterType);
					enProfile eProfile1 = mProfileMap[strProfile1];
					enProfile eProfile2 = mProfileMap[strProfile2];
					mFrequencies[eProfile1][eProfile2][uiTypeId] = fFreq;
				}
			}
		}
		else
		{
			// TODO : This should be an exception
			std::cout << "Error opening file " << astrFile << std::endl;
		}
		infile.close();
	}
}

float ScoringTermTypeKnowingNCMs::score(
	const enProfile& aeProfile1,
	const enProfile& aeProfile2,
	unsigned int auiPos1,
	unsigned int auiPos2,
	unsigned int auiTypeId,
	const enNucleotide& aeNuc1,
	const enNucleotide& aeNuc2) const
{
	float fScore = mFrequencies[aeProfile1][aeProfile2][auiTypeId];
	return fScore;
}

unsigned int ScoringTermTypeKnowingNCMs::interactionTypeStringToId(const std::string& astrType) const
{
	const char* szTypes[] = {
		"cWW", "tWW",   "cWH", "tWH", "cWS", "tWS", "uWR", "uWP",
		"cHW", "tHW",   "cHH", "tHH", "cHS", "tHS", "uHR", "uHP",
		"cSW", "tSW",   "cSH", "tSH", "cSS", "tSS", "uSR", "uSP",
		"uRW", "uRH",   "uRS", "uRR", "uRP",
		"uPW", "uPH",   "uPS", "uPR", "uPP"
		};
	unsigned int i = 0;
	for(i = 0; i < mInteractionTypeMap.size(); ++ i)
	{
		if(astrType == szTypes[i])
		{
			break;
		}
	}
	assert(i < mInteractionTypeMap.size());
	return i;
}

