/*
 * ScoringTermPositionKnowingNCMs.cc
 *
 *  Created on: Jun 2, 2010
 *      Author: blanchmf
 */

#include "ScoringTermPositionKnowingNCMs.h"

#include "StringUtil.h"

#include <cstdlib>
#include <fstream>
#include <iostream>

ScoringTermPositionKnowingNCMs::ScoringTermPositionKnowingNCMs()
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
}

void ScoringTermPositionKnowingNCMs::read(const std::string& astrFile)
{
	unsigned int uiNbProfiles = mProfileMap.size();
	mFrequencies.resize(uiNbProfiles);
	for(unsigned int i = 0; i < uiNbProfiles; ++ i)
	{
		// By default assume everything is probable
		mFrequencies[i].resize(uiNbProfiles);
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
					unsigned int uiProfile1 = mProfileMap[strProfile1];
					unsigned int uiProfile2 = mProfileMap[strProfile2];
					unsigned int uiPos1 = atol(fields[2].c_str());
					unsigned int uiPos2 = atol(fields[3].c_str());
					float fFreq = atof(fields[4].c_str());

					interaction_coordinate coords(uiPos1, uiPos2);
					mFrequencies[uiProfile1][uiProfile2][coords] = fFreq;
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

float ScoringTermPositionKnowingNCMs::score(
	const enProfile& aeProfile1,
	const enProfile& aeProfile2,
	unsigned int auiPos1,
	unsigned int auiPos2,
	unsigned int auiTypeId,
	const enNucleotide& aeNuc1,
	const enNucleotide& aeNuc2) const
{
	float fScore = 0.0f;
	interaction_coordinate coords(auiPos1, auiPos2);
	std::map<interaction_coordinate, float>::const_iterator it;
	it = mFrequencies[aeProfile1][aeProfile2].find(coords);
	if(it != mFrequencies[aeProfile1][aeProfile2].end())
	{
		fScore = it->second;
	}
	return fScore;
}
