/*
 * t3tablebuilder.cc
 *
 *  Created on: May 1, 2010
 *      Author: blanchmf
 */

#include "t3tablebuilder.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include "StringUtil.h"
#include "StringTable.h"

T3TableBuilder::T3TableBuilder(
	const std::string& astrGroupFile,
	const std::string& astrInteractions)
{
	// Initialize the profile map
	mProfileMap["3"] = 0;
	mProfileMap["4"] = 1;
	mProfileMap["5"] = 2;
	mProfileMap["6"] = 3;
	mProfileMap["2_2"] = 4;
	mProfileMap["2_3"] = 5;
	mProfileMap["2_4"] = 6;
	mProfileMap["2_5"] = 7;
	mProfileMap["2_6"] = 8;
	mProfileMap["3_3"] = 9;
	mProfileMap["3_4"] = 10;
	mProfileMap["3_5"] = 11;
	mProfileMap["4_4"] = 12;

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

	// Read the PDB groups
	readPDBGroups(astrGroupFile);

	// Read the interactions
	readInteractions(astrInteractions);
}

void T3TableBuilder::computeTable()
{
	// Compute the statistics
	computeStatistics();

	// Display the frequencies
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		for(unsigned int j = 0; j < mProfileMap.size(); ++ j)
		{
			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				for(unsigned int l = 0; l < mFrequencies[i][j][k].size(); ++ l)
				{
					for(unsigned int m = 0; m < mFrequencies[i][j][k][l].size(); ++ m)
					{
						std::cout << profileString(i);
						std::cout << " : " << profileString(j);
						std::cout << " : " << interactionTypeString(k);
						std::cout << " : " << l;
						std::cout << " : " << m;
						std::cout << " : " << mFrequencies[i][j][k][l][m] << std::endl;
					}
				}
			}
		}
	}
}

void T3TableBuilder::readPDBGroups(const std::string& astrGroupFile)
{
	// Read the file
	mGroupFile.read(astrGroupFile);
	std::map<unsigned int, std::list<annotate::RNAGroupFileEntry> >::const_iterator it;
	for(it = mGroupFile.groups().begin(); it != mGroupFile.groups().end(); ++ it)
	{
		mGroupModelMap[it->first] = annotate::RNAGroupModel(it->second);
	}
}

unsigned int T3TableBuilder::profileId(const std::string& astrProfile) const
{
	unsigned int uiProfId = 0;
	std::map<std::string, unsigned int>::const_iterator it;
	it = mProfileMap.find(astrProfile);
	if(it != mProfileMap.end())
	{
		uiProfId = it->second;
	}else
	{
		std::cout << "Could not find profile : " << astrProfile << std::endl;
		assert(false);
	}
	return uiProfId;
}

void T3TableBuilder::readInteractions(const std::string& astrInteractionsFile)
{
	std::ifstream infile;
	infile.open (astrInteractionsFile.c_str(), std::ifstream::in);
	std::string strLine;

	mInteractions.clear();

	while (std::getline(infile, strLine).good())
	{
		annotate::cleanString(strLine, ' ');
		std::vector<std::string> fields = annotate::splitStringFields(strLine, ";");
		if(3 == fields.size())
		{
			std::vector<std::string> interactionFields;
			interactionFields = annotate::splitStringFields(fields[0], ":");
			std::string strPDB = interactionFields[0];
			unsigned int uiModel = atol(interactionFields[1].c_str());
			mccore::ResId resId1 = mccore::ResId(interactionFields[2].c_str());
			mccore::ResId resId2 = mccore::ResId(interactionFields[3].c_str());
			annotate::InteractionInfo::enFace eFace1 = getFace(interactionFields[4]);
			annotate::InteractionInfo::enFace eFace2 = getFace(interactionFields[5]);
			annotate::BasePair::enOrientation eOrientation = getOrientation(interactionFields[6]);

			// Read the interaction
			annotate::InteractionInfo inter(
				strPDB,
				uiModel,
				resId1,
				resId2,
				eFace1,
				eFace2,
				eOrientation);

			// Read the cycles associated
			std::vector<std::string> cycle1Fields;
			std::vector<std::string> cycle2Fields;
			cycle1Fields =  annotate::splitStringFields(fields[1], ":");
			cycle2Fields =  annotate::splitStringFields(fields[2], ":");

			annotate::CycleInfo cycle1 = getCycleInfo(strPDB, uiModel, cycle1Fields[0], cycle1Fields[1], cycle1Fields[2]);
			annotate::CycleInfo cycle2 = getCycleInfo(strPDB, uiModel, cycle2Fields[0], cycle2Fields[1], cycle2Fields[2]);

			if(cycle1.contains(inter.resId1()) && cycle2.contains(inter.resId2()))
			{
				mInteractions.insert(interaction_cycle_pair(inter, cycle1, cycle2));
			}else if(cycle2.contains(inter.resId1()) && cycle1.contains(inter.resId2()))
			{
				mInteractions.insert(interaction_cycle_pair(inter, cycle2, cycle1));
			}else
			{
				std::cout << "The interaction is invalid" << std::endl;
				assert(false);
			}
		}
	}
	infile.close();
}

annotate::BasePair::enOrientation T3TableBuilder::getOrientation(
	const std::string& astrOrientation) const
{
	annotate::BasePair::enOrientation eOrientation = annotate::BasePair::eUnknown;
	if(astrOrientation == "Cis")
	{
		eOrientation = annotate::BasePair::eCis;
	}else if(astrOrientation == "Trans")
	{
		eOrientation = annotate::BasePair::eTrans;
	}else if(astrOrientation == "Unknown")
	{
		eOrientation = annotate::BasePair::eUnknown;
	}else
	{
		assert(false); // There is a problem with the file
	}
	return eOrientation;
}

annotate::InteractionInfo::enFace T3TableBuilder::getFace(
	const std::string& astrGeneralFaces) const
{
	annotate::InteractionInfo::enFace eFace = annotate::InteractionInfo::eWatson;
	if(astrGeneralFaces == "Watson")
	{
		eFace = annotate::InteractionInfo::eWatson;
	}else if(astrGeneralFaces == "Hoogsteen")
	{
		eFace = annotate::InteractionInfo::eHoogsteen;
	}else if(astrGeneralFaces == "Sugar")
	{
		eFace = annotate::InteractionInfo::eSugar;
	}else if(astrGeneralFaces == "Ribose")
	{
		eFace = annotate::InteractionInfo::eRibose;
	}
	else if(astrGeneralFaces == "Phosphate")
	{
		eFace = annotate::InteractionInfo::ePhosphate;
	}else
	{
		std::cout << "Failed to convert face : " << astrGeneralFaces << std::endl;
		assert(false); // There is a problem with the file
	}
	return eFace;
}


annotate::CycleInfo T3TableBuilder::getCycleInfo(
	const std::string& astrPDB,
	unsigned int uiModel,
	const std::string& astrProfile,
	const std::string& astrResIds,
	const std::string& astrSequence) const
{
	annotate::CycleProfile profile(astrProfile);
	std::vector<std::vector<mccore::ResId> > resIds;
	resIds.resize(profile.strandProfile().size());
	std::list<unsigned int>::const_iterator it = profile.strandProfile().begin();
	unsigned int iCursor = 0;
	unsigned int iStrand = 0;
	// std::vector<mccore::ResId> resIdsFields = annotate::splitStringFields(astrResIds, "-");
	std::vector<mccore::ResId> resIdsFields = getResIds(astrResIds);
	for(; it != profile.strandProfile().end(); ++ it, ++iStrand)
	{
		resIds[iStrand].resize(*it);
		for(unsigned int i = 0; i < *it; ++ i)
		{
			resIds[iStrand][i] = resIdsFields[i + iCursor];
		}
		iCursor += *it;
	}
	std::vector<std::string> sequenceFields = annotate::splitStringFields(astrResIds, "-");
	return annotate::CycleInfo(astrPDB, uiModel, profile, profile, resIds, sequenceFields);
}

annotate::CycleInfo::residue_strand T3TableBuilder::getResIds(
	const std::string& aResidues) const
{
	annotate::CycleInfo::residue_strand residues;
	std::string strResidues = aResidues;

	annotate::cleanString(strResidues, ' ');
	std::vector<std::string> residuesString = annotate::splitStringFields(strResidues, "-");
	std::vector<std::string>::const_iterator it = residuesString.begin();
	for(it = residuesString.begin(); it != residuesString.end(); ++it)
	{
		mccore::ResId res(it->c_str());
		residues.push_back(res);
	}
	return residues;
}

unsigned int T3TableBuilder::profileId(const annotate::CycleInfo& aCycle) const
{
	std::string strProfile = aCycle.getProfile().toString();
	return profileId(strProfile);
}

unsigned int T3TableBuilder::interactionTypeId(const annotate::InteractionInfo& aInteraction) const
{
	face_pair faces(aInteraction.face1(), aInteraction.face2());
	interaction_type interType(faces, aInteraction.orientation());

	unsigned int uiTypeId = 0;
	std::map<interaction_type, unsigned int>::const_iterator it;
	it = mInteractionTypeMap.find(interType);
	if(it != mInteractionTypeMap.end())
	{
		uiTypeId = it->second;
	}else
	{
		std::cout << "Could not find interaction type : (" << aInteraction.face1();
				std::cout << "," << aInteraction.face2() << "," << aInteraction.orientation();
				std::cout << ")" << std::endl;
		assert(false);
	}
	return uiTypeId;
}

unsigned int T3TableBuilder::flipInteractionTypeId(const annotate::InteractionInfo& aInteraction) const
{
	face_pair faces(aInteraction.face2(), aInteraction.face1());
	interaction_type interType(faces, aInteraction.orientation());

	unsigned int uiTypeId = 0;
	std::map<interaction_type, unsigned int>::const_iterator it;
	it = mInteractionTypeMap.find(interType);
	if(it != mInteractionTypeMap.end())
	{
		uiTypeId = it->second;
	}else
	{
		std::cout << "Could not find interaction type : (" << aInteraction.face2();
		std::cout << "," << aInteraction.face1() << "," << aInteraction.orientation();
		std::cout << ")" << std::endl;
		assert(false);
	}
	return uiTypeId;
}

std::string T3TableBuilder::profileString(unsigned int auiProfile) const
{
	std::string strProfile;

	switch(auiProfile)
	{
	case 0:
		strProfile = "3";
		break;
	case 1:
		strProfile = "4";
		break;
	case 2:
		strProfile = "5";
		break;
	case 3:
		strProfile = "6";
		break;

	case 4:
		strProfile = "2_2";
		break;
	case 5:
		strProfile = "2_3";
		break;
	case 6:
		strProfile = "2_4";
		break;
	case 7:
		strProfile = "2_5";
		break;
	case 8:
		strProfile = "2_6";
		break;
	case 9:
		strProfile = "3_3";
		break;
	case 10:
		strProfile = "3_4";
		break;
	case 11:
		strProfile = "3_5";
		break;
	case 12:
		strProfile = "4_4";
		break;
	}
	return strProfile;
}

std::string T3TableBuilder::interactionTypeString(unsigned int auiTypeId) const
{
	const char* szTypes[] = {
		"cWW", "tWW",   "cWH", "tWH", "cWS", "tWS", "uWR", "uWP",
		"cHW", "tHW",   "cHH", "tHH", "cHS", "tHS", "uHR", "uHP",
		"cSW", "tSW",   "cSH", "tSH", "cSS", "tSS", "uSR", "uSP",
		"uRW", "uRH",   "uRS", "uRR", "uRP",
		"uPW", "uPH",   "uPS", "uPR", "uPP"
		};
	std::string strType = szTypes[auiTypeId];
	return strType;
}

T3TableBuilder::interaction_coordinate T3TableBuilder::getInteractionCoordinates(
	const interaction_cycle_pair& aInteraction) const
{
	interaction_coordinate coords;
	coords.first = getNucleotidePosition(aInteraction.second, aInteraction.first.resId1());
	coords.second = getNucleotidePosition(aInteraction.third, aInteraction.first.resId2());
	return coords;
}

unsigned int T3TableBuilder::getNucleotidePosition(
	const annotate::CycleInfo& aCycle,
	const mccore::ResId& aResId) const
{
	unsigned int uiPosition = 0;
	unsigned int i = 0;
	std::vector<mccore::ResId> resIds = aCycle.getResIds();
	for(i = 0; i < resIds.size(); ++ i)
	{
		if(aResId == resIds[i])
		{
			break;
		}
	}
	assert(i < resIds.size());
	unsigned int uiNbStrands = aCycle.getNbStrands();
	assert(0 < uiNbStrands && uiNbStrands < 3);
	uiPosition = i;
	if(2 == uiNbStrands)
	{
			annotate::CycleInfo::residue_profile resProf = aCycle.getStrandResIds();
		if(resProf[0].size() == resProf[1].size())
		{
			// Symmetric cycle
			uiPosition %= resProf[0].size();
		}
	}
	return uiPosition;
}

unsigned int T3TableBuilder::profileSize(unsigned int auiProfileId) const
{
	unsigned int uiSize = 0;

	switch(auiProfileId)
	{
	case 0:
		uiSize = 3;
		break;
	case 1:
		uiSize = 4;
		break;
	case 2:
		uiSize = 5;
		break;
	case 3:
		uiSize = 6;
		break;
	case 4:
		uiSize = 4;
		break;
	case 5:
		uiSize = 5;
		break;
	case 6:
		uiSize = 6;
		break;
	case 7:
		uiSize = 7;
		break;
	case 8:
		uiSize = 8;
		break;
	case 9:
		uiSize = 6;
		break;
	case 10:
		uiSize = 7;
		break;
	case 11:
		uiSize = 8;
		break;
	case 12:
		uiSize = 8;
		break;
	}
	return uiSize;
}
