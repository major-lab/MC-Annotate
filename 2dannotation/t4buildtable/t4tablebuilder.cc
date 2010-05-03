/*
 * t4tablebuilder.cc
 *
 *  Created on: May 2, 2010
 *      Author: blanchmf
 */

#include "t4tablebuilder.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include "StringUtil.h"
#include "StringTable.h"

T4TableBuilder::T4TableBuilder(
	const std::string& astrGroupFile,
	const std::string& astrInteractions)
{
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

void T4TableBuilder::computeTable()
{
	// Compute the statistics
	computeStatistics();

	// Display the frequencies
	for(unsigned int i = 0; i < 4; ++ i)
	{
		for(unsigned int j = 0; j < 4; ++ j)
		{
			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				std::cout << nucleotideString(i);
				std::cout << " : " << nucleotideString(j);
				std::cout << " : " << interactionTypeString(k);
				std::cout << " : " << mFrequencies[i][j][k] << std::endl;
			}
		}
	}
}

void T4TableBuilder::readPDBGroups(const std::string& astrGroupFile)
{
	// Read the file
	mGroupFile.read(astrGroupFile);
	std::map<unsigned int, std::list<annotate::RNAGroupFileEntry> >::const_iterator it;
	for(it = mGroupFile.groups().begin(); it != mGroupFile.groups().end(); ++ it)
	{
		mGroupModelMap[it->first] = annotate::RNAGroupModel(it->second);
	}
}

void T4TableBuilder::readInteractions(const std::string& astrInteractionsFile)
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

annotate::BasePair::enOrientation T4TableBuilder::getOrientation(
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

annotate::InteractionInfo::enFace T4TableBuilder::getFace(
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


annotate::CycleInfo T4TableBuilder::getCycleInfo(
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
	std::vector<std::string> sequenceFields = annotate::splitStringFields(astrSequence, "-");
	return annotate::CycleInfo(astrPDB, uiModel, profile, profile, resIds, sequenceFields);
}

annotate::CycleInfo::residue_strand T4TableBuilder::getResIds(
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

unsigned int T4TableBuilder::interactionTypeId(const annotate::InteractionInfo& aInteraction) const
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

unsigned int T4TableBuilder::flipInteractionTypeId(const annotate::InteractionInfo& aInteraction) const
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

std::string T4TableBuilder::interactionTypeString(unsigned int auiTypeId) const
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

std::string T4TableBuilder::nucleotideString(unsigned int auiNucId) const
{
	const char* szNucleotides[] = {"A", "C", "G", "U"};
	std::string strNucleotide = szNucleotides[auiNucId];
	return strNucleotide;
}
