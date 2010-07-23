//                              -*- Mode: C++ -*-
// t1buildtable.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Mar 17 11:14:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t1buildtable.h"

#include "StringUtil.h"
#include "StringTable.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

static const char* shortopts = "Vhlvc:s:g:";

T1BuildTable::T1BuildTable(int argc, char * argv [])
{
	// Read the command line options
	readOptions(argc, argv);

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

	// Read the PDB groups
	readPDBGroups();

	// Read the cycle count
	readCyclesCount();

	// Read the interacting pairs
	readCyclesPairs();

	// Output the frequencies
	outputFrequencies();
}

void T1BuildTable::version () const
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void T1BuildTable::usage () const
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-hlvV] -g <PDB Groups File> -s <cycle count file> -c <3D interacting cycle pairs file>"
		<< std::endl;
}

void T1BuildTable::help () const
{
	mccore::gOut (0)
		<< "This program computes a table of frequency of interacting cycle pairs." << std::endl
		<< "  -g    file containing the PDB File group for normalization" << std::endl
		<< "  -s	file containing the count of each type of NCM" << std::endl
		<< "  -c	file containing the identified 3D interacting cycles pairs" << std::endl
		<< "  -h	print this help" << std::endl
		<< "  -l	be more verbose (log)" << std::endl
		<< "  -v	be verbose" << std::endl
		<< "  -V	print the software version info" << std::endl;
}

void T1BuildTable::readOptions (int argc, char* argv[])
{
	int c;

	while ((c = getopt (argc, argv, shortopts)) != EOF)
	{
		switch (c)
		{
		case 'c':
		{
			mstrCyclesPairsFile = optarg;
			break;
		}
		case 'g':
		{
			mstrPDBGroupsFile = optarg;
			break;
		}
		case 's':
		{
			mstrCyclesCountFile = optarg;
			break;
		}
		case 'V':
			version ();
			exit (EXIT_SUCCESS);
			break;
		case 'h':
			usage ();
			help ();
			exit (EXIT_SUCCESS);
			break;
		case 'l':
			mccore::gErr.setVerboseLevel (mccore::gErr.getVerboseLevel () + 1);
			break;
		case 'v':
			mccore::gOut.setVerboseLevel (mccore::gOut.getVerboseLevel () + 1);
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if(0 == mstrCyclesCountFile.size() || 0 == mstrCyclesPairsFile.size() || 0 == mstrPDBGroupsFile.size())
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

void T1BuildTable::readPDBGroups()
{
	// Read the file
	mGroupFile.read(mstrPDBGroupsFile);
	std::map<unsigned int, std::list<annotate::RNAGroupFileEntry> >::const_iterator it;
	for(it = mGroupFile.groups().begin(); it != mGroupFile.groups().end(); ++ it)
	{
		mGroupModelMap[it->first] = annotate::RNAGroupModel(it->second);
	}
}

unsigned int T1BuildTable::profileId(const std::string& astrProfile) const
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

void T1BuildTable::readCyclesCount()
{
	std::ifstream infile;
	mCycleCount.clear();
	mCycleCount.resize(mProfileMap.size(), 0);

	infile.open (mstrCyclesCountFile.c_str(), std::ifstream::in);
	std::string strLine;
	while (std::getline(infile, strLine).good())
	{
		annotate::cleanString(strLine, ' ');
		std::vector<std::string> fields = annotate::splitStringFields(strLine, ":");
		if(1 < fields.size())
		{
			std::string strPDB = fields[0];
			unsigned int uiModel = atol(fields[1].c_str());
			annotate::CycleInfo cycle = getCycleInfo(strPDB, uiModel, fields[2], fields[3], fields[4]);
			unsigned int uiGroup = mGroupFile.getGroup(cycle);
			if(mGroupModelMap[uiGroup].isNew(cycle))
			{
				mGroupModelMap[uiGroup].addCycle(cycle);
				unsigned int uiProfile = profileId(fields[2]);
				mCycleCount[uiProfile] = 1 + mCycleCount[uiProfile];
			}
		}
	}
	infile.close();
}

void T1BuildTable::readCyclesPairs()
{
	std::ifstream infile;
	mCyclesPairsCount.clear();

	mCyclesPairsCount.resize(mProfileMap.size());
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		mCyclesPairsCount[i].resize(mProfileMap.size(), 0);
	}

	infile.open (mstrCyclesPairsFile.c_str(), std::ifstream::in);
	std::string strLine;

	while (std::getline(infile, strLine).good())
	{
		annotate::cleanString(strLine, ' ');
		std::vector<std::string> fields = annotate::splitStringFields(strLine, ";");

		std::vector<std::string> modelFields = annotate::splitStringFields(fields[0], ":");
		std::string strPDB = modelFields[0];
		unsigned int uiModel = atol(modelFields[1].c_str());

		if(2 < fields.size())
		{
			std::vector<std::string> fields1 =  annotate::splitStringFields(fields[1], ":");
			std::vector<std::string> fields2 =  annotate::splitStringFields(fields[2], ":");

			annotate::CycleInfo cycle1 = getCycleInfo(strPDB, uiModel, fields1[0], fields1[1], fields1[2]);
			annotate::CycleInfo cycle2 = getCycleInfo(strPDB, uiModel, fields2[0], fields2[1], fields2[2]);

			unsigned int uiProfile1 = profileId(fields1[0]);
			unsigned int uiProfile2 = profileId(fields2[0]);
			unsigned int uiGroup = mGroupFile.getGroup(cycle1);

			if(mGroupModelMap[uiGroup].isNew(cycle1, cycle2))
			{
				// Interaction has not been considered yet
				mGroupModelMap[uiGroup].addCyclePair(cycle1, cycle2);
				mCyclesPairsCount[uiProfile1][uiProfile2] += 1;
				mCyclesPairsCount[uiProfile2][uiProfile1] += 1;
			}
		}
	}
	infile.close();

	computeNormalizedFrequencies();
}

annotate::CycleInfo T1BuildTable::getCycleInfo(
	const std::string& astrPDB,
	unsigned int uiModel,
	const std::string& astrProfile,
	const std::string& astrResIds,
	const std::string& astrSequence) const
{
	annotate::CycleProfile profile(astrProfile);
	std::vector<std::vector<mccore::ResId> > resIds;
	resIds.resize(profile.strandProfile().size());
	std::vector<unsigned int>::const_iterator it = profile.strandProfile().begin();
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
	std::vector<std::string> sequenceFields = annotate::splitStringFields(astrResIds, "-");
	return annotate::CycleInfo(astrPDB, uiModel, profile, profile, resIds, sequenceFields);
}

annotate::CycleInfo::residue_strand T1BuildTable::getResIds(
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

void T1BuildTable::computeNormalizedFrequencies()
{
	unsigned int uiNbProfiles = mProfileMap.size();
	assert(0 < uiNbProfiles);
	mInteractionFrequencies.resize(uiNbProfiles);
	for(unsigned int k = 0; k < uiNbProfiles; ++ k)
	{
		mInteractionFrequencies[k].resize(uiNbProfiles, 0.0f);
	}
	for(unsigned int i = 0; i < uiNbProfiles; ++ i)
	{
		float fCount1 = (float)mCycleCount[i];
		for(unsigned int j = 0; j < uiNbProfiles; ++ j)
		{
			float fCount2 = (float)mCycleCount[j];
			float fInterCount = (float)mCyclesPairsCount[i][j];
			float fTotalCount = fCount1 * fCount2;
			float fFreq = fInterCount / fTotalCount;
			mInteractionFrequencies[i][j] = fFreq;
		}
	}
}

void T1BuildTable::outputFrequencies() const
{
	unsigned int uiNbProfiles = mProfileMap.size();

	std::vector<std::string> profiles;
	profiles.resize(uiNbProfiles);

	std::map<std::string, unsigned int>::const_iterator it;
	for(it = mProfileMap.begin(); it != mProfileMap.end(); ++ it)
	{
		profiles[it->second] = it->first;
	}

	annotate::StringTable stringTable(uiNbProfiles + 1);
	for(unsigned int i = 0; i < uiNbProfiles; ++ i)
	{
		std::vector<string>& tableRow = stringTable.addRow();
		tableRow[0] = profiles[i];
		for(unsigned int j = 0; j < uiNbProfiles; ++ j)
		{
			std::ostringstream oss;
			oss << mInteractionFrequencies[i][j];
			tableRow[j + 1] = oss.str();
		}
	}

	std::cout << stringTable.toString(" : ");
}

int main (int argc, char *argv[])
{
	T1BuildTable theApp(argc, argv);

	return EXIT_SUCCESS;
}
