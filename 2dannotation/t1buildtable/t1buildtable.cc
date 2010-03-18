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

static const char* shortopts = "Vhlvc:s:";

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
		<< " [-hlvV] -s <cycle count> -c <3D interacting cycle pairs>"
		<< std::endl;
}

void T1BuildTable::help () const
{
	mccore::gOut (0)
		<< "This program computes a table of frequency of interacting cycle pairs." << std::endl
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

	if(0 == mstrCyclesCountFile.size() || 0 == mstrCyclesPairsFile.size())
	{
		usage ();
		exit (EXIT_FAILURE);
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
			unsigned int uiProfile = profileId(fields[0]);
			mCycleCount[uiProfile] = atol(fields[1].c_str());
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

		if(2 < fields.size())
		{
			std::vector<std::string> fields1 =  annotate::splitStringFields(fields[1], ":");
			std::vector<std::string> fields2 =  annotate::splitStringFields(fields[2], ":");

			unsigned int uiProfile1 = profileId(fields1[0]);
			unsigned int uiProfile2 = profileId(fields2[0]);

			mCyclesPairsCount[uiProfile1][uiProfile2] += 1;
			mCyclesPairsCount[uiProfile2][uiProfile1] += 1;
		}
	}
	infile.close();

	computeNormalizedFrequencies();
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
		unsigned int uiCount1 = mCycleCount[i];
		for(unsigned int j = 0; j < uiNbProfiles; ++ j)
		{
			unsigned int uiCount2 = mCycleCount[j];
			unsigned int uiInterCount = mCyclesPairsCount[i][j];
			unsigned int uiTotalCount = uiCount1 * uiCount2;
			float fFreq = ((float)(uiInterCount)) / ((float)(uiTotalCount));
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
