/*
 * mc3dihg.cc
 *
 *  Created on: May 3, 2010
 *      Author: blanchmf
 */

#include "mcsincm.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>

#include "StringTable.h"
#include "StringUtil.h"

static const char* gszShortopts = "Vhi:a:b:";

MCScoreInteractingNCM::MCScoreInteractingNCM(int argc, char * argv [])
: mfAppVersion(1.0f),
mstrAppName("MC3MCScoreInteractingNCM")
{
	read_options (argc, argv);
	mInteractions = readInteractionsFile(mstrInteractionsFile);
	mCycles1 = readIndexedCycleFile(mstrCyclesFile1);
	mCycles2 = readIndexedCycleFile(mstrCyclesFile2);

	initializeTables();

	computeScores();

	std::map<std::pair<std::string, std::string>, float>::const_iterator it2 = mScores.begin();
	for(; it2 != mScores.end(); ++ it2)
	{
		std::cout << it2->second;
		std::cout << " : " << it2->first.first;
		std::cout << " : " << it2->first.second;
		std::cout << std::endl;
	}
}


void MCScoreInteractingNCM::version () const
{
	std::cout << mstrAppName << " "; // Nom du logiciel
	std::cout << mfAppVersion << " ";			// Version du logiciel
	std::cout << "(" << __DATE__ << ")";	// Date de la compilation
	std::cout << std::endl;
}

void MCScoreInteractingNCM::usage () const
{
	std::cout << "usage: "
		<< " [-hV] -i <scored interaction file> -a <cycles from first domain> -b <cycles from second domain>"
		<< std::endl;
}


void MCScoreInteractingNCM::help () const
{
	std::cout
		<< "This program scores the NCMs pairs according to the contributing interactions." << std::endl
		<< "  -i	       File containing the scored hypothetical interactions." << std::endl
		<< "  -a <file> File containing the cycles from the first domain." << std::endl
		<< "  -b <file> File containing the cycles from the second domain." << std::endl;
}

void MCScoreInteractingNCM::read_options (int argc, char* argv[])
{
	int c;

	while ((c = getopt (argc, argv, gszShortopts)) != -1)
	{
		switch (c)
		{
		case 'V':
			version ();
			exit (EXIT_SUCCESS);
			break;
		case 'h':
			help();
			exit(EXIT_SUCCESS);
			break;
		case 'i':
			mstrInteractionsFile = optarg;
			break;
		case 'a':
			mstrCyclesFile1 = optarg;
			break;
		case 'b':
			mstrCyclesFile2 = optarg;
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if (0 == mstrInteractionsFile.size())
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

std::list<MCScoreInteractingNCM::interaction_entry> MCScoreInteractingNCM::readInteractionsFile(const std::string astrFileName) const
{
	std::list<MCScoreInteractingNCM::interaction_entry> interactions;

	std::ifstream infile;
	infile.open(astrFileName.c_str(), std::ios_base::in);
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
				assert(4 == fields.size());
				float fScore = atof(fields[0].c_str());
				std::string strFaces = fields[1];
				unsigned int uiRes1 = atol(fields[2].c_str());
				unsigned int uiRes2 = atol(fields[3].c_str());

				resid_pair resIds(uiRes1, uiRes2);
				interactions.push_back(interaction_entry(std::pair<float, std::string>(fScore, strFaces), resIds));
			}
		}
	}
	else
	{
		// TODO : This should be an exception
		std::cout << "Error opening file " << astrFileName << std::endl;
	}
	infile.close();

	return interactions;
}

MCScoreInteractingNCM::indexed_cycles MCScoreInteractingNCM::readIndexedCycleFile(
	const std::string& astrFileName) const
{
	std::list<annotate::CycleInfo> indexedCycleFile;
	std::ifstream infile;
	infile.open(astrFileName.c_str(), std::ios_base::in);
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
				assert(5 == fields.size());
				std::string strIndex = fields[0].c_str();
				std::string strProfile = fields[1];
				std::string strSequence = fields[2];
				std::string strIds = fields[3];
				annotate::CycleInfo cycle = getCycleInfo(strIndex, strProfile, strIds, strSequence);
				indexedCycleFile.push_back(cycle);
			}
		}
	}
	else
	{
		// TODO : This should be an exception
		std::cout << "Error opening file " << astrFileName << std::endl;
	}
	infile.close();
	indexed_cycles cycles;
	cycles.resize(indexedCycleFile.size());
	std::copy(indexedCycleFile.begin(), indexedCycleFile.end(), cycles.begin());
	return cycles;
}

annotate::CycleInfo MCScoreInteractingNCM::getCycleInfo(
	const std::string& astrIdentifier,
	const std::string& astrProfile,
	const std::string& astrResIds,
	const std::string& astrSequence) const
{
	annotate::CycleProfile fileProfile(astrProfile);
	std::vector<std::vector<mccore::ResId> > resIds;
	resIds.resize(fileProfile.strandProfile().size());
	std::vector<unsigned int>::const_iterator it = fileProfile.strandProfile().begin();
	unsigned int iCursor = 0;
	unsigned int iStrand = 0;
	std::vector<mccore::ResId> resIdsFields = getResIds(astrResIds);
	for(; it != fileProfile.strandProfile().end(); ++ it, ++iStrand)
	{
		resIds[iStrand].resize(*it);
		for(unsigned int i = 0; i < *it; ++ i)
		{
			resIds[iStrand][i] = resIdsFields[i + iCursor];
		}
		iCursor += *it;
	}
	std::vector<std::string> sequenceFields;

	annotate::CycleProfile profile = fileProfile;
	if(resIds.size() == 2 && resIds[1].size() < resIds[0].size())
	{
		// We need to flip the profile, resids and sequence
		std::vector<mccore::ResId> strand1 = resIds[1];
		std::vector<mccore::ResId> strand2 = resIds[0];
		resIds[0] = strand1;
		resIds[1] = strand2;

		sequenceFields = annotate::splitStringFields(flipSequence(resIds[1].size(), resIds[0].size(), astrSequence), "-");
		profile = annotate::CycleProfile::Rotate(profile);
	}else
	{
		sequenceFields = annotate::splitStringFields(astrSequence, "-");
	}
	return annotate::CycleInfo(astrIdentifier, 0, fileProfile, profile, resIds, sequenceFields);
}

// TODO : This is duplicate from table builders
annotate::CycleInfo::residue_strand MCScoreInteractingNCM::getResIds(
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

std::string MCScoreInteractingNCM::flipSequence(
	unsigned int auiStrand1,
	unsigned int auiStrand2,
	const std::string& astrSequence) const
{
	std::ostringstream oss;
	std::vector<std::string> fields = annotate::splitStringFields(astrSequence, "-");
	unsigned int i;
	// Copy strand 2 at beginning
	for(i = auiStrand1; i < fields.size(); i ++)
	{
		if(i != auiStrand1)
		{
			oss << "-";
		}
		oss << fields[i];
	}
	// Copy strand 1 after
	for(i = 0; i < auiStrand1; i ++)
	{
		oss << "-" << fields[i];
	}
	return oss.str();
}

void MCScoreInteractingNCM::initializeTables()
{
	// Initialize the score table
	indexed_cycles::const_iterator it1;
	for(it1 = mCycles1.begin(); it1 != mCycles1.end(); ++ it1)
	{
		std::string strCycle = it1->toString();
		std::vector<mccore::ResId>::const_iterator itRes = it1->getResIds().begin();
		for(;itRes != it1->getResIds().end(); ++ itRes)
		{
			mResCycleMap1[*itRes].insert(strCycle);
		}
	}

	indexed_cycles::const_iterator it2;
	for(it2 = mCycles2.begin(); it2 != mCycles2.end(); ++ it2)
	{
		std::string strCycle = it2->toString();
		std::vector<mccore::ResId>::const_iterator itRes = it2->getResIds().begin();
		for(;itRes != it2->getResIds().end(); ++ itRes)
		{
			mResCycleMap2[*itRes].insert(strCycle);
		}
	}

	for(it1 = mCycles1.begin(); it1 != mCycles1.end(); ++ it1)
	{
		std::string strCycle1 = it1->toString();
		for(it2 = mCycles2.begin(); it2 != mCycles2.end(); ++ it2)
		{
			std::string strCycle2 = it2->toString();
			std::pair<std::string, std::string> entry(strCycle1, strCycle2);
			mScores[entry] = 0.0f;
		}
	}
}

void MCScoreInteractingNCM::computeScores()
{
	std::list<interaction_entry>::const_iterator it;
	for(it = mInteractions.begin(); it != mInteractions.end(); ++ it)
	{
		float fScore = it->first.first;
		mccore::ResId res1 = it->second.first;
		mccore::ResId res2 = it->second.second;

		std::set<std::string>::const_iterator itCycle1;
		std::set<std::string>::const_iterator itCycle1End = mResCycleMap1[res1].end();

		std::set<std::string>::const_iterator itCycle2;
		std::set<std::string>::const_iterator itCycle2End = mResCycleMap2[res2].end();

		for(itCycle1 = mResCycleMap1[res1].begin(); itCycle1 != itCycle1End; ++ itCycle1)
		{
			for(itCycle2 = mResCycleMap2[res2].begin(); itCycle2 != itCycle2End; ++ itCycle2)
			{
				std::pair<std::string, std::string> entry(*itCycle1, *itCycle2);
				mScores[entry] += fScore;
			}
		}
	}
}

int main(int argc, char* argv[])
{
	MCScoreInteractingNCM theApp(argc, argv);
	return EXIT_SUCCESS;
}
