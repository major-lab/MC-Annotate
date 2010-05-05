/*
 * mc3dihf.cc
 *
 *  Created on: May 4, 2010
 *      Author: blanchmf
 */

#include "mc3dihg.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>

#include "StringTable.h"
#include "StringUtil.h"

static const char* gszShortopts = "Vha:b:1:2:3:4:";

MC3DInteractionHypothesisGenerator::MC3DInteractionHypothesisGenerator(int argc, char * argv [])
: mfAppVersion(1.0f),
mstrAppName("MC3DInteractionHypothesisGenerator")
{
	read_options (argc, argv);

	initializeMaps();

	read_database();
	mIndexedCycleFile1 = readIndexedCycleFile(mstrFirstCycles);
	mIndexedCycleFile2 = readIndexedCycleFile(mstrSecondCycles);

	computeScores();

	displayHypothesis();
}

void MC3DInteractionHypothesisGenerator::initializeMaps()
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

	mNucleotideStringIdMap.insert(std::pair<std::string, unsigned int>("A", 0));
	mNucleotideStringIdMap.insert(std::pair<std::string, unsigned int>("C", 1));
	mNucleotideStringIdMap.insert(std::pair<std::string, unsigned int>("G", 2));
	mNucleotideStringIdMap.insert(std::pair<std::string, unsigned int>("U", 3));
	mNucleotideIdStringMap.insert(std::pair<unsigned int, std::string>(0, "A"));
	mNucleotideIdStringMap.insert(std::pair<unsigned int, std::string>(1, "C"));
	mNucleotideIdStringMap.insert(std::pair<unsigned int, std::string>(2, "G"));
	mNucleotideIdStringMap.insert(std::pair<unsigned int, std::string>(3, "U"));
}

void MC3DInteractionHypothesisGenerator::version () const
{
	std::cout << mstrAppName << " "; // Nom du logiciel
	std::cout << mfAppVersion << " ";			// Version du logiciel
	std::cout << "(" << __DATE__ << ")";	// Date de la compilation
	std::cout << std::endl;
}

void MC3DInteractionHypothesisGenerator::usage () const
{
	std::cout << "usage: "
		<< " [-hV] -a <cycles from first domain> -b <cycle from second domain> [-1 <t1 file>] [-2 <t2 file>] [-3 <t3 file>] [-4 <t4 file>]"
		<< std::endl;
}


void MC3DInteractionHypothesisGenerator::help () const
{
	std::cout
		<< "This program computes the probability of a junctions between a sets of NCMs." << std::endl
		<< "  -h	print this information" << std::endl
		<< "  -V    print the software version info" << std::endl
		<< "  -a    List of cycles in the form 'id : profile : sequence'" << std::endl
		<< "  -b    List of cycles in the form 'id : profile : sequence' against which to test" << std::endl
		<< "  -1	File containing the term 1 table." << std::endl
		<< "  -2	File containing the term 2 table." << std::endl
		<< "  -3	File containing the term 3 table." << std::endl
		<< "  -4	File containing the term 4 table." << std::endl;
}

void MC3DInteractionHypothesisGenerator::read_options (int argc, char* argv[])
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
		case 'a':
			mstrFirstCycles = optarg;
			break;
		case 'b':
			mstrSecondCycles = optarg;
			break;
		case '1':
			mstrTerm1File = optarg;
			break;
		case '2':
			mstrTerm2File = optarg;
			break;
		case '3':
			mstrTerm3File = optarg;
			break;
		case '4':
			mstrTerm4File = optarg;
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if ((0 == mstrFirstCycles.size()) || (0 == mstrSecondCycles.size()))
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

void MC3DInteractionHypothesisGenerator::read_database()
{
	// Read first term
	readTerm1();

	// Read second term
	readTerm2();

	// Read third term
	readTerm3();

	// Read fourth term
	readTerm4();
}

MC3DInteractionHypothesisGenerator::cycle_key MC3DInteractionHypothesisGenerator::cycleKey(const std::string& astrKey) const
{
	std::string strKey = astrKey;
	annotate::cleanString(strKey, ' ');
	std::vector<std::string> fields = annotate::splitStringFields(strKey, ":");
	return cycle_key(fields[0], fields[1]);
}

std::string MC3DInteractionHypothesisGenerator::cycleKeyToString(const MC3DInteractionHypothesisGenerator::cycle_key& aKey) const
{
	std::ostringstream oss;
	oss << aKey.first << " : " << aKey.second;
	return oss.str();
}

MC3DInteractionHypothesisGenerator::indexed_cycles MC3DInteractionHypothesisGenerator::readIndexedCycleFile(const std::string& astrFileName) const
{
	std::map<std::string, cycle_key> indexedCycleFile;
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
				std::string strIndex = fields[0].c_str();
				std::string strProfile = fields[1];
				std::string strSequence = fields[2];
				cycle_key key(strProfile, strSequence);
				key = forceSymmetric(key);
				indexedCycleFile.insert(std::pair<std::string, cycle_key>(strIndex, key));
			}
		}
	}
	else
	{
		// TODO : This should be an exception
		std::cout << "Error opening file " << astrFileName << std::endl;
	}
	infile.close();
	return indexedCycleFile;
}

void MC3DInteractionHypothesisGenerator::readTerm1()
{
	unsigned int uiNbProfiles = mProfileMap.size();
	mTerm1.resize(uiNbProfiles);
	for(unsigned int i = 0; i < uiNbProfiles; ++ i)
	{
		// By default assume everything is probable
		mTerm1[i].resize(uiNbProfiles, 1.0f);
	}

	if(0 < mstrTerm1File.size())
	{
		std::ifstream infile;
		infile.open(mstrTerm1File.c_str(), std::ios_base::in);
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
					assert(fields.size() == 1 + mProfileMap.size());
					std::string strProfile = fields[0];
					unsigned int uiProfileId = mProfileMap[strProfile]; // TODO : Check this exist
					for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
					{
						mTerm1[uiProfileId][i] = atof(fields[i + 1].c_str());
					}
				}
			}
		}
		else
		{
			// TODO : This should be an exception
			std::cout << "Error opening file " << mstrTerm1File << std::endl;
		}
		infile.close();
	}
}

void MC3DInteractionHypothesisGenerator::readTerm2()
{
	unsigned int uiNbProfiles = mProfileMap.size();
	mTerm2.resize(uiNbProfiles);
	for(unsigned int i = 0; i < uiNbProfiles; ++ i)
	{
		// By default assume everything is probable
		mTerm2[i].resize(uiNbProfiles);
		for(unsigned int j = 0; j < uiNbProfiles; ++ j)
		{
			mTerm2[i][j].resize(mInteractionTypeMap.size(), 1.0f);
		}
	}

	if(0 < mstrTerm2File.size())
	{
		std::ifstream infile;
		infile.open(mstrTerm2File.c_str(), std::ios_base::in);
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
					unsigned int uiProfile1 = mProfileMap[strProfile1];
					unsigned int uiProfile2 = mProfileMap[strProfile2];
					mTerm2[uiProfile1][uiProfile2][uiTypeId] = fFreq;
				}
			}
		}
		else
		{
			// TODO : This should be an exception
			std::cout << "Error opening file " << mstrTerm2File << std::endl;
		}
		infile.close();
	}
}

void MC3DInteractionHypothesisGenerator::readTerm3()
{
	float fDefault = 1.0f; // In case we want to ignore term 3 effect
	if(0 < mstrTerm3File.size())
	{
		fDefault = 0.0f;
	}

	unsigned int uiNbProfiles = mProfileMap.size();
	mTerm3.resize(uiNbProfiles);
	for(unsigned int i = 0; i < uiNbProfiles; ++ i)
	{
		unsigned int uiProfSize1 = profileSize(i);
		mTerm3[i].resize(uiNbProfiles);
		for(unsigned int j = 0; j < uiNbProfiles; ++ j)
		{
			unsigned int uiProfSize2 = profileSize(j);
			mTerm3[i][j].resize(mInteractionTypeMap.size());
			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				mTerm3[i][j][k].resize(uiProfSize1);
				for(unsigned int l = 0; l < uiProfSize1; ++ l)
				{
					mTerm3[i][j][k][l].resize(uiProfSize2, fDefault);
				}
			}
		}
	}

	if(0 < mstrTerm3File.size())
	{
		std::ifstream infile;
		infile.open(mstrTerm3File.c_str(), std::ios_base::in);
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
					unsigned int uiPos1 = atol(fields[3].c_str());
					unsigned int uiPos2 = atol(fields[4].c_str());
					float fFreq = atof(fields[5].c_str());

					unsigned int uiTypeId = interactionTypeStringToId(strInterType);
					unsigned int uiProfile1 = mProfileMap[strProfile1];
					unsigned int uiProfile2 = mProfileMap[strProfile2];
					mTerm3[uiProfile1][uiProfile2][uiTypeId][uiPos1][uiPos2] = fFreq;
				}
			}
		}
		else
		{
			// TODO : This should be an exception
			std::cout << "Error opening file " << mstrTerm3File << std::endl;
		}
		infile.close();
	}
}

void MC3DInteractionHypothesisGenerator::readTerm4()
{
	float fDefault = 1.0f; // In case we want to ignore term 3 effect
	if(0 < mstrTerm4File.size())
	{
		fDefault = 0.0f;
	}

	const unsigned int uiNbNucleotides = 4;
	mTerm4.resize(uiNbNucleotides);
	for(unsigned int i = 0; i < uiNbNucleotides; ++ i)
	{
		mTerm4[i].resize(uiNbNucleotides);
		for(unsigned int j = 0; j < uiNbNucleotides; ++ j)
		{
			mTerm4[i][j].resize(mInteractionTypeMap.size(), fDefault);
		}
	}

	if(0 < mstrTerm4File.size())
	{
		std::ifstream infile;
		infile.open(mstrTerm4File.c_str(), std::ios_base::in);
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
					std::string strNucleotide1 = fields[0];
					std::string strNucleotide2 = fields[1];
					std::string strInterType = fields[2];
					float fFreq = atof(fields[3].c_str());

					unsigned int uiTypeId = interactionTypeStringToId(strInterType);
					unsigned int uiNucId1 = nucleotideId(strNucleotide1);
					unsigned int uiNucId2 = nucleotideId(strNucleotide2);
					mTerm4[uiNucId1][uiNucId2][uiTypeId] = fFreq;
				}
			}
		}
		else
		{
			// TODO : This should be an exception
			std::cout << "Error opening file " << mstrTerm4File << std::endl;
		}
		infile.close();
	}
}

MC3DInteractionHypothesisGenerator::cycle_key MC3DInteractionHypothesisGenerator::forceSymmetric(const cycle_key& aKey) const
{
	cycle_key symKey = aKey;
	std::vector<std::string> fields = annotate::splitStringFields(aKey.first, "_");
	if(2 == fields.size())
	{
		unsigned int uiStrand1 = atol(fields[0].c_str());
		unsigned int uiStrand2 = atol(fields[1].c_str());
		if(uiStrand1 == uiStrand2)
		{
			// May be symmetric, try inversing strands
			std::string strFlippedSeq = flipSequence(uiStrand1, uiStrand2, aKey.second);
			if(strFlippedSeq < aKey.second)
			{
				symKey.second = strFlippedSeq;
			}
		}else if(uiStrand1 > uiStrand2)
		{
			std::ostringstream oss;
			oss << uiStrand2 << "_" << uiStrand1;
			symKey.first = oss.str();
			std::string strFlippedSeq = flipSequence(uiStrand1, uiStrand2, aKey.second);
			symKey.second = strFlippedSeq;
		}
	}
	return symKey;
}

std::string MC3DInteractionHypothesisGenerator::flipSequence(
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

void MC3DInteractionHypothesisGenerator::computeScores()
{
	indexed_cycles::const_iterator it1;
	indexed_cycles::const_iterator it2;
	for(it1 = mIndexedCycleFile1.begin(); it1 != mIndexedCycleFile1.end(); ++ it1)
	{
		// float fProbI1 = freqTertiary(it1->second);
		for(it2 = mIndexedCycleFile2.begin(); it2 != mIndexedCycleFile2.end(); ++ it2)
		{
			unsigned int uiProfile1 = mProfileMap[it1->second.first];
			unsigned int uiProfile2 = mProfileMap[it2->second.first];
			unsigned int uiProfSize1 = profileSize(uiProfile1);
			unsigned int uiProfSize2 = profileSize(uiProfile2);
			float fTerm1 = mTerm1[uiProfile1][uiProfile2];
			std::vector<float>& vTerm2 = mTerm2[uiProfile1][uiProfile2];
			for(unsigned int uiTypeId = 0; uiTypeId < this->mInteractionTypeMap.size(); ++ uiTypeId)
			{
				float fTerm2 = vTerm2[uiTypeId];
				float fTerm1_2 = fTerm1 * fTerm2;
				std::vector<std::vector<float> >& vTerm3 = mTerm3[uiProfile1][uiProfile2][uiTypeId];
				for(unsigned int uiPos1 = 0; uiPos1 < uiProfSize1; ++ uiPos1)
				{
					unsigned int uiNuc1 = nucleotideId(getNucleotide(it1->second, uiPos1));
					std::vector<float>& vTerm3a = vTerm3[uiPos1];
					for(unsigned int uiPos2 = 0; uiPos2 < uiProfSize2; ++ uiPos2)
					{
						unsigned int uiNuc2 = nucleotideId(getNucleotide(it2->second, uiPos2));
						float fTerm3 = vTerm3a[uiPos2];
						float fTerm1_2_3 = fTerm1_2 * fTerm3;

						float fTerm4 = mTerm4[uiNuc1][uiNuc2][uiTypeId];
						float fScore = fTerm4 * fTerm1_2_3;

						InteractionHypothesis hyp(it1->first,it2->first, uiTypeId, uiPos1, uiPos2);
						mScoredHypothesis.insert(std::pair<float, InteractionHypothesis>(fScore, hyp));
					}
				}
			}
		}
	}
}

void MC3DInteractionHypothesisGenerator::displayHypothesis() const
{
	std::multimap<float, InteractionHypothesis>::const_iterator it;
	for(it = mScoredHypothesis.begin(); it != mScoredHypothesis.end(); ++ it)
	{
		std::cout << it->first;
		std::cout << " : " << interactionTypeIdToString(it->second.muiType);
		std::cout << " : " << it->second.muiPos1;
		std::cout << " : " << it->second.muiPos2;
		std::cout << " : " << it->second.mstrCycle1 << " : " << it->second.mstrCycle2;
		std::cout << std::endl;
	}
}

unsigned int MC3DInteractionHypothesisGenerator::interactionTypeStringToId(const std::string& astrType) const
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

std::string MC3DInteractionHypothesisGenerator::interactionTypeIdToString(unsigned int auiTypeId) const
{
	const char* szTypes[] = {
		"cWW", "tWW",   "cWH", "tWH", "cWS", "tWS", "uWR", "uWP",
		"cHW", "tHW",   "cHH", "tHH", "cHS", "tHS", "uHR", "uHP",
		"cSW", "tSW",   "cSH", "tSH", "cSS", "tSS", "uSR", "uSP",
		"uRW", "uRH",   "uRS", "uRR", "uRP",
		"uPW", "uPH",   "uPS", "uPR", "uPP"
		};
	assert(auiTypeId < mInteractionTypeMap.size());
	std::string strType = szTypes[auiTypeId];
	return strType;
}

unsigned int MC3DInteractionHypothesisGenerator::profileSize(unsigned int auiProfileId) const
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

std::string MC3DInteractionHypothesisGenerator::nucleotideString(unsigned int auiNucId) const
{
	std::map<unsigned int, std::string>::const_iterator it;
	it = mNucleotideIdStringMap.find(auiNucId);
	assert(it != mNucleotideIdStringMap.end());
	return it->second;
}

unsigned int MC3DInteractionHypothesisGenerator::nucleotideId(const std::string& astrNucleotide) const
{
	std::map<std::string, unsigned int>::const_iterator it;
	it = mNucleotideStringIdMap.find(astrNucleotide);
	assert(it != mNucleotideStringIdMap.end());
	return it->second;
}

std::string MC3DInteractionHypothesisGenerator::getNucleotide(
	const cycle_key& aCycle,
	unsigned int auiPosition) const
{
	std::vector<std::string> fields = annotate::splitStringFields(aCycle.second, "-");
	return fields[auiPosition];
}

int main(int argc, char* argv[])
{
	MC3DInteractionHypothesisGenerator theApp(argc, argv);
	return EXIT_SUCCESS;
}
