/*
 * ExtractCyclesApp.cc
 *
 *  Created on: Jan 28, 2010
 *      Author: blanchmf
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "ExtractCyclesApp.h"

#include "StringUtil.h"

#include <cstdio>
#include <fstream>

const char* shortopts = "Vhuc:d:";

ExtractCyclesApp::ExtractCyclesApp(int argc, char * argv [])
{
	readOptions(argc, argv);

	// Read the cycle files
	readCyclesFile();
}

void ExtractCyclesApp::version () const
{
	std::cout
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")"
		<< std::endl;
}

void ExtractCyclesApp::usage () const
{
	std::cout
	 	<< "usage: " << PACKAGE
		<< " [-fhuvV] [-d <output directory>] -p <distant pairs file> -c <cycles file>"
		<< std::endl;
}

void ExtractCyclesApp::help () const
{
	std::cout
		<< "This program read cycle structures and return the corresponding residue ids."
		<< std::endl
		<< "	-i directory where the PDB files are stored." << std::endl
		<< "  -d	specify an output directory" << std::endl
		<< "  -f	filter out composite cycles" << std::endl
		<< "  -h	print this help" << std::endl
		<< "  -c	file containing the identified cycles" << std::endl
		<< "  -V	print the software version info" << std::endl;
}

void ExtractCyclesApp::readOptions (int argc, char* argv[])
{
	int c;

	while ((c = getopt (argc, argv, shortopts)) != EOF)
	{
		switch (c)
		{
		case 'V':
        	version ();
			exit (EXIT_SUCCESS);
			break;
		case 'c':
		{
			mstrCyclesFile = optarg;
			break;
		}
		case 'd':
		{
			mstrOutputDirectory = optarg;
		}
		case 'i':
		{
			mstrInputDirectory = optarg;
		}
		case 'h':
			usage ();
			help ();
			exit (EXIT_SUCCESS);
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if(mstrInputDirectory.empty())
	{
		usage ();
		exit (EXIT_FAILURE);
	} else if(mstrCyclesFile.empty())
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

/**
 * readCyclesFile
 * @brief Read the cycles from a file.
 */
void ExtractCyclesApp::readCyclesFile()
{
	std::ifstream infile;
	mCycles.clear();
	infile.open (mstrCyclesFile.c_str(), std::ifstream::in);
	std::string strLine;
	// std::vector<std::string> fields;

	while(std::getline(infile, strLine).good())
	{
		annotate::CycleInfo cInfo = readCycleFileLine(strLine);
		/*
		annotate::cleanString(strLine, ' ');

		fields = annotate::splitStringFields(strLine, ":");

		std::string strPDBFile = fields[0];
		unsigned int uiModel = atol(fields[1].c_str());
		std::string strProfile = fields[2];
		std::string strPredProfile = fields[3];
		std::string strResIds = fields[4];
		std::string strSeq = fields[5];

		annotate::CycleProfile prof(strProfile);
		annotate::CycleProfile predProf(strPredProfile);

		CycleInfo::residue_profile resProfile;
		resProfile = getStrandResidues(strResIds, predProf);

		std::vector<std::string> residues = annotate::splitStringFields(strSeq, "-");
		CycleInfo cInfo(strPDBFile, uiModel, predProf, prof, resProfile, residues);*/

		mCycles.insert(cInfo);
	}

	infile.close();
}

annotate::CycleInfo ExtractCyclesApp::readCycleFileLine(const std::string& astrLine) const
{
	std::string strLine = astrLine;
	annotate::cleanString(strLine, ' ');

	std::vector<std::string> fields = annotate::splitStringFields(strLine, ":");

	std::string strPDBFile = fields[0];
	unsigned int uiModel = atol(fields[1].c_str());
	std::string strFileProfile = fields[2];
	std::string strProfile = fields[3];
	std::string strResIds = fields[4];
	std::string strSeq = fields[5];

	annotate::CycleProfile prof(strProfile);
	annotate::CycleProfile fileProf(strFileProfile);

	annotate::CycleInfo::residue_profile resProfile;
	resProfile = getStrandResidues(strResIds, fileProf);

	std::vector<std::string> residues = annotate::splitStringFields(strSeq, "-");
	return annotate::CycleInfo(strPDBFile, uiModel, fileProf, prof, resProfile, residues);
}

std::vector<std::vector<std::string> > ExtractCyclesApp::getStrandResidues(
	const std::string& aResidues,
	const annotate::CycleProfile& aProfile) const
{
	std::list<std::string> residues = getResidues(aResidues);
	std::vector<std::vector<std::string> > strandResidues;

	std::list<std::string>::const_iterator itRes = residues.begin();
	std::list<unsigned int>::const_iterator itProf;
	for(itProf = aProfile.strandProfile().begin();
		itProf != aProfile.strandProfile().end();
		++ itProf)
	{
		std::vector<std::string> strand;
		unsigned int iRes = 0;
		for(iRes = 0; iRes < *itProf; ++ iRes)
		{
			strand.push_back(*itRes);
			++ itRes;
		}
		strandResidues.push_back(strand);
		strand.clear();
	}
	return strandResidues;
}

std::list<std::string> ExtractCyclesApp::getResidues(
	const std::string& aResidues) const
{
	std::list<std::string> residues;
	std::string strResidues = aResidues;

	annotate::cleanString(strResidues, ' ');
	while(0 < strResidues.size())
	{
		std::string strRes;
		std::size_t sep = strResidues.rfind('-');
		if(sep != std::string::npos)
		{
			strRes = strResidues.substr(sep + 1, strResidues.size() - (sep + 1));
			strResidues.erase(sep);
		}
		else
		{
			strRes = strResidues;
			strResidues.clear();
		}
		residues.push_front(strRes);
	}
	return residues;
}
