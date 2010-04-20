/*
 * RNAGroups.cc
 *
 *  Created on: Mar 27, 2010
 *      Author: blanchmf
 */

#include "RNAGroups.h"
#include <list>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

#include <mccore/Binstream.h>
#include <mccore/Messagestream.h>
#include <mccore/Molecule.h>
#include <mccore/Pdbstream.h>
#include <mccore/ResidueFactoryMethod.h>

#include "AnnotateModel.h"
#include "StringTable.h"
#include "StringUtil.h"

// Command line options
static const char* gszShortopts = "Vhi:p:";

/**
 * Constructor
 */
RNAGroups::RNAGroups(int argc, char * argv [])
: mfAppVersion(1.0f),
mstrAppName("RNAGroups")
{
	read_options (argc, argv);
	readRNASelectOutput();
	std::cout << toString();
	for(unsigned int i = 0; i < mGroups.size(); ++ i)
	{
		validateGroup(i);
	}
}

/**
 * Display information on software version.
 */
void RNAGroups::version () const
{
	std::cout << mstrAppName << " ";	// Software name
	std::cout << mfAppVersion << " ";	// Software version
	std::cout << "(" << __DATE__ << ")";// Compilation date
	std::cout << std::endl;
}

/**
 * Display information on command line usage for the software.
 */
void RNAGroups::usage () const
{
	std::cout << "usage: [-hV] [-p <PDB Database Path>] -i <rnaselect output file>";
	std::cout << std::endl;
}


/**
 * Display detailed help on the parameters taken by the software.
 */
void RNAGroups::help () const
{
	std::cout
		<< "This program extract the groups from rnaselect output." << std::endl
		<< "  -h                print this information" << std::endl
		<< "  -V                print the software version info" << std::endl
		<< "  -i                rnaselect output file" << std::endl
		<< "  -p                PDB Database path" << std::endl;
}

/**
 * Read the parameters from the command line.  The parameters are kept in class
 * members.
 */
void RNAGroups::read_options (int argc, char* argv[])
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
			mstrRNASelectFile = optarg;
			break;
		case 'p':
			mstrMoleculePath = optarg;
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if(0 == mstrRNASelectFile.size())
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

/*
 * @brief Read the RNA Select output file for groups.
 */
void RNAGroups::readRNASelectOutput()
{
	std::ifstream infile;
	infile.open(mstrRNASelectFile.c_str(), std::ios_base::in);
	mGroups.clear();
	if(infile.good())
	{
		std::string strLine;
		while(std::getline(infile, strLine).good())
		{
			annotate::cleanString(strLine, ' ');
			std::vector<std::string> fields = annotate::splitStringFields(strLine, ":");
			if(0 < fields.size() && fields[0] == "Group")
			{
				// We've found a group, create an entry for it
				unsigned int uiGroup = atol(fields[1].c_str());
				std::list<std::string> content = readGroupContent(infile);
				mGroups[uiGroup] = content;
				std::list<std::string>::const_iterator it;
				for(it = content.begin(); it != content.end(); ++ it)
				{
					mFileToGroup[*it] = uiGroup;
				}
			}
		}
	}
	else
	{
		// TODO : This should be an exception
		std::cout << "Error opening file " << mstrRNASelectFile << std::endl;
	}
	infile.close();
}

/*
 * @brief Read the content of a single group from the file.
 */
std::list<std::string> RNAGroups::readGroupContent(std::ifstream& aInputFile)
{
	assert(aInputFile.good());
	std::string strLine;
	bool bGetLine = std::getline(aInputFile, strLine).good();
	assert(bGetLine);
	annotate::cleanString(strLine, '\t');
	annotate::cleanStringStart(strLine, ' ');
	std::vector<std::string> fields = annotate::splitStringFields(strLine, " ");
	std::list<std::string> content;
	for(std::vector<std::string>::const_iterator it = fields.begin(); it != fields.end(); ++ it)
	{
		content.push_back(*it);
	}
	return content;
}

/*
 * Convert the content of the group to a ':' separated string, starting with the
 * group name.
 */
std::string RNAGroups::toString() const
{
	std::ostringstream oss;
	std::map<unsigned int, std::list<std::string> >::const_iterator itGroup;
	for(itGroup = mGroups.begin(); itGroup != mGroups.end(); ++ itGroup)
	{
		std::cout << itGroup->first << " : ";
		std::list<std::string>::const_iterator itFile;
		for(itFile = itGroup->second.begin(); itFile != itGroup->second.end(); ++ itFile)
		{
			if(itFile != itGroup->second.begin())
			{
				std::cout<< " : ";
			}
			std::cout << *itFile;
		}
		std::cout << std::endl;
	}
	return oss.str();
}

mccore::Molecule RNAGroups::loadMolecule(const string &filename, bool abBinary)
{
	mccore::ResidueFM rFM;
	annotate::AnnotateModelFM aFM (mccore::ResIdSet(), 0, &rFM);
	mccore::Molecule molecule(&aFM);

	if (abBinary)
	{
		mccore::izfBinstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			std::ostringstream oss;
			oss << PACKAGE << ": cannot open binary file '" << filename << "'.";
			throw mccore::FileNotFoundException(oss.str(), __FILE__, __LINE__);
		}
		in >> molecule;
		in.close ();
	}
	else
	{
		mccore::izfPdbstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			std::ostringstream oss;
			oss << PACKAGE << ": cannot open pdb file '" << filename << "'.";
			throw mccore::FileNotFoundException(oss.str(), __FILE__, __LINE__);
		}
		in >> molecule;
		in.close ();
	}
	return molecule;
}


std::list<mccore::GraphModel> RNAGroups::getRNAModels(const mccore::Molecule& aMolecule) const
{
	std::list<mccore::GraphModel> models;
	mccore::Molecule::const_iterator it;
	for(it = aMolecule.begin(); it != aMolecule.end(); ++ it)
	{
		mccore::GraphModel model = *it;
		model.keepRNA();
		if(0 < model.size())
		{
			models.push_back(model);
		}
	}
	return models;
}

std::list<mccore::Residue> RNAGroups::getResidues(const mccore::GraphModel& aModel) const
{
	assert(0 < aModel.size());
	std::list<mccore::Residue> residues;
	mccore::GraphModel::const_iterator it;
	for(it = aModel.begin(); it != aModel.end(); ++it)
	{
		assert(it->getType()->isRNA());
		residues.push_back(*it);
	}
	return residues;
}

void RNAGroups::validateGroup(unsigned int auiGroup)
{
	std::cout << "validating group : " << auiGroup << std::endl;
	std::list<std::string>& groupFiles = mGroups[auiGroup];
	std::list<std::string>::const_iterator it;
	for(it = groupFiles.begin(); it != groupFiles.end(); ++ it)
	{
		std::cout << "loading molecule : " << *it << std::endl;
		mccore::Molecule mol = loadMolecule(mstrMoleculePath + *it + ".pdb.gz", false);
		std::list<mccore::GraphModel> models = getRNAModels(mol);
		if(1 < models.size())
		{
			std::cout << *it << " contains multiple RNA chains!" << std::endl;
		}
		else
		{
			std::cout << ".";
		}
	}
}

int main(int argc, char* argv[])
{
	RNAGroups theApp(argc, argv);
	return EXIT_SUCCESS;
}
