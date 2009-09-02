//                              -*- Mode: C++ -*-
// mcsc2resids.cc
// Copyright © 2001-09 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Jul 29 10:35:00 2009

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mcsc2resids.h"

// libmcannotate.a
#include "CycleProfile.h"
#include "StringTable.h"
#include "StringUtil.h"

#include "mccore/Binstream.h"
#include "mccore/GraphModel.h"
#include "mccore/Messagestream.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/Relation.h"
#include "mccore/ResidueFactoryMethod.h"
#include "mccore/ResIdSet.h"
#ifdef HAVE_LIBRNAMLC__
#include "mccore/RnamlReader.h"
#endif
#include "mccore/Version.h"

#include <algorithm>
#include <cassert>
#include <set>
#include <string>
#include <sstream>

bool binary = false;
bool oneModel = false;
const char* shortopts = "Vbhlv";

typedef std::list<const mccore::Residue*> residue_strand;
typedef std::vector<residue_strand > residue_profile;

struct stCycleInfo
{
	std::string strFile;
	std::string strModelId;
	std::string strFileProfile;
	std::string strProfile;
	std::string strResIds;
	std::string strNucleotides;
};

// PROTOTYPES ------------------------------------------------------------------
std::string getResiduesString(const residue_profile& aProfile);

void version ()
{
	mccore::Version mccorev;
	mccore::gOut (0) << PACKAGE << " " << VERSION << " (" << __DATE__ << ")"
		<< std::endl
		<< "  using " << mccorev << endl;
}


void usage ()
{
	mccore::gOut (0) << "usage: " << PACKAGE
		<< " [-bhlvV] <structure file> ..."
		<< std::endl;
}

void help ()
{
	mccore::gOut (0)
		<< "This program read cycle structures and return the corresponding residue ids." << std::endl
		<< "  -b                read binary files instead of pdb files" << std::endl
		<< "  -h                print this help" << std::endl
		<< "  -l                be more verbose (log)" << std::endl
		<< "  -v                be verbose" << std::endl
		<< "  -V                print the software version info" << std::endl;
}


void read_options (int argc, char* argv[])
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
		case 'b':
			binary = true;
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

	if (argc - optind < 1)
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}


mccore::Molecule* loadFile (const std::string &filename)
{
	mccore::ResIdSet residueSelection;
	mccore::Molecule *molecule;
	mccore::ResidueFM rFM;

	molecule = 0;
	if (binary)
	{
		mccore::izfBinstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			mccore::gErr (0) << PACKAGE << ": cannot open binary file '" << filename << "'." << endl;
			return 0;
		}
		molecule = new mccore::Molecule ();
		in >> *molecule;
		in.close ();
	}
	else
	{
#ifdef HAVE_LIBRNAMLC__
		RnamlReader reader (filename.c_str (), &aFM);

		if (0 == (molecule = reader.read ()))
		{
#endif
		mccore::izfPdbstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			mccore::gErr (0) << PACKAGE << ": cannot open pdb file '"
				<< filename << "'." << std::endl;
			return 0;
		}
		molecule = new mccore::Molecule ();
		in >> *molecule;
		in.close ();
#ifdef HAVE_LIBRNAMLC__
		}
#endif
	}
	return molecule;
}

std::string getModelIndex(const std::string& aFileName)
{
	std::string strModel;
	std::string::size_type index;
	std::string filename = aFileName;
	if (std::string::npos != (index = filename.rfind ("/")))
    {
		filename.erase (0, index + 1);
    }
	if (string::npos != (index = filename.find (".")))
    {
		filename.erase (index, filename.size ());
    }
    if(std::string::npos != (index = filename.rfind("_")))
    {
    	strModel = filename.substr(index + 1);
    }

	return strModel;
}

bool isProfileParallel(const std::string& aProfile)
{
	bool bParallel = false;
	std::string aDelim = "_";
	std::vector<std::string> fields = annotate::splitStringFields(aProfile, aDelim);

	if(0 < fields.size())
	{
		bParallel = (fields.back() == "p");
	}

	return bParallel;
}

std::string getModelProfile(const std::string& aFileName)
{
	std::string::size_type index;
	std::string strProfile = aFileName;
	if (std::string::npos != (index = strProfile.rfind ("/")))
    {
		strProfile.erase (0, index + 1);
    }
	if (string::npos != (index = strProfile.find (".")))
    {
		strProfile.erase (index, strProfile.size ());
    }
    if(std::string::npos != (index = strProfile.find("_")))
    {
    	strProfile.erase (0, index + 1);
    }
    if(std::string::npos != (index = strProfile.find("_model")))
    {
    	strProfile.erase (index, strProfile.size());
    }

	return strProfile;
}

std::string getFilePrefix(const std::string& aFileName)
{
	std::string::size_type index;
	std::string filename = aFileName;
	if (std::string::npos != (index = filename.rfind ("/")))
    {
		filename.erase (0, index + 1);
    }
	if (string::npos != (index = filename.find (".")))
    {
		filename.erase (index, filename.size ());
    }
	return filename;
}

std::string getPdbFileName(const std::string& aFileName)
{
	std::string::size_type index;
	std::string filename = getFilePrefix(aFileName);
    if(std::string::npos != (index = filename.find("_")))
    {
    	filename.erase(index, filename.size());
    }
	return filename;
}

residue_profile getResProfile(
	const mccore::GraphModel& aModel,
	annotate::CycleProfile& aProfile)
{
	// Compute the profile
	residue_profile profile;
	std::list<unsigned int>::const_iterator itProf;
	mccore::GraphModel::const_iterator itRes = aModel.begin();
	for(itProf = aProfile.strandProfile().begin();
		itProf != aProfile.strandProfile().end();
		++ itProf)
	{
		residue_strand strand;
		for(unsigned int iRes = 0; iRes < *itProf; ++ iRes)
		{
			strand.push_back(&(*itRes));
			++ itRes;
		}
		profile.push_back(strand);
		strand.clear();
	}
	return profile;
}

residue_profile orderProfile(const residue_profile& aProfile)
{
	residue_profile profile;
	if(1 == aProfile.size())
	{
		// Loop, nothing to do
		profile = aProfile;
	}
	else if(2 == aProfile.size())
	{
		if(aProfile[1].front()->getResId() < aProfile[0].front()->getResId())
		{
			profile.push_back(aProfile[1]);
			profile.push_back(aProfile[0]);
		}
		else
		{
			// Already in order
			profile = aProfile;
		}
	}else
	{
		mccore::gOut(0) << "Only single and dual stranded cycle supported" << std::endl;
		exit(0);
	}

	return profile;
}

std::string getProfileString(
	residue_profile& aProfile,
	const std::string& astrProfile)
{
	std::ostringstream oss;
	residue_profile::const_iterator it;
	for(it = aProfile.begin(); it != aProfile.end(); ++ it)
	{
		if(it != aProfile.begin())
		{
			oss << "_";
		}
		oss << it->size();
	}

	std::vector<std::string> fields = annotate::splitStringFields(astrProfile, "_");
	if(1 == fields.back().size() && std::isalpha(fields.back()[0]))
	{
		oss << "_" << fields.back();
	}
	return oss.str();
}

std::string getResiduesString(const residue_profile& aProfile)
{
	std::ostringstream oss;
	residue_profile::const_iterator it;
	for(it = aProfile.begin(); it != aProfile.end(); ++ it)
	{
		if(it != aProfile.begin())
		{
			oss << "-";
		}
		residue_strand::const_iterator itRes;
		for(itRes = it->begin(); itRes != it->end(); ++ itRes)
		{
			if(itRes != it->begin())
			{
				oss << "-";
			}
			oss << (*itRes)->getResId();
		}
	}

	return oss.str();
}

std::string getNucleotideString(residue_profile& aProfile)
{
	std::ostringstream oss;
	residue_profile::const_iterator it;
	for(it = aProfile.begin(); it != aProfile.end(); ++ it)
	{
		if(it != aProfile.begin())
		{
			oss << "-";
		}
		residue_strand::const_iterator itRes;
		for(itRes = it->begin(); itRes != it->end(); ++ itRes)
		{
			if(itRes != it->begin())
			{
				oss << "-";
			}
			oss << mccore::Pdbstream::stringifyResidueType((*itRes)->getType());
		}
	}

	return oss.str();
}

void makeParallel(residue_profile& aProfile)
{
	if(2 == aProfile.size())
	{
		std::reverse(aProfile[1].begin(), aProfile[1].end());
	}
}

int main (int argc, char *argv[])
{
	annotate::StringTable stringTable(6);

	read_options (argc, argv);

	while (optind < argc)
    {
		mccore::Molecule *molecule;
		mccore::Molecule::iterator molIt;
		std::string filename = (std::string) argv[optind];

		molecule = loadFile (filename);
		if (0 != molecule)
		{
			for (molIt = molecule->begin (); molecule->end () != molIt; ++molIt)
			{
				mccore::GraphModel &am = (mccore::GraphModel&) *molIt;

				std::string strFileProfile = getModelProfile(filename);
				annotate::CycleProfile expectedProfile(strFileProfile);
				residue_profile resProfile = getResProfile(am, expectedProfile);
				bool bIsParallel = isProfileParallel(strFileProfile);
				if(bIsParallel)
				{
					makeParallel(resProfile);
				}
				resProfile = orderProfile(resProfile);

				std::vector<string>& tableRow = stringTable.addRow();
				tableRow[0] = getPdbFileName(filename);
				tableRow[1] = getModelIndex(filename);
				tableRow[2] = strFileProfile;
				tableRow[3] = getProfileString(resProfile, strFileProfile);
				tableRow[4] = getResiduesString(resProfile);
				tableRow[5] = getNucleotideString(resProfile);
			}
			delete molecule;
		}
		++optind;
	}

	mccore::gOut(0) << stringTable.toString(" : ");

	return EXIT_SUCCESS;
}
