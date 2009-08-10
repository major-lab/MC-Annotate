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

std::string getModelProfile(const std::string& aFileName)
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
    if(std::string::npos != (index = filename.find("_")))
    {
    	filename.erase (0, index + 1);
    }
    if(std::string::npos != (index = filename.find("_")))
    {
    	filename.erase (index, filename.size());
    }
    
	return filename;
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
	std::list<unsigned int>& aProfile)
{
	// Compute the profile
	residue_profile profile;
	std::list<unsigned int>::const_iterator itProf;
	mccore::GraphModel::const_iterator itRes = aModel.begin();
	for(itProf = aProfile.begin(); itProf != aProfile.end(); ++ itProf)
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
		if(aProfile[1].front() < aProfile[0].front())
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

std::list<unsigned int> getExpectedProfile(const std::string& astrProfile)
{
	std::list<unsigned int> profile;
	std::string::const_iterator it;
	for(it = astrProfile.begin(); it != astrProfile.end(); ++ it)
	{
		std::string strNumber(1, *it);
		profile.push_back(atol (strNumber.c_str()));
	}	
	return profile;
}

std::string getProfileString(residue_profile& aProfile)
{
	std::ostringstream oss;
	residue_profile::const_iterator it;
	for(it = aProfile.begin(); it != aProfile.end(); ++ it)
	{
		oss << it->size();
	}
	return oss.str();
}

std::string getResiduesString(residue_profile& aProfile)
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

int main (int argc, char *argv[])
{
	std::vector<std::size_t> colSizes;
	colSizes.resize(6, 0);
	std::list<stCycleInfo> infos;
	
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
				std::list<unsigned int> expectedProfile = getExpectedProfile(strFileProfile);
				residue_profile resProfile = getResProfile(am, expectedProfile);
				resProfile = orderProfile(resProfile);
				
				stCycleInfo info;
				info.strFile = getPdbFileName(filename);
				colSizes[0] = std::max(info.strFile.size(), colSizes[0]);
				info.strModelId = getModelIndex(filename);
				colSizes[1] = std::max(info.strModelId.size(), colSizes[1]);
				info.strFileProfile = getModelProfile(filename);
				colSizes[2] = std::max(info.strFileProfile.size(), colSizes[2]);
				info.strProfile =  getProfileString(resProfile);
				colSizes[3] = std::max(info.strProfile.size(), colSizes[3]);
				info.strResIds = getResiduesString(resProfile);
				colSizes[4] = std::max(info.strResIds.size(), colSizes[4]);
				info.strNucleotides = getNucleotideString(resProfile);
				colSizes[5] = std::max(info.strNucleotides.size(), colSizes[5]);
				infos.push_back(info);
			}
			delete molecule;
		}
		++optind;
	}
	
	std::list<stCycleInfo>::iterator infoIt;
	for(infoIt = infos.begin(); infoIt != infos.end(); ++infoIt)
	{
		infoIt->strFile.resize(colSizes[0], ' ');
		infoIt->strModelId.resize(colSizes[1], ' ');
		infoIt->strFileProfile.resize(colSizes[2], ' ');
		infoIt->strProfile.resize(colSizes[3], ' ');
		infoIt->strResIds.resize(colSizes[4], ' ');
		infoIt->strNucleotides.resize(colSizes[5], ' ');
	}
	for(infoIt = infos.begin(); infoIt != infos.end(); ++infoIt)
	{
		mccore::gOut(0) << infoIt->strFile << " : ";
		mccore::gOut(0) << infoIt->strModelId << " : ";
		mccore::gOut(0) << infoIt->strFileProfile << " : ";
		mccore::gOut(0) << infoIt->strProfile << " : ";
		mccore::gOut(0) << infoIt->strResIds << " : ";
		mccore::gOut(0) << infoIt->strNucleotides;
		mccore::gOut(0) << std::endl;
	}
	
	return EXIT_SUCCESS;	
}
