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

#include <cassert>
#include <set>
#include <string>
#include <sstream>

bool binary = false;
bool oneModel = false;
const char* shortopts = "Vbf:hlv";

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

unsigned int getModelIndex(const std::string& aFileName)
{
	unsigned int uiModel = 0;
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
    	std::string strModel = filename.substr(index + 1);
    	uiModel = strtol (strModel.c_str(), 0, 10);   	
    }
    
	return uiModel;
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

std::vector<std::list<mccore::ResId> > getResProfile(
	const mccore::GraphModel& aModel,
	std::list<unsigned int>& aProfile)
{
	// Compute the profile
	std::vector<std::list<mccore::ResId> > profile;
	std::list<unsigned int>::const_iterator itProf;
	mccore::GraphModel::const_iterator itRes = aModel.begin();
	for(itProf = aProfile.begin(); itProf != aProfile.end(); ++ itProf)
	{
		std::list<mccore::ResId> strand; 
		for(unsigned int iRes = 0; iRes < *itProf; ++ iRes)
		{
			strand.push_back(itRes->getResId());
			++ itRes;
		}
		profile.push_back(strand);
		strand.clear();
	}
	return profile;
}

std::vector< std::list<mccore::ResId> > orderProfile(
	const std::vector< std::list<mccore::ResId> >& aProfile)
{
	std::vector< std::list<mccore::ResId> > profile;
	if(1 == aProfile.size())
	{
		// Loop, nothing to do
		profile = aProfile;
	}
	else if(2 == aProfile.size())
	{
		if(aProfile[0].front() < aProfile[1].front())
		{
			// Already in order
			profile = aProfile;
		}
		else
		{
			profile.push_back(aProfile[1]);
			profile.push_back(aProfile[0]);			
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

std::string getProfileString(std::vector<std::list<mccore::ResId> >& aProfile)
{
	std::ostringstream oss;
	std::vector<std::list<mccore::ResId> >::const_iterator it;
	for(it = aProfile.begin(); it != aProfile.end(); ++ it)
	{
		oss << it->size();
	}
	return oss.str();
}

std::string getResiduesString(std::vector<std::list<mccore::ResId> >& aProfile)
{
	std::ostringstream oss;
	std::vector<std::list<mccore::ResId> >::const_iterator it;
	for(it = aProfile.begin(); it != aProfile.end(); ++ it)
	{
		if(it != aProfile.begin())
		{
			oss << "-";
		}
		std::list<mccore::ResId>::const_iterator itRes;
		for(itRes = it->begin(); itRes != it->end(); ++ itRes)
		{
			if(itRes != it->begin())
			{
				oss << "-";
			}
			oss << *itRes;
		}
	}
	
	return oss.str();
}

int main (int argc, char *argv[])
{
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
				
				mccore::gOut(0) << getPdbFileName(filename) << " : ";
				mccore::gOut(0) << getModelIndex(filename) << "\t: ";
				std::string strFileProfile = getModelProfile(filename);
				std::list<unsigned int> expectedProfile = getExpectedProfile(strFileProfile);
				std::vector<std::list<mccore::ResId> > resProfile = getResProfile(am, expectedProfile);
				resProfile = orderProfile(resProfile);
				
				mccore::gOut(0) << strFileProfile << "\t: ";
				
				// Output the real profile
				mccore::gOut(0) << getProfileString(resProfile) << "\t: ";
				
				// Output the residues
				mccore::gOut(0) << getResiduesString(resProfile);
				
				mccore::gOut(0) << std::endl;
			}
			delete molecule;
		}
		++optind;
	}
	return EXIT_SUCCESS;	
}
