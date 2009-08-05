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
#include "AnnotateModel.h"
#include "AnnotationCycles.h"

#include "mccore/Binstream.h"
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

#include <set>
#include <string>

bool binary = false;
bool oneModel = false;
unsigned int environment = 0;
unsigned char gucRelationMask = 
	mccore::Relation::adjacent_mask
	| mccore::Relation::pairing_mask 
	| mccore::Relation::stacking_mask 
	| mccore::Relation::backbone_mask;
mccore::ResIdSet residueSelection;
const char* shortopts = "Vbf:hlvm:";

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
		<< " [-bhlvV] [-m <relation mask>] <structure file> ..."
		<< std::endl;
}

void help ()
{
	mccore::gOut (0)
		<< "This program read cycle structures and return the corresponding residue ids." << std::endl
		<< "  -b                read binary files instead of pdb files" << std::endl
		<< "  -m <mask>         annotation mask string: any combination of 'A' (adjacent), 'S' (stacking), 'P' (pairing) and 'B' (backbone). (default: all)" << std::endl
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
		case 'm':
    	{
    		std::string strMask = optarg;
    		gucRelationMask = 0;
    		if(strMask.find("A") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::adjacent_mask;
    		}
    		if(strMask.find("S") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::stacking_mask;
    		}
    		if(strMask.find("P") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::pairing_mask;
    		}
    		if(strMask.find("B") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::backbone_mask;
    		}
			break;
    	}
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
	mccore::Molecule *molecule;
	mccore::ResidueFM rFM;
	annotate::AnnotateModelFM aFM (residueSelection, environment, &rFM);

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
		molecule = new mccore::Molecule (&aFM);
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
		molecule = new mccore::Molecule (&aFM);
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

annotate::Cycle mergeCycles(const std::set<annotate::Cycle>& aCycles)
{
	std::set< annotate::Cycle >::const_iterator it = aCycles.begin();
	annotate::Cycle cycle(*it);
	it ++;
	for(;it != aCycles.end(); ++it)
	{
		cycle = cycle.merge(*it);
	}
	return cycle;	
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
				
				annotate::AnnotateModel &am = (annotate::AnnotateModel&) *molIt;
				am.name(getFilePrefix(filename));
					
				annotate::AnnotationInteractions annInteractions;
				annotate::AnnotationCycles annCycles(0);
				am.addAnnotation(annInteractions);
				am.addAnnotation(annCycles);
				am.annotate (gucRelationMask);
				mccore::gOut(0) << getPdbFileName(filename) << " : ";
				mccore::gOut(0) << getModelIndex(filename) << " : ";
				if(0 < annCycles.getCycles().size())
				{
					annotate::Cycle cycle = mergeCycles(annCycles.getCycles());
					std::list<mccore::ResId>::const_iterator itRes;
					for(itRes = cycle.getResidues().begin(); 
						itRes != cycle.getResidues().end(); 
						++ itRes)
					{
						if(itRes != cycle.getResidues().begin())
						{
							mccore::gOut(0) << "-";
						}
						mccore::gOut(0) << *itRes;
					}
				}
				mccore::gOut(0) << std::endl;				
			}
			delete molecule;
		}
		++optind;
	}
	return EXIT_SUCCESS;	
}
