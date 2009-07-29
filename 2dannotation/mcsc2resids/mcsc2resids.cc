//                              -*- Mode: C++ -*- 
// mcsc2resids.cc
// Copyright © 2001-09 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Jul 29 10:35:00 2009


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
#include "AnnotateModel.h"
#include "AnnotationCycles.h"

bool binary = false;
bool oneModel = false;
unsigned int environment = 0;
unsigned int modelNumber = 0;  // 1 based vector identifier, 0 means all
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
		<< " [-bhlvV] [-f <model number>] [-m <relation mask>] <structure file> ..."
		<< endl;
}

void help ()
{
	mccore::gOut (0)
		<< "This program annotate structures (and more)." << endl
		<< "  -b                read binary files instead of pdb files" << endl
		<< "  -m <mask>         annotation mask string: any combination of 'A' (adjacent), 'S' (stacking), 'P' (pairing) and 'B' (backbone). (default: all)" << endl
		<< "  -f model number   model to print" << endl
		<< "  -h                print this help" << endl
		<< "  -l                be more verbose (log)" << endl
		<< "  -v                be verbose" << endl
		<< "  -V                print the software version info" << endl;    
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
		case 'f':
		{
			long int tmp;

			tmp = strtol (optarg, 0, 10);
			if (ERANGE == errno	|| EINVAL == errno || 0 > tmp)
			{
				mccore::gErr (0) << PACKAGE << ": invalid model value." << std::endl;
				exit (EXIT_FAILURE);
			}
			modelNumber = tmp;
			oneModel = true;
			break;
		}
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


mccore::Molecule* loadFile (const string &filename)
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
				if (0 != modelNumber)
				{
					--modelNumber;
				}
				else
				{
					annotate::AnnotateModel &am = (annotate::AnnotateModel&) *molIt;
					am.name(getFilePrefix(filename));
					annotate::AnnotationCycles annCycles(0);
					am.addAnnotation(annCycles);														
					mccore::gOut (0) << filename << std::endl;
					am.annotate (gucRelationMask);
					mccore::gOut(0) << am;
					
					if (oneModel)
					{
						break;
					}
				}
			}
			delete molecule;
		}
		++optind;
	}
	return EXIT_SUCCESS;	
}
