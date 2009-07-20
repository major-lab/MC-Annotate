//                              -*- Mode: C++ -*- 
// annotate.cc
// Copyright © 2001-06 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Patrick Gendron
// Created On       : Fri May 18 09:38:07 2001
// $Revision: 58 $
// $Id: MC-Annotate.cc 58 2006-11-15 21:09:19Z larosem $


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cerrno>
#include <cstdlib>
#include <string>
#include <sstream>
#include <unistd.h>

#include "mccore/Binstream.h"
#include "mccore/Exception.h"
#include "mccore/Messagestream.h"
#include "mccore/ModelFactoryMethod.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"
#include "mccore/ResidueFactoryMethod.h"
#include "mccore/ResIdSet.h"
#ifdef HAVE_LIBRNAMLC__
#include "mccore/RnamlReader.h"
#endif
#include "mccore/Version.h"

#include "AnnotateModel.h"
#include "AnnotationCycles.h"
#include "AnnotationInteractions.h"
#include "AnnotationLinkers.h"
#include "AnnotationLoops.h"
#include "AnnotationResSecondaryStructures.h"
#include "AnnotationStems.h"
#include "AnnotationTertiaryCycles.h"
#include "AnnotationTertiaryPairs.h"
#include "AnnotationTertiaryStacks.h"
#include "AnnotationTertiaryStructures.h"

using namespace mccore;
using namespace std;
using namespace annotate;

bool binary = false;
unsigned int environment = 0;
bool oneModel = false;
unsigned int modelNumber = 0;  // 1 based vector identifier, 0 means all
std::string gstrCycleOutputDir = "";
std::string gstrStructureOutputDir = "";
unsigned int guiMaxCycleSize = 0;
unsigned int guiMax3DCycleSize = 0;
unsigned char gucRelationMask = 
	mccore::Relation::adjacent_mask
	| mccore::Relation::pairing_mask 
	| mccore::Relation::stacking_mask 
	| mccore::Relation::backbone_mask;
ResIdSet residueSelection;
const char* shortopts = "Vbc:s:z:x:e:f:hlr:vm:";

struct stCycleInfoEntry
{
	std::string strProfile;
	std::string strResIds;
	std::string strSequence;
};

struct stCycleInformation
{
	std::string strSignature;
	std::string strModel;
	std::list<stCycleInfoEntry> connects;
};

void
version ()
{
  mccore::Version mccorev;
  
  gOut (0) << PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << endl
	   << "  using " << mccorev << endl;
}


void
usage ()
{
  gOut (0) << "usage: " << PACKAGE
	   << " [-bhlvV] [-e num] [-z <cycle size>] [-x <3D cycle size>] [-s <structure directory>] [-f <model number>] [-r <residue ids>] [-m <relation mask>] [-c <cycle directory>] <structure file> ..."
	   << endl;
}

void
help ()
{
  gOut (0)
    << "This program annotate structures (and more)." << endl
    << "  -b                read binary files instead of pdb files" << endl
    << "  -c <directory>    directory in which to output tertiary cycle files" << endl
	<< "  -s <directory>    directory in which to output tertiary structure files" << endl
    << "  -m <mask>         annotation mask string: any combination of 'A' (adjacent), 'S' (stacking), 'P' (pairing) and 'B' (backbone). (default: all)" << endl
    << "  -z <cycle size>   maximum size of cycles" << endl
	<< "  -x <3Dcycle size> maximum size of non-adjacent cycles" << endl
    << "  -e num            number of surrounding layers of connected residues to annotate" << endl
    << "  -f model number   model to print" << endl
    << "  -h                print this help" << endl
    << "  -l                be more verbose (log)" << endl
    << "  -r sel            extract these residues from the structure" << endl 
    << "  -v                be verbose" << endl
    << "  -V                print the software version info" << endl;    
}


void
read_options (int argc, char* argv[])
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
	case 'e':
	  {
	    long int tmp;

	    tmp = strtol (optarg, 0, 10);
	    if (ERANGE == errno
		|| EINVAL == errno
		|| 0 > tmp)
	      {
		gErr (0) << PACKAGE << ": invalid environment value." << endl;
		exit (EXIT_FAILURE);
	      }
	    environment = tmp;
	    break;
	  }
	case 'c':
	{
		gstrCycleOutputDir = optarg;
		break;		
	}
	case 's':
	{
		gstrStructureOutputDir = optarg;
		break;		
	}
	case 'f':
	  {
	    long int tmp;

	    tmp = strtol (optarg, 0, 10);
	    if (ERANGE == errno
		|| EINVAL == errno
		|| 0 > tmp)
	      {
		gErr (0) << PACKAGE << ": invalid model value." << endl;
		exit (EXIT_FAILURE);
	      }
	    modelNumber = tmp;
	    oneModel = true;
	    break;
	  }
	 case 'z':
	  {
	    guiMaxCycleSize = atol (optarg);
	    break;
	  }
	  case 'x':
	  {
	    guiMax3DCycleSize = atol (optarg);
	    break;
	  }
        case 'h':
          usage ();
          help ();
          exit (EXIT_SUCCESS);
          break;
        case 'l':
          gErr.setVerboseLevel (gErr.getVerboseLevel () + 1);
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
	case 'r':
	  try
	    {
	      residueSelection.insert (optarg);
	    }
	  catch (IntLibException &e)
	    {
	      gErr (0) << PACKAGE << ": invalid residue selection." << endl;
	      exit (EXIT_FAILURE);
	    }
	  break;
	case 'v':
	  gOut.setVerboseLevel (gOut.getVerboseLevel () + 1);
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


mccore::Molecule*
loadFile (const string &filename)
{
  Molecule *molecule;
  ResidueFM rFM;
  AnnotateModelFM aFM (residueSelection, environment, &rFM);

  molecule = 0;
  if (binary)
    {
      izfBinstream in;

      in.open (filename.c_str ());
      if (in.fail ())
	{
	  gErr (0) << PACKAGE << ": cannot open binary file '" << filename << "'." << endl;
	  return 0;
	}
      molecule = new Molecule (&aFM);
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
	  izfPdbstream in;
	  
	  in.open (filename.c_str ());
	  if (in.fail ())
	    {
	      gErr (0) << PACKAGE << ": cannot open pdb file '" << filename << "'." << endl;
	      return 0;
	    }
	  molecule = new Molecule (&aFM);
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

void dumpCyclesFiles(
	const std::string& aFilePrefix, 
	const AnnotationTertiaryCycles& aAnnotationCycles)
{
	int i = 1;
	std::set<annotate::Cycle>::const_iterator it;
	for(it = aAnnotationCycles.getCycles().begin();
		it != aAnnotationCycles.getCycles().end();
		++ it)
	{		
		std::ostringstream oss;
		if(!gstrCycleOutputDir.empty())
		{
			oss << gstrCycleOutputDir << "/";
		}
		oss << aFilePrefix << "_l" << it->getModel().size() << "_c" << i << ".pdb.gz";
		ozfPdbstream ops;
		ops.open ((oss.str ()).c_str ());
		if (! ops)
		{
			gErr (0) << "Cannot open file " << oss.str() << endl;
			continue;
		}
		ops << it->getModel();
		ops.close ();
		++ i;
	}	
}

void dumpStructuresFiles(
	const std::string& aFilePrefix, 
	const AnnotationTertiaryStructures& aAnnotationStructures)
{
	int i = 1;
	for(std::list<TertiaryStructure>::const_iterator it = aAnnotationStructures.getStructures().begin();
		it != aAnnotationStructures.getStructures().end();
		++ it)
	{		
		std::ostringstream oss;
		if(!gstrStructureOutputDir.empty())
		{
			oss << gstrStructureOutputDir << "/";
		}
		oss << it->name() << "-l" << it->getModel().size() << ".pdb.gz";
		ozfPdbstream ops;
		ops.open ((oss.str ()).c_str ());
		if (! ops)
		{
			gErr (0) << "Cannot open file " << oss.str() << endl;
			continue;
		}
		ops << it->getModel();
		ops.close ();
		++ i;
	}	
}

std::string cycleResidues(const annotate::Cycle& aCycle)
{
	std::ostringstream oss;
	std::list<mccore::ResId>::const_iterator it;
	for(it = aCycle.getResidues().begin(); it != aCycle.getResidues().end(); ++ it)
	{
		if(it != aCycle.getResidues().begin())
		{
			oss << ", ";
		}
		oss << *it;
	}
	return oss.str();
}

std::string cycleProfile(const annotate::Cycle& aCycle)
{
	std::ostringstream oss;
	std::vector<unsigned int>::const_iterator it;
	for(it = aCycle.profile().begin(); it != aCycle.profile().end(); ++ it)
	{
		oss << *it;
	}
	return oss.str();
}

std::string outputCycleInformation(const stCycleInformation& info)
{
	std::ostringstream oss;
	oss << "index:" << std::endl;
	oss << "\t" << info.strSignature << std::endl;
	oss << "info:" << std::endl;
	oss << "\t" << info.strModel << std::endl;
	std::list<stCycleInfoEntry>::const_iterator it;
	for(it = info.connects.begin(); it != info.connects.end(); ++it)
	{
		oss << "\t" << it->strProfile << "\t : ";;
		oss << it->strResIds << ";\t";
		oss << it->strSequence << ";" << std::endl;
	}
	return oss.str();
}

stCycleInformation cycleInformation(
	const annotate::Cycle& aCycle,
	const annotate::AnnotationTertiaryCycles& annot)
{
	stCycleInformation info;
	info.strModel = aCycle.modelName();
	info.strSignature = cycleProfile(aCycle) + " = ";
	
	stCycleInfoEntry cycleEntry;
	cycleEntry.strProfile = cycleProfile(aCycle);
	cycleEntry.strResIds = cycleResidues(aCycle);
	cycleEntry.strSequence = aCycle.getSequence();
	
	info.connects.push_back(cycleEntry);

	std::set<Cycle> connections = annot.getConnections(aCycle);
	std::set<Cycle>::const_iterator it;
	for(it = connections.begin(); it != connections.end(); ++it)
	{
		stCycleInfoEntry entry;
		entry.strProfile = cycleProfile(*it);
		entry.strResIds = cycleResidues(*it);
		entry.strSequence = it->getSequence();
		if(it != connections.begin())
		{
			info.strSignature += " + ";
		}
		info.strSignature += entry.strProfile;
		info.connects.push_back(entry);	
	}

	return info;
}

int
main (int argc, char *argv[])
{
	std::list<stCycleInformation> cyclesInformations;
	std::list<Cycle> non_adjacent_cycles;
	read_options (argc, argv);

	while (optind < argc)
    {
		Molecule *molecule;
		Molecule::iterator molIt;
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
					AnnotateModel &am = (AnnotateModel&) *molIt;
					am.name(getFilePrefix(filename));
					AnnotationInteractions annInteractions;
					AnnotationStems annStems;
					AnnotationLinkers annLinkers;
					AnnotationLoops annLoops;
					AnnotationTertiaryPairs annTertiaryPairs;
					AnnotationTertiaryStacks annTertiaryStacks;
					AnnotationCycles annCycles(guiMaxCycleSize);
					AnnotationTertiaryCycles annTertiaryCycles(guiMax3DCycleSize);
					AnnotationResSecondaryStructures annResSecondaryStructures;
					AnnotationTertiaryStructures annTertiaryStructures;
		  
		  			am.addAnnotation(annInteractions);
					am.addAnnotation(annStems);
					am.addAnnotation(annLinkers);
					am.addAnnotation(annLoops);
					am.addAnnotation(annTertiaryPairs);
					am.addAnnotation(annTertiaryStacks);
					am.addAnnotation(annCycles);
					am.addAnnotation(annTertiaryCycles);
					am.addAnnotation(annResSecondaryStructures);
					am.addAnnotation(annTertiaryStructures);
					
					gOut (0) << "Annotating Model ------------------------------------------------" << endl;										
					gOut (0) << filename << std::endl;
					am.annotate (gucRelationMask);
					gOut(0) << am;
					
					if(!gstrCycleOutputDir.empty())
					{
						dumpCyclesFiles(getFilePrefix(filename), annTertiaryCycles);
					}
					
					if(!gstrStructureOutputDir.empty())
					{
						dumpStructuresFiles(getFilePrefix(filename), annTertiaryStructures);
					}
										
					for(std::set<Cycle>::const_iterator itCycle = annTertiaryCycles.getCycles().begin();
						itCycle != annTertiaryCycles.getCycles().end();
						++ itCycle)
					{
						non_adjacent_cycles.push_back(*itCycle);
						cyclesInformations.push_back(cycleInformation(*itCycle, annTertiaryCycles));
					}
					
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
	
	mccore::gOut (0) << "------------------------------------------------------------" << std::endl;
	std::list<stCycleInformation>::const_iterator itInfo;
	for(itInfo = cyclesInformations.begin(); 
		itInfo != cyclesInformations.end(); 
		++ itInfo)
	{
		mccore::gOut(0) << "------------------------------" << std::endl;
		mccore::gOut(0) << outputCycleInformation(*itInfo);
		mccore::gOut(0) << std::endl;
	}
	
	return EXIT_SUCCESS;	
}
