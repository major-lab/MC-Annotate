//                              -*- Mode: C++ -*- 
// mcssc.cc
// Copyright © 2009 Laboratoire de Biologie Informatique et Théorique.
//                  Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Aug 18 11:38:07 2009
// $Revision: 1 $
// $Id: mcssc.cc 1 2009-09-18 11:38:07 blanchmf $



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cerrno>
#include <cstdlib>
#include <string>
#include <sstream>
#include <unistd.h>
#include <cassert>

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
#include "AnnotationChains.h"
#include "AnnotationInteractions.h"
#include "AnnotationLinkers.h"
#include "AnnotationLoops.h"
#include "AnnotationStems.h"
#include "BaseLink.h"
#include "BasePair.h"
#include "Cycle.h"

bool binary = false;
bool oneModel = false;
unsigned int guiModelNumber = 0;  // 1 based vector identifier, 0 means all
unsigned char gucRelationMask = 
	mccore::Relation::adjacent_mask
	| mccore::Relation::pairing_mask 
	| mccore::Relation::stacking_mask 
	| mccore::Relation::backbone_mask;
const char* shortopts = "Vbf:hlvm:";

void version ()
{
	mccore::Version mccorev;

	mccore::gOut (0) << PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
	   << "  using " << mccorev << std::endl;
}


void usage ()
{
	mccore::gOut (0) << "usage: " << PACKAGE
	   << " [-bhlvV] [-f <model number>] [-m <relation mask>] <structure file> ..."
	   << std::endl;
}

void help ()
{
	mccore::gOut(0)
		<< "This program annotate secondary structures and produces output the cycles found." << endl
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
			if (ERANGE == errno || EINVAL == errno || 0 > tmp)
			{
				mccore::gErr (0) << PACKAGE << ": invalid model value." << std::endl;
				exit (EXIT_FAILURE);
			}
			guiModelNumber = tmp;
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
			mccore::gOut.setVerboseLevel (gOut.getVerboseLevel () + 1);
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
	annotate::AnnotateModelFM aFM (residueSelection, 0, &rFM);

	molecule = 0;
	if (binary)
	{
		mccore::izfBinstream in;

		in.open (filename.c_str());
		if (in.fail ())
		{
			mccore::gErr(0) << PACKAGE << ": cannot open binary file '";
			mccore::gErr(0) << filename << "'." << endl;
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
			izfPdbstream in;
	
			in.open (filename.c_str ());
			if (in.fail ())
			{
				mccore::gErr(0) << PACKAGE << ": cannot open pdb file '";
				mccore::gErr(0) << filename << "'." << std::endl;
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

std::string getPdbFileName(const std::string& aFileName)
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

void clearInteractionsSet(annotate::Cycle::interactions_set& aSet)
{
	// Clean up the allocated memory
	annotate::Cycle::interactions_set::iterator it;
	for(it = aSet.begin(); it != aSet.end(); ++ it)
	{
		delete *it;	
	}
	aSet.clear();
}

std::set<annotate::Cycle> computeStemCycles(
	const mccore::GraphModel& aModel, 
	const annotate::Stem& aStem)
{
	std::set<annotate::Cycle> cycles;
	
	std::vector< annotate::BasePair >::const_iterator it;
	std::vector< annotate::BasePair >::const_iterator itPrev;
	for(it = aStem.basePairs().begin();	it != aStem.basePairs().end();	++ it)
	{
		if(it != aStem.basePairs().begin())
		{
			annotate::Cycle::interactions_set interactions;
			
			// First pair
			interactions.insert(new annotate::BasePair(
				itPrev->first, itPrev->fResId, 
				itPrev->second, itPrev->rResId)); 
				
			interactions.insert(new annotate::BasePair(
				it->first, it->fResId, 
				it->second, it->rResId)); 
				
			interactions.insert(new annotate::BaseLink(
				itPrev->first, itPrev->fResId, 
				it->first, it->fResId));
				
			interactions.insert(new annotate::BaseLink(
				it->second, it->rResId, 
				itPrev->second, itPrev->rResId));
			
			// Insert the residues in the model
			annotate::Cycle cycle(interactions);
			
			clearInteractionsSet(interactions);
			
			cycles.insert(cycle);
		}
		itPrev = it;
	}
	return cycles;
}

annotate::Cycle computeLoopCycle(
	const annotate::AnnotateModel &aModel, 
	const annotate::Loop& aLoop)
{
	assert(0 < aLoop.linkers().size());

	std::vector< annotate::Linker >::const_iterator it;
	
	annotate::Cycle::interactions_set interactions;
	
	std::pair<std::set<annotate::BaseLink>, std::set<annotate::BasePair> > interPair;
	interPair = aLoop.getInteractions();
	
	std::set<annotate::BaseLink>::const_iterator itLink;
	for(itLink = interPair.first.begin(); 
		itLink != interPair.first.end(); 
		++ itLink)
	{
		interactions.insert(new annotate::BaseLink(
			itLink->first, itLink->fResId, 
			itLink->second, itLink->rResId));
	}
	
	std::set<annotate::BasePair>::const_iterator itPair;
	for(itPair = interPair.second.begin(); 
		itPair != interPair.second.end(); 
		++ itPair)
	{
		interactions.insert(new annotate::BasePair(
			itPair->first, itPair->fResId, 
			itPair->second, itPair->rResId));
	}
	
	// Insert the residues in the model
	assert(0 < interactions.size());
	annotate::Cycle cycle(interactions);
	
	return cycle;
}

std::set<annotate::Cycle> computeSecondaryStructureCycles(
	const annotate::AnnotateModel &aModel)
{
	std::set<annotate::Cycle> cycles;
	const annotate::AnnotationStems* pStems = NULL;
	const annotate::AnnotationLoops* pLoops = NULL;
	
	pStems = aModel.getAnnotation<annotate::AnnotationStems>();
	pLoops = aModel.getAnnotation<annotate::AnnotationLoops>();
	
	std::vector<annotate::Stem>::const_iterator itStem;
	for(itStem = pStems->getStems().begin(); 
		itStem != pStems->getStems().end(); 
		++ itStem)
	{
		std::set<annotate::Cycle> stemCycles = computeStemCycles(aModel, *itStem);
		cycles.insert(stemCycles.begin(), stemCycles.end());
	}
	
	std::vector<annotate::Loop>::const_iterator itLoop;
	for(itLoop = pLoops->getLoops().begin(); 
		itLoop != pLoops->getLoops().end(); 
		++ itLoop)
	{
		annotate::Cycle cycle = computeLoopCycle(aModel, *itLoop);
		cycles.insert(cycle);
	}
	
	return cycles;
}

std::string cycleProfileStrandString(const annotate::Cycle& aCycle)
{
	std::ostringstream oss;
	
	// Describe the profile
	std::vector<unsigned int>::const_iterator it = aCycle.profile().begin();
	for(; it != aCycle.profile().end(); ++ it)
	{
		if(it != aCycle.profile().begin())
		{
			oss << "_";
		}
		oss << *it;
	}
	return oss.str();
}

std::string getResiduesString(
	const annotate::AnnotateModel &aModel,
	const annotate::Cycle& aCycle)
{
	std::ostringstream oss;
	std::list<mccore::ResId>::const_iterator it;
	for(it = aCycle.residues().begin(); it != aCycle.residues().end(); ++ it)
	{
		if(it != aCycle.residues().begin())
		{
			oss << "-";
		}
		GraphModel::const_iterator itRes = aModel.find(*it);
		assert(itRes != aModel.end());
		oss << mccore::Pdbstream::stringifyResidueType(itRes->getType());
	}
	
	return oss.str();
}

std::string cycleProfileString(const annotate::Cycle& aCycle)
{
	std::ostringstream oss;
	
	// Describe the profile
	annotate::Cycle::enType eCycleType = aCycle.getType();
	switch(eCycleType)
	{
		case annotate::Cycle::eLOOSE:
			oss << aCycle.residues().size();
			oss << "_L";
			break;
		case annotate::Cycle::e2STRANDS_PARALLEL:
			oss << cycleProfileStrandString(aCycle) << "_p";
			break;
		case annotate::Cycle::eMULTIBRANCH:
			oss << aCycle.residues().size();
			oss << "_M";
			break;
		default:
			oss << cycleProfileStrandString(aCycle);
			break;		
	}
	
	return oss.str();
}

std::string cycleString(
	const annotate::AnnotateModel &aModel, 
	const annotate::Cycle& aCycle)
{
	std::ostringstream oss;
	
	// Describe the profile
	oss << cycleProfileString(aCycle) << " : ";
	oss << cycleProfileString(aCycle) << " : ";
	
	std::list<mccore::ResId>::const_iterator it;
	for(it = aCycle.residues().begin(); it != aCycle.residues().end(); ++ it)
	{
		if(it != aCycle.residues().begin())
		{
			oss << "-";
		}
		oss << *it;
	}
	oss << " : " << getResiduesString(aModel, aCycle);
	
	return oss.str();
}

int main (int argc, char *argv[])
{
	std::set<annotate::Cycle> cycles;
	read_options (argc, argv);

	while (optind < argc)
	{
		Molecule *molecule;
		Molecule::iterator molIt;
		std::string filename = (std::string) argv[optind];
      
		molecule = loadFile (filename);
		if (0 != molecule)
		{
			unsigned int uiModel = 1;
			for (molIt = molecule->begin (); molecule->end () != molIt; ++molIt)
			{
				if (guiModelNumber != 0 && uiModel != guiModelNumber)
				{
					++ uiModel;
				}
				else
				{
					annotate::AnnotateModel &am = (annotate::AnnotateModel&) *molIt;
					am.name(getPdbFileName(filename));
					annotate::AnnotationInteractions annInteractions;
					annotate::AnnotationChains annChains;
					annotate::AnnotationStems annStems;
					annotate::AnnotationLinkers annLinkers;
					annotate::AnnotationLoops annLoops;
		  
		  			am.addAnnotation(annInteractions);
		  			am.addAnnotation(annChains);
					am.addAnnotation(annStems);
					am.addAnnotation(annLinkers);
					am.addAnnotation(annLoops);
					
					am.annotate (gucRelationMask);
					
					cycles = computeSecondaryStructureCycles(am);
					
					std::set<annotate::Cycle>::const_iterator itCycle;
					for( itCycle = cycles.begin(); itCycle != cycles.end(); ++ itCycle)
					{
						mccore::gOut(0) << getPdbFileName(filename) << " : ";
						mccore::gOut(0) << uiModel << " : ";
						mccore::gOut(0) << cycleString(am, *itCycle);
						mccore::gOut(0) << std::endl;
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
	
	return EXIT_SUCCESS;	
}
