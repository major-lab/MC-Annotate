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

#include <cassert>
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
	mccore::ResIdSet residueSelection;
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
	const annotate::AnnotateModel& aModel)
{
	// Compute the profile
	std::vector<std::list<mccore::ResId> > profile;
	const annotate::AnnotationInteractions* pInteractions = 
		aModel.getAnnotation<annotate::AnnotationInteractions>();
	
	// Adds the structure
	std::list<mccore::ResId> strand;
	GraphModel::const_iterator itPrev = aModel.end();
	GraphModel::const_iterator it;
	for(it = aModel.begin(); it != aModel.end(); ++ it)
	{
		if(	it == aModel.begin() 
			|| pInteractions->areContiguous(
				itPrev->getResId(), 
				it->getResId()))
		{
			strand.push_back(it->getResId());
		}
		else
		{
			profile.push_back(strand);
			strand.clear();
			strand.push_back(it->getResId());
		}
		itPrev = it;
	}
	profile.push_back(strand);
	strand.clear();
	return profile;
}

std::list<unsigned int> getProfile(const annotate::AnnotateModel& aModel)
{
	// Compute the profile
	std::list<unsigned int> profile;
	std::vector<std::list<mccore::ResId> > resProfile = getResProfile(aModel);
	
	std::vector<std::list<mccore::ResId> >::const_iterator itStrand;
	for(itStrand = resProfile.begin(); 
		itStrand != resProfile.end(); 
		++ itStrand)
	{
		profile.push_back(itStrand->size());
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

/** 
 * @brief Creates a cycle from the annotate model, and using provided profile.
 * @details It is necessary to provide the profile as the same residues may 
 * form different cycles depending on what interactions constitute them.
 */
annotate::Cycle getCycle(
	const annotate::AnnotateModel& aModel, 
	const std::list<unsigned int>& aExpectedProfile)
{
	std::vector<std::vector<mccore::ResId> > strands;

	// Try and make a cycle from the given information
	GraphModel::const_iterator itResidue = aModel.begin();
	std::list<unsigned int>::const_iterator itStrand;
	for(itStrand = aExpectedProfile.begin(); 
		itStrand != aExpectedProfile.end(); 
		++ itStrand)
	{
		std::vector<mccore::ResId> strand;
		for(unsigned int iRes = 0; iRes < *itStrand; ++ iRes)
		{
			strand.push_back(itResidue->getResId());
			++ itResidue;
		}
		strands.push_back(strand);
		strand.clear();
	}

	return createCycleFromStrands(aModel, strands);	
}

annotate::Cycle createCycleFromStrands(
	const annotate::AnnotateModel& aModel, 
	std::vector<std::vector<mccore::ResId> > aStrands)
{
	// Try and make a cycle from the given information
	const annotate::AnnotationInteractions* pInteractions = 
		aModel.getAnnotation<annotate::AnnotationInteractions>();
	
	std::set<annotate::BaseInteraction> interactions;
	std::vector<std::vector<mccore::ResId> >::const_iterator itStrand;
	
	// Connects adjacency
	for(itStrand = aStrands.begin(); itStrand != aStrands.end(); ++ itStrand)
	{
		std::vector<mccore::ResId>::const_iterator itRes;
		std::vector<mccore::ResId>::const_iterator itPrevRes = itStrand->end();
		for(itRes = itStrand->begin(); itRes != itStrand->end(); ++ itRes)
		{
			if(itRes != itStrand->begin())
			{
				std::list<const annotate::BaseInteraction*> inters;
				inters = pInteractions->getInteractions(*itPrevRes, *itRes);
				
				const annotate::BaseLink* pLink = findFirstBaseLink(inters);
				
				assert(NULL != pLink); 
				annotate::BaseInteraction inter1(
					pLink->first,
					pLink->fResId, 
					pLink->second,
					pLink->rResId);
				interactions.insert(inter1);
			}
			itPrevRes = itRes;
		}
	}
	
	// Connects the strands
	if(1 == aStrands.size())
	{
		std::list<const annotate::BaseInteraction*> inters;
		inters = pInteractions->getInteractions(
			aStrands[0].front(), 
			aStrands[0].back());
			
		const annotate::BasePair* pPair = findFirstBasePair(inters);
		annotate::BaseInteraction inter1(
					pPair->first,
					pPair->fResId, 
					pPair->second,
					pPair->rResId);
		interactions.insert(inter1);
	}
	else if(2 == aStrands.size())
	{
		std::list<const annotate::BaseInteraction*> inters;
		inters = pInteractions->getInteractions(
			aStrands[0].front(), 
			aStrands[1].back());
			
		const annotate::BasePair* pPair = findFirstBasePair(inters);
		assert(NULL != pPair);
		annotate::BaseInteraction inter1(
					pPair->first,
					pPair->fResId, 
					pPair->second,
					pPair->rResId);
					
		inters = pInteractions->getInteractions(
			aStrands[0].back(), 
			aStrands[1].front());
			
		pPair = findFirstBasePair(inters);
		assert(NULL != pPair);
		annotate::BaseInteraction inter2(
			pPair->first, 	pPair->fResId, 
			pPair->second,	pPair->rResId);
		
		interactions.insert(inter1);
		interactions.insert(inter2);
		
	} else
	{
		mccore::gOut(0) << "1 or 2 strands cycle supported" << std::endl;
		assert(false);
	}
	return annotate::Cycle(aModel, interactions, gucRelationMask);	
}

const annotate::BaseLink* findFirstBaseLink(
	std::list<const annotate::BaseInteraction*>& aInteractions)
{
	const annotate::BaseLink* pFound = NULL;
	std::list<const annotate::BaseInteraction*>::const_iterator it;
	for(it = aInteractions.begin(); 
		it != aInteractions.end() && NULL == pFound; 
		++ it)
	{
		pFound = dynamic_cast<const annotate::BaseLink*>(*it);
	}
	return pFound;
}

const annotate::BasePair* findFirstBasePair(
	std::list<const annotate::BaseInteraction*>& aInteractions)
{
	const annotate::BasePair* pFound = NULL;
	std::list<const annotate::BaseInteraction*>::const_iterator it;
	for(it = aInteractions.begin(); 
		it != aInteractions.end() && NULL == pFound; 
		++ it)
	{
		pFound = dynamic_cast<const annotate::BasePair*>(*it);
	}
	return pFound;
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
				
				std::list<unsigned int> profile = getProfile(am);
				mccore::gOut(0) << std::endl;
				if(2 < profile.size())
				{
					mccore::gOut(0) << "Maximum of 2 strand supported" << std::endl;
					exit(0);
				}
				
				mccore::gOut(0) << getPdbFileName(filename) << " : ";
				mccore::gOut(0) << getModelIndex(filename) << " : ";
				std::list<unsigned int>::const_iterator itProfile;
				for(itProfile = profile.begin(); 
					itProfile != profile.end(); 
					++ itProfile)
				{
					mccore::gOut(0) << *itProfile;
				}
				mccore::gOut(0) << " : ";
				
				std::string strFileProfile = getModelProfile(filename);
				
				mccore::gOut(0) << strFileProfile << "\t: ";
				std::list<unsigned int> fileProfile = 
					getExpectedProfile(strFileProfile);
				
				annotate::Cycle cycle = getCycle(am, fileProfile);
				mccore::gOut(0) << " cp ";
				std::vector<unsigned int>::const_iterator itCProf;
				for(itCProf = cycle.profile().begin(); 
					itCProf != cycle.profile().end(); 
					++ itCProf)
				{
					mccore::gOut(0) << *itCProf;
				}
				mccore::gOut(0) << " : ";
				
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
				mccore::gOut(0) << std::endl;
			}
			delete molecule;
		}
		++optind;
	}
	return EXIT_SUCCESS;	
}
