//                              -*- Mode: C++ -*-
// mc3dssc.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Mar 18 14:04:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mc3dssc.h"

#include "CycleInfo.h"

#include "CycleInfoFile.h"
#include "InteractionInfoFile.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

static const char* shortopts = "Vhlvc:p:";

MC3DSecondaryStructureCycles::MC3DSecondaryStructureCycles(int argc, char * argv [])
{
	// Read the command line options
	readOptions(argc, argv);

	// Read the interactions
	annotate::InteractionInfoFile infileInteractions;
	infileInteractions.read(mstrInteractionsFile.c_str());
	std::set<annotate::InteractionInfo>::const_iterator itInter;
	for(itInter = infileInteractions.interactions().begin();
		itInter != infileInteractions.interactions().end();
		++ itInter)
	{
		annotate::ModelInfo model = itInter->getModelInfo();
		mInteractions[model].insert(*itInter);
		mModels.insert(model);
	}

	// Read the cycles
	annotate::CycleInfoFile infileCycles;
	infileCycles.read(mstrCyclesFile.c_str());
	std::set<annotate::CycleInfo>::const_iterator itCycle;
	for(itCycle = infileCycles.cycles().begin();
		itCycle != infileCycles.cycles().end();
		++ itCycle)
	{
		annotate::ModelInfo model = itCycle->getModelInfo();
		mCycles[model].insert(*itCycle);
		mModels.insert(model);
	}

	identifyCycles();
}

void MC3DSecondaryStructureCycles::version () const
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void MC3DSecondaryStructureCycles::usage () const
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-hlvV] -c <cycles file> -p <distant pairs file>"
		<< std::endl;
}

void MC3DSecondaryStructureCycles::help () const
{
	mccore::gOut (0)
		<< "This program identifies cycles belonging to the secondary structure." << std::endl
		<< "  -c	file containing the identified cycles" << std::endl
		<< "  -p	file containing the non-adjacent interacting pairs" << std::endl
		<< "  -h	print this help" << std::endl
		<< "  -l	be more verbose (log)" << std::endl
		<< "  -v	be verbose" << std::endl
		<< "  -V	print the software version info" << std::endl;
}


void MC3DSecondaryStructureCycles::readOptions (int argc, char* argv[])
{
	int c;

	while ((c = getopt (argc, argv, shortopts)) != EOF)
	{
		switch (c)
		{
		case 'c':
		{
			mstrCyclesFile = optarg;
			break;
		}
		case 'p':
		{
			mstrInteractionsFile = optarg;
			break;
		}
		case 'V':
			version ();
			exit (EXIT_SUCCESS);
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

	if(0 == mstrInteractionsFile.size() || 0 == mstrCyclesFile.size())
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

void MC3DSecondaryStructureCycles::identifyCycles()
{
	std::set<annotate::ModelInfo>::const_iterator itModel;
	for(itModel = mModels.begin(); itModel != mModels.end(); ++ itModel)
	{
		identifyModelCycles(*itModel);
	}
}

void MC3DSecondaryStructureCycles::identifyModelCycles(const annotate::ModelInfo& aModel)
{
	std::map<annotate::ModelInfo, interaction_set >::iterator itInterSet = mInteractions.find(aModel);
	std::map<annotate::ModelInfo, cycle_set >::iterator itCycleSet = mCycles.find(aModel);
	if(itInterSet != mInteractions.end() && itCycleSet != mCycles.end())
	{
		cycle_set cycleSet = itCycleSet->second;
		cycle_set::const_iterator itCycle;

		// Remove cycles with parallel profile or happening between different
		// molecules
		cycle_set toRemove;
		for(itCycle = cycleSet.begin(); itCycle != cycleSet.end(); ++ itCycle)
		{
			if( itCycle->getProfile().type() == annotate::Cycle::e2STRANDS_PARALLEL
							    || !singleChain(*itCycle))
			{
				toRemove.insert(*itCycle);
			}
		}
		cycleSet = annotate::SetDifference(cycleSet, toRemove);
		toRemove.clear();

		// Remove cycles containing tertiary interaction
		interaction_set::const_iterator itInter;
		for(itInter = itInterSet->second.begin();
			itInter != itInterSet->second.end();
			++ itInter)
		{
			for(itCycle = cycleSet.begin(); itCycle != cycleSet.end(); ++ itCycle)
			{
				if(itCycle->contains(*itInter))
				{
					toRemove.insert(*itCycle);
				}
			}
		}
		cycleSet = annotate::SetDifference(cycleSet, toRemove);
		outputCycles(cycleSet);
	}
}

bool MC3DSecondaryStructureCycles::singleChain(const annotate::CycleInfo& aCycle) const
{
	bool bSingle = true;
	char chChainId = ' ';
	std::vector<mccore::ResId> ids = aCycle.getResIds();
	std::vector<mccore::ResId>::const_iterator it;
	for(it = ids.begin(); it != ids.end() && bSingle; ++ it)
	{
		try
		{
		if(it == ids.begin())
		{
			chChainId = it->getChainId();
		}else if(it->getChainId() != chChainId)
		{
			bSingle = false;
		}
		} catch(mccore::FatalLibException except)
		{
			std::cerr << except << " : " << *it << std::endl;
		}
	}
	return bSingle;
}

void MC3DSecondaryStructureCycles::outputCycles(
	std::set<annotate::CycleInfo>& aCycles) const
{
	std::set<annotate::CycleInfo>::const_iterator it = aCycles.begin();
	for(it = aCycles.begin(); it != aCycles.end(); ++ it)
	{
		std::cout << it->getPDBFile() << " : " << it->getModel();
		std::cout << " : " << it->getProfile().toString();
		std::cout << " : " << it->resIdsString("-");
		std::cout << " : " << it->residuesString("-") << std::endl;
	}
}

int main (int argc, char *argv[])
{
	MC3DSecondaryStructureCycles theApp(argc, argv);

	return EXIT_SUCCESS;
}
