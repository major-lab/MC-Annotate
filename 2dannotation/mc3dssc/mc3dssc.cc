//                              -*- Mode: C++ -*-
// mc3dicp.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Mar 22 10:22:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mc3dicp.h"

#include "CycleInfo.h"

#include "CycleInfoFile.h"
#include "InteractionInfoFile.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

static const char* shortopts = "Vhlvc:p:";

MC3DInteractingCyclePairs::MC3DInteractingCyclePairs(int argc, char * argv [])
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
		annotate::CycleInfo cycle = *itCycle;
		if(2 == cycle.getNbStrands() && cycle.getStrandResIds()[0].size() > cycle.getStrandResIds()[1].size())
		{
			std::cout << "Flipping strand" << std::endl;
			assert(false);
			cycle = annotate::CycleInfo::flipStrand(cycle);
		}
		mCycles[model].insert(cycle);
		mModels.insert(model);
	}

	identifyPairs();
}

void MC3DInteractingCyclePairs::version () const
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void MC3DInteractingCyclePairs::usage () const
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-hlvV] -c <cycles file> -p <distant pairs file>"
		<< std::endl;
}

void MC3DInteractingCyclePairs::help () const
{
	mccore::gOut (0)
		<< "This program identifies pairs of cycles interacting together." << std::endl
		<< "  -c	file containing the identified cycles" << std::endl
		<< "  -p	file containing the non-adjacent interacting pairs" << std::endl
		<< "  -h                print this help" << std::endl
		<< "  -l                be more verbose (log)" << std::endl
		<< "  -v                be verbose" << std::endl
		<< "  -V                print the software version info" << std::endl;
}


void MC3DInteractingCyclePairs::readOptions (int argc, char* argv[])
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

void MC3DInteractingCyclePairs::identifyPairs()
{
	std::set<annotate::ModelInfo>::const_iterator itModel;
	for(itModel = mModels.begin(); itModel != mModels.end(); ++ itModel)
	{
		identifyModelPairs(*itModel);
	}
}

void MC3DInteractingCyclePairs::identifyModelPairs(const annotate::ModelInfo& aModel)
{
	std::map<annotate::ModelInfo, interaction_set >::const_iterator itInterSet = mInteractions.find(aModel);
	std::map<annotate::ModelInfo, cycle_set >::const_iterator itCycleSet = mCycles.find(aModel);
	if(itInterSet != mInteractions.end() && itCycleSet != mCycles.end())
	{
		interaction_set::const_iterator itInter;
		for(itInter = itInterSet->second.begin();
			itInter != itInterSet->second.end();
			++ itInter)
		{
			cycle_set::const_iterator itCycle;
			std::set<annotate::CycleInfo> left;
			std::set<annotate::CycleInfo> right;
			for(itCycle = itCycleSet->second.begin();
				itCycle != itCycleSet->second.end();
				++ itCycle)
			{
				if(itCycle->contains(itInter->getRes1()) && !itCycle->contains(*itInter))
				{
					left.insert(*itCycle);
				}
				if(itCycle->contains(itInter->getRes2()) && !itCycle->contains(*itInter))
				{
					right.insert(*itCycle);
				}
			}
			outputPairs(left, right);
		}
	}
}

void MC3DInteractingCyclePairs::outputPairs(
	std::set<annotate::CycleInfo>& aLeft,
	std::set<annotate::CycleInfo>& aRight)
{
	std::set<std::pair<annotate::CycleInfo, annotate::CycleInfo> > pairs;
	std::set<annotate::CycleInfo>::const_iterator it1 = aLeft.begin();
	std::set<annotate::CycleInfo>::const_iterator it2 = aRight.begin();
	for(it1 = aLeft.begin(); it1 != aLeft.end(); ++ it1)
	{
		for(it2 = aRight.begin(); it2 != aRight.end(); ++ it2)
		{
			std::pair<annotate::CycleInfo, annotate::CycleInfo> interactingPair;
			if(*it1 < *it2)
			{
				interactingPair.first = *it1;
				interactingPair.second = *it2;
			}else
			{
				interactingPair.first = *it2;
				interactingPair.second = *it1;
			}
			pairs.insert(interactingPair);
		}
	}
	std::set<std::pair<annotate::CycleInfo, annotate::CycleInfo> >::const_iterator itPair;
	for(itPair = pairs.begin(); itPair != pairs.end(); ++ itPair)
	{
		// Model
		std::cout << itPair->first.getPDBFile() << " : " << itPair->first.getModel();

		// First Cycle
		std::cout << " ; " << itPair->first.getProfile().toString();
		std::cout << " : " << itPair->first.resIdsString("-");
		std::cout << " : " << itPair->first.residuesString("-");

		// Second cycle
		std::cout << " ; " << itPair->second.getProfile().toString();
		std::cout << " : " << itPair->second.resIdsString("-");
		std::cout << " : " << itPair->second.residuesString("-");

		std::cout << std::endl;
	}
}

int main (int argc, char *argv[])
{
	MC3DInteractingCyclePairs theApp(argc, argv);

	return EXIT_SUCCESS;
}
