//                              -*- Mode: C++ -*-
// db2ncm.cc
// Copyright © 2001-10 Laboratoire d'ingénierie des ARN.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Mon May 17 10:53:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "db2ncm.h"

#include "CycleInfo.h"

#include "CycleInfoFile.h"
#include "InteractionInfoFile.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <sstream>

static const char* shortopts = "Vhlvi:s:d:p:";

DotBracketToNCM::DotBracketToNCM(int argc, char * argv [])
{
	// Read the command line options
	readOptions(argc, argv);

	identifyPairs();

	std::vector<std::pair<unsigned int, unsigned int> > pairs = identifyPairs();

	mNCMs = identifyNCMs(1, pairs);
	std::set<annotate::CycleInfo>::const_iterator it;
	for(it = mNCMs.begin(); it != mNCMs.end(); ++ it)
	{
		std::cout << it->getPDBFile();
		std::cout << " : " << it->getFileProfile().toString();
		std::cout << " : " << it->residuesString("-");
		std::cout << " : " << it->resIdsString("-");
		std::cout.precision(10);
		std::cout << " : " << mfProbability;
		std::cout << std::endl;
	}
}

void DotBracketToNCM::version () const
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void DotBracketToNCM::usage () const
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-hlvV] -i <identifier> -c <cycles file> -p <distant pairs file>"
		<< std::endl;
}

void DotBracketToNCM::help () const
{
	mccore::gOut (0)
		<< "This read an identifier, a sequence, a list of resids and a " << std::endl
		<< "dot bracket notation and outputs a list of NCMs." << std::endl
		<< "  -i    Identifier" << std::endl
		<< "  -s	Sequence" << std::endl
		<< "  -r    Residue IDs" << std::endl
		<< "  -d	Dot Brackets" << std::endl
		<< "  -p    Probability associated" << std::endl
		<< "  -h    print this help" << std::endl
		<< "  -l    be more verbose (log)" << std::endl
		<< "  -v    be verbose" << std::endl
		<< "  -V    print the software version info" << std::endl;
}


void DotBracketToNCM::readOptions (int argc, char* argv[])
{
	int c;

	while ((c = getopt (argc, argv, shortopts)) != EOF)
	{
		switch (c)
		{
		case 'i':
		{
			mstrIdentifier = optarg;
			break;
		}
		case 's':
		{
			mstrSequence = optarg;
			break;
		}
		case 'd':
		{
			mstrDotBracket = optarg;
			break;
		}
		case 'p':
		{
			mfProbability = atof(optarg);
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

	if(0 == mstrSequence.size() || 0 == mstrDotBracket.size())
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

std::vector<DotBracketToNCM::index_pair> DotBracketToNCM::identifyPairs() const
{
	std::pair<unsigned int, unsigned int> def(0,0);
	std::vector<std::pair<unsigned int, unsigned int> > pairs;
	pairs.resize(mstrDotBracket.size(), def);
	std::list<unsigned int> opened;
	for(unsigned int i = 0; i < mstrDotBracket.size(); ++ i)
	{
		char c = mstrDotBracket[i];
		if(c == '(')
		{
			opened.push_back(i);
		}
		else if(c == ')')
		{
			unsigned int iOpen = opened.back();
			opened.pop_back();
			index_pair entry(iOpen, i);
			pairs[iOpen] = entry;
			pairs[i] = entry;
		}
	}
	return pairs;
}

std::set<annotate::CycleInfo> DotBracketToNCM::identifyNCMs(
	unsigned int auiDBId,
	const std::vector<DotBracketToNCM::index_pair>& aPairs) const
{
	std::set<annotate::CycleInfo> cycles;
	for(unsigned int i = 0; i < aPairs.size(); ++ i)
	{
		if(aPairs[i].first == i)
		{
			// Opening bracket
			std::vector<mccore::ResId> strand1;
			bool bFound = false;
			mccore::ResId res1;
			res1.setResNo(i + 1);
			strand1.push_back(res1);
			for(unsigned int j = i + 1; j < aPairs.size() && !bFound; ++ j)
			{
				if(mstrSequence[j] == 'X' || mstrSequence[j] == 'x')
				{
					bFound = true;
				}
				else if(aPairs[j].first == aPairs[j].second)
				{
					// Unpaired
					res1.setResNo(j + 1);
					strand1.push_back(res1);
				}else if(i == aPairs[j].first)
				{
					// Closing bracket, we've found a loop
					mccore::ResId res1;
					res1.setResNo(j + 1);
					strand1.push_back(res1);

					annotate::CycleInfo cycle = getCycleInfo(auiDBId, cycles.size() + 1, strand1);
					cycles.insert(cycle);
					bFound = true;
				}else if(j < aPairs[j].second)
				{
					// Another opening bracket, we have a 2 stranded NCM
					res1.setResNo(j + 1);
					strand1.push_back(res1);

					std::vector<mccore::ResId> strand2;
					mccore::ResId res2;
					res2.setResNo(aPairs[j].second + 1);
					strand2.push_back(res2);
					for(unsigned int k = aPairs[j].second + 1; k < aPairs.size() && !bFound; ++k)
					{
						if(aPairs[k].first == aPairs[k].second)
						{
							if('X' == mstrSequence[k] || 'x' == mstrSequence[k])
							{
								// Break out of this
								break;
							}else
							{
								// Unpaired
								res2.setResNo(k + 1);
								strand2.push_back(res2);
							}
						}else if(k == aPairs[i].second)
						{
							// We've closed the circle
							res2.setResNo(k + 1);
							strand2.push_back(res2);
							bFound = true;
						}else
						{
							assert(false);
						}
					}
					if(bFound)
					{
						annotate::CycleInfo cycle = getCycleInfo(auiDBId, cycles.size() + 1, strand1, strand2);
						cycles.insert(cycle);
					}
					else
					{
						bFound = true;
					}
				}
				else
				{
					assert(false);
				}
			}
		}
	}
	return cycles;
}

annotate::CycleInfo DotBracketToNCM::getCycleInfo(
	unsigned int auiDBId,
	unsigned int auiCycleId,
	const std::vector<mccore::ResId>& aStrand1,
	const std::vector<mccore::ResId>& aStrand2) const
{
	std::ostringstream oss;
	oss << aStrand1.size() << "_" << aStrand2.size();
	annotate::CycleProfile fileProfile(oss.str());
	annotate::CycleProfile profile = fileProfile;
	std::vector<mccore::ResId> residues;

	if(aStrand2.size() < aStrand1.size())
	{
		profile = annotate::CycleProfile::Rotate(fileProfile);
	}
	residues.insert(residues.end(), aStrand1.begin(), aStrand1.end());
	residues.insert(residues.end(), aStrand2.begin(), aStrand2.end());

	std::vector<std::string> sequence;
	sequence.resize(residues.size());
	for(unsigned int i = 0; i < residues.size(); ++ i)
	{
		sequence[i] = mstrSequence[residues[i].getResNo() - 1];
	}

	return annotate::CycleInfo(
		cycleName(auiDBId, auiCycleId),
		1,
		fileProfile,
		profile,
		residues,
		sequence);
}

annotate::CycleInfo DotBracketToNCM::getCycleInfo(
	unsigned int auiDBId,
	unsigned int auiCycleId,
	const std::vector<mccore::ResId>& aStrand1) const
{
	std::ostringstream oss;
	oss << aStrand1.size();
	annotate::CycleProfile profile(oss.str());
	std::vector<std::string> sequence;
	sequence.resize(aStrand1.size());
	for(unsigned int i = 0; i < aStrand1.size(); ++ i)
	{
		sequence[i] = mstrSequence[aStrand1[i].getResNo() - 1];
	}

	return annotate::CycleInfo(
		cycleName(auiDBId, auiCycleId),
		1,
		profile,
		profile,
		aStrand1,
		sequence);
}

std::string DotBracketToNCM::cycleName(
	unsigned int auiDBId,
	unsigned int auiCycleId) const
{
	std::ostringstream oss;
	oss << mstrIdentifier;
	oss << "." << setw(2) << setfill('0') << auiDBId;
	oss << "." << setw(2) << setfill('0') << auiCycleId;
	return oss.str();
}

int main (int argc, char *argv[])
{
	DotBracketToNCM theApp(argc, argv);

	return EXIT_SUCCESS;
}
