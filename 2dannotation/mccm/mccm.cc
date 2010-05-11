//                              -*- Mode: C++ -*-
// mccm.cc
// Copyright © 2001-09 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Jul 29 10:35:00 2009

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mccm.h"

// libmcannotate.a
#include "Cycle.h"
#include "CycleProfile.h"
#include "StringTable.h"
#include "StringUtil.h"

#include "CycleInfo.h"
#include "CycleInfoFile.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>

bool binary = false;
bool oneModel = false;
std::string gstrPrimaryCyclesFile;
std::string gstrSecondaryCyclesFile;
const char* shortopts = "Vbhlv";

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
		<< " [-bhlvV] <first cycles file> <second cycles file>"
		<< std::endl;
}

void help ()
{
	mccore::gOut (0)
		<< "Merge two files containing cycles into one, removing duplicates." << std::endl
		<< "  -h                print this help" << std::endl
		<< "  -c				First cycles file" << std::endl
		<< "  -s				Second cycles file" << std::endl
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
		case 'c':
		{
			gstrPrimaryCyclesFile = optarg;
			break;
		}
		case 's':
		{
			gstrSecondaryCyclesFile = optarg;
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
		case 'v':
			mccore::gOut.setVerboseLevel (mccore::gOut.getVerboseLevel () + 1);
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if (argc - optind != 2)
	{
		usage ();
		exit (EXIT_FAILURE);
	}
	else
	{
		gstrPrimaryCyclesFile = argv[optind ++];
		gstrSecondaryCyclesFile = argv[optind];
	}
}

std::set<annotate::CycleInfo> readCyclesFile(const std::string& astrFile)
{
	annotate::CycleInfoFile infile;
	infile.read(astrFile);

	return infile.cycles();
}

std::string getModelString(const annotate::CycleInfo& aCycle)
{
	std::ostringstream oss;
	oss << aCycle.getModel();
	return oss.str();
}

std::string getResIdsString(const annotate::CycleInfo& aCycle)
{
	std::ostringstream oss;
	std::vector<mccore::ResId> resIds = aCycle.getResIds();
	std::vector<mccore::ResId>::const_iterator it;
	for(it = resIds.begin(); it != resIds.end(); ++ it)
	{
		if(it != resIds.begin())
		{
			oss << "-";
		}
		oss << *it;
	}
	return oss.str();
}

std::string getSequenceString(const annotate::CycleInfo& aCycle)
{
	std::ostringstream oss;
	std::vector<std::string> res = aCycle.getSequence();
	std::vector<std::string>::const_iterator it;
	for(it = res.begin(); it != res.end(); ++ it)
	{
		if(it != res.begin())
		{
			oss << "-";
		}
		oss << *it;
	}
	return oss.str();
}

int main (int argc, char *argv[])
{
	annotate::StringTable stringTable(6);

	read_options (argc, argv);

	std::set<annotate::CycleInfo> cycles = readCyclesFile(gstrPrimaryCyclesFile);
	std::set<annotate::CycleInfo> cycles2 = readCyclesFile(gstrSecondaryCyclesFile);
	cycles.insert(cycles2.begin(), cycles2.end());

	std::set<annotate::CycleInfo>::const_iterator it;
	for(it = cycles.begin(); it != cycles.end(); ++ it)
	{
		std::vector<string>& tableRow = stringTable.addRow();
		tableRow[0] = it->getPDBFile();
		tableRow[1] = getModelString(*it);
		tableRow[2] = it->getProfile().toString();
		tableRow[3] = it->getProfile().toString();
		tableRow[4] = getResIdsString(*it);
		tableRow[5] = getSequenceString(*it);
	}

	mccore::gOut(0) << stringTable.toString(" : ");

	return EXIT_SUCCESS;
}
