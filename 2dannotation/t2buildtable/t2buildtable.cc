//                              -*- Mode: C++ -*-
// t1buildtable.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Mar 17 11:14:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t2buildtable.h"
#include "t2buildtableh1.h"
#include "t2tablebuilderh2.h"
#include "t2tablebuilderh3.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

static const char* shortopts = "Vhlvi:g:s:";

T2BuildTable::T2BuildTable(int argc, char * argv [])
{
	// Read the command line options
	readOptions(argc, argv);

	// Create the build table
	switch(miHypothesis)
	{
	case 1:
		{
			T2BuildTableH1 tableBuilder(mstrPDBGroupsFile, mstrInteractionsFile);
			tableBuilder.computeTable();
		}
		break;
	case 2:
		{
			T2TableBuilderH2 tableBuilder(mstrPDBGroupsFile, mstrInteractionsFile);
			tableBuilder.computeTable();
		}
		break;
	case 3:
		{
			T2TableBuilderH3 tableBuilder(mstrPDBGroupsFile, mstrInteractionsFile);
			tableBuilder.computeTable();
		}
		break;
	default:
		usage ();
		exit (EXIT_FAILURE);
	}
}

void T2BuildTable::version () const
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void T2BuildTable::usage () const
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-hlvV] -s <hypothese> -g <PDB Groups File> -i <interactions file>"
		<< std::endl;
}

void T2BuildTable::help () const
{
	mccore::gOut (0)
		<< "This program computes a table of frequency of interacting cycle pairs." << std::endl
		<< "  -s	identifier of the hypothesis to test" << std::endl
		<< "  -g    file containing the PDB File group for normalization" << std::endl
		<< "  -c	file containing the identified 3D interactions" << std::endl
		<< "  -h	print this help" << std::endl
		<< "  -l	be more verbose (log)" << std::endl
		<< "  -v	be verbose" << std::endl
		<< "  -V	print the software version info" << std::endl;
}

void T2BuildTable::readOptions (int argc, char* argv[])
{
	int c;

	while ((c = getopt (argc, argv, shortopts)) != EOF)
	{
		switch (c)
		{
		case 'g':
		{
			mstrPDBGroupsFile = optarg;
			break;
		}
		case 'i':
		{
			mstrInteractionsFile = optarg;
			break;
		}
		case 's':
		{
			miHypothesis = atoi(optarg);
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

	if(0 == mstrInteractionsFile.size() || 0 == mstrPDBGroupsFile.size())
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

int main (int argc, char *argv[])
{
	T2BuildTable theApp(argc, argv);

	return EXIT_SUCCESS;
}
