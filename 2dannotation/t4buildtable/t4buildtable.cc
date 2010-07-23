//                              -*- Mode: C++ -*-
// t4buildtable.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Sat May 2 12:09:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t4buildtable.h"
#include "t4tablebuilderh1.h"

#include <mccore/Messagestream.h>
#include <mccore/Version.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>

static const char* shortopts = "Vhlvi:g:s:";

T4BuildTable::T4BuildTable(int argc, char * argv [])
{
	// Read the command line options
	readOptions(argc, argv);

	// Create the build table
	switch(miHypothesis)
	{
	case 1:
		{
			T4BuildTableH1 tableBuilder(mstrPDBGroupsFile, mstrInteractionsFile);
			tableBuilder.computeTable();
		}
		break;
	default:
		usage ();
		exit (EXIT_FAILURE);
	}
}

void T4BuildTable::version () const
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void T4BuildTable::usage () const
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-hlvV] -s <hypothese> -g <PDB Groups File> -i <interactions file>"
		<< std::endl;
}

void T4BuildTable::help () const
{
	mccore::gOut (0)
		<< "This program computes a table of frequency of tertiary interaction type given sequence." << std::endl
		<< "  -s	identifier of the hypothesis to test" << std::endl
		<< "  -g    file containing the PDB File group for normalization" << std::endl
		<< "  -c	file containing the identified 3D interactions" << std::endl
		<< "  -h	print this help" << std::endl
		<< "  -l	be more verbose (log)" << std::endl
		<< "  -v	be verbose" << std::endl
		<< "  -V	print the software version info" << std::endl;
}

void T4BuildTable::readOptions (int argc, char* argv[])
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
	T4BuildTable theApp(argc, argv);

	return EXIT_SUCCESS;
}
