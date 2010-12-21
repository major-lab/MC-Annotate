/*
 * PDB2DBConsoleApp.cc
 *
 *  Created on: Jul 21, 2010
 *      Author: blanchmf
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "PDB2DBConsoleApp.h"

#include <mccore/Messagestream.h>
#include <mccore/Version.h>

#include <cstdio>
#include <cstdlib>

static const char* shortopts = "Vbc:s:e:f:ghlr:vm:";

PDB2DotBracketParamsConsole::PDB2DotBracketParamsConsole(int argc, char * argv []) : PDB2DotBracketParams()
{
	read_options (argc, argv);
}

void PDB2DotBracketParamsConsole::version () const
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void PDB2DotBracketParamsConsole::usage () const
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-ghlsvV] [-f <model number>] [-c <max layer>] [-s <max layer>] <structure file> ..."
		<< std::endl;
}

void PDB2DotBracketParamsConsole::help () const
{
	mccore::gOut (0)
		<< "This program annotate structures (and more)." << std::endl
		<< "  -c max layer      specify the number of layers to display as combined notation." << std::endl
		<< "  -s max layer	    display the number of layers to display separately" << std::endl
		<< "  -f model number   model to annotate (default : all models) " << std::endl
		<< "  -g                Insert \'X\' in sequence and \'.\' in dot bracket notation where nucleotides are missing from the model." << std::endl
		<< "  -h                print this help" << std::endl
		<< "  -l                be more verbose (log)" << std::endl
		<< "  -v                be verbose" << std::endl
		<< "  -V                print the software version info" << std::endl;
}

void PDB2DotBracketParamsConsole::read_options (int argc, char* argv[])
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
		case 'f':
		{
			long int tmp;

			tmp = strtol (optarg, 0, 10);
			if (ERANGE == errno	|| EINVAL == errno || 0 > tmp)
			{
				mccore::gErr(0) << PACKAGE << ": invalid model value.";
				mccore::gErr(0) << std::endl;
				exit (EXIT_FAILURE);
			}
			muiModelNumber = tmp;
			mbOneModel = true;
			break;
		}
		case 'g':
			mbCompleteGaps = true;
			break;
		case 's':
			muiSplitLayers = strtol(optarg, 0, 10);
			break;
		case 'c':
			muiCombinedLayers = strtol(optarg, 0, 10);
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

	if (argc - optind < 1)
	{
		usage ();
		exit (EXIT_FAILURE);
	}else
	{
		while (optind < argc)
		{
			std::string strFileName = argv[optind];
			mFiles.push_back(strFileName);
			++optind;
		}
	}
}

int main(int argc, char* argv[])
{
	PDB2DotBracketParamsConsole params(argc, argv);
	PDB2DotBracket theApp(params, params.mFiles);
	return EXIT_SUCCESS;
}
