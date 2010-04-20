/*
 * molsepRNA.cc
 *
 *  Created on: Apr 07, 2010
 *      Author: blanchmf
 */

#include "molsepRNA.h"
#include <list>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

#include <mccore/Binstream.h>
#include <mccore/Messagestream.h>
#include <mccore/Molecule.h>
#include <mccore/Pdbstream.h>
#include <mccore/ResidueFactoryMethod.h>

#include "AnnotateModel.h"
#include "StringTable.h"
#include "StringUtil.h"

// Options pour la lecture de la ligne de commande
static const char* gszShortopts = "Vhcd:";

/**
 * Constructeur de la classe.
 */
molsepRNA::molsepRNA(int argc, char * argv [])
: mfAppVersion(1.0f),
mstrAppName("molsepRNA"),
mbSplitChains(false)
{
	read_options (argc, argv);
	while (optind < argc)
	{
		std::string strFile = (std::string) argv[optind];
		mccore::gOut(0) << "Loading " << strFile << std::endl;
		mccore::Molecule mol = loadMolecule(strFile, false); // TODO : Support binary
		unsigned int uiModel = 1;
		mccore::Molecule::iterator it;
		for(it = mol.begin();it != mol.end(); ++ it, ++uiModel)
		{
			it->keepRNA();
			if(0 < it->size())
			{
				if(mbSplitChains)
				{
					writeChains(strFile, uiModel, *it);
				}else
				{
					writeModel(strFile, uiModel, *it);
				}

			}
		}
		++optind;
	}
}

/**
 * Affiche l'information de version du logiciel
 */
void molsepRNA::version () const
{
	std::cout << mstrAppName << " "; // Nom du logiciel
	std::cout << mfAppVersion << " ";			// Version du logiciel
	std::cout << "(" << __DATE__ << ")";	// Date de la compilation
	std::cout << std::endl;
}

/**
 * Affiche l'information sur la ligne de commande correspondant au logiciel.
 */
void molsepRNA::usage () const
{
	std::cout << "usage: [-hVc] [-d <Output directory>] <structure file> ...";
	std::cout << std::endl;
}


/**
 * Affiche une description detaillee des differents parametres acceptes par
 * l'application.
 */
void molsepRNA::help () const
{
	std::cout
		<< "Separates a PDB file in individual files for each model." << std::endl
		<< "Files are saved in \"directory\" (defaults to current) with the model number appended to the input file name." << std::endl
		<< "  -c				Also separate chains in different files"
		<< "  -h                print this information" << std::endl
		<< "  -V                print the software version info" << std::endl
		<< "  -d                output directory" << std::endl;
}

/**
 * Effectue la lecture des options de la ligne de commande.  Les parametres sont
 * conserves dans des variables globales.
 */
void molsepRNA::read_options (int argc, char* argv[])
{
	int c;

	while ((c = getopt (argc, argv, gszShortopts)) != -1)
	{
		switch (c)
		{
		case 'V':
			version ();
			exit (EXIT_SUCCESS);
			break;
		case 'h':
			help();
			exit(EXIT_SUCCESS);
			break;
		case 'c':
			mbSplitChains = true;
			break;
		case 'd':
			mstrOutputPath = optarg;
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if(optind >= argc)
	{
		usage ();
		exit (EXIT_FAILURE);
	}

	if(0 == mstrOutputPath.size())
	{
		mstrOutputPath = ".";
	}
}

std::string molsepRNA::getFilePrefix(const std::string& astrFile) const
{
	std::string::size_type index;
	std::string filename = astrFile;
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

mccore::Molecule molsepRNA::loadMolecule(const string &filename, bool abBinary)
{
	mccore::ResidueFM rFM;
	annotate::AnnotateModelFM aFM (mccore::ResIdSet(), 0, &rFM);
	mccore::Molecule molecule(&aFM);

	if (abBinary)
	{
		mccore::izfBinstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			std::ostringstream oss;
			oss << PACKAGE << ": cannot open binary file '" << filename << "'.";
			throw mccore::FileNotFoundException(oss.str(), __FILE__, __LINE__);
		}
		in >> molecule;
		in.close ();
	}
	else
	{
		mccore::izfPdbstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			std::ostringstream oss;
			oss << PACKAGE << ": cannot open pdb file '" << filename << "'.";
			throw mccore::FileNotFoundException(oss.str(), __FILE__, __LINE__);
		}
		in >> molecule;
		in.close ();
	}
	return molecule;
}

void molsepRNA::writeModel(
	const std::string astrOriginalFile,
	unsigned int auiModel,
	const mccore::AbstractModel& aModel) const
{
	std::ostringstream oss;
	oss << mstrOutputPath << "/" << getFilePrefix(astrOriginalFile);
	oss << "_" << auiModel << ".pdb.gz";
	std::string strOutputFile = oss.str ();

	mccore::ozfPdbstream out;
	out.open(strOutputFile.c_str());
	if (out.fail ())
	{
		std::ostringstream ossErr;
		ossErr << PACKAGE << ": cannot open pdb file '" << strOutputFile << "'.";
		throw mccore::FileNotFoundException(ossErr.str(), __FILE__, __LINE__);
	}
	mccore::gOut(0) << "Writing " << strOutputFile << std::endl;
	out << aModel;
	out.close ();
}

void molsepRNA::writeChains(
	const std::string astrOriginalFile,
	unsigned int auiModel,
	const mccore::AbstractModel& aModel) const
{
	std::map<char, mccore::GraphModel>::iterator itModel;
	std::map<char, mccore::GraphModel> models;
	mccore::GraphModel::const_iterator it = aModel.begin();
	for(; it != aModel.end(); ++ it)
	{
		char cChain = it->getResId().getChainId();
		itModel = models.find(cChain);
		if(itModel == models.end())
		{
			std::pair<std::map<char, mccore::GraphModel>::iterator, bool> retInsert;
			retInsert = models.insert(std::pair<char, mccore::GraphModel>(cChain, mccore::GraphModel()));
			assert(retInsert.second);
			itModel = retInsert.first;
		}
		itModel->second.insert(*it);
	}

	for(itModel = models.begin(); itModel != models.end(); ++itModel)
	{
		std::ostringstream oss;
		oss << mstrOutputPath << "/" << getFilePrefix(astrOriginalFile);
		oss << "_" << auiModel;
		oss << "_" << itModel->first << ".pdb.gz";
		std::string strOutputFile = oss.str ();

		mccore::ozfPdbstream out;
		out.open(strOutputFile.c_str());
		if (out.fail ())
		{
			std::ostringstream ossErr;
			ossErr << PACKAGE << ": cannot open pdb file '" << strOutputFile << "'.";
			throw mccore::FileNotFoundException(ossErr.str(), __FILE__, __LINE__);
		}
		mccore::gOut(0) << "Writing " << strOutputFile << std::endl;
		out << itModel->second;
		out.close ();
	}
}

int main(int argc, char* argv[])
{
	molsepRNA theApp(argc, argv);
	return EXIT_SUCCESS;
}
