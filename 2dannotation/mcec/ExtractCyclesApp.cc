/*
 * ExtractCyclesApp.cc
 *
 *  Created on: Jan 28, 2010
 *      Author: blanchmf
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "ExtractCyclesApp.h"

#include "CycleInfoFile.h"

// libmcannotate
#include "AnnotateModel.h"

// libmccore
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/ResidueFactoryMethod.h"

#include <cstdio>
#include <fstream>
#include <sstream>

const char* shortopts = "Vhd:i:";

ExtractCyclesApp::ExtractCyclesApp(int argc, char * argv [])
{
	readOptions(argc, argv);

	// Read the cycle files
	readCyclesFile();
}

ExtractCyclesApp::~ExtractCyclesApp()
{
	mCyclesModels.clear();
	mCycles.clear();
}

void ExtractCyclesApp::version () const
{
	std::cout
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")"
		<< std::endl;
}

void ExtractCyclesApp::usage () const
{
	std::cout
	 	<< "usage: " << PACKAGE
		<< " [-hV] [-d <output directory>] -i <pdb directory> <cycles file> ..."
		<< std::endl;
}

void ExtractCyclesApp::help () const
{
	std::cout
		<< "This program read cycle structures and return the corresponding residue ids."
		<< std::endl
		<< "  -i	directory where the PDB files are stored." << std::endl
		<< "  -d	specify an output directory" << std::endl
		<< "  -h	print this help" << std::endl
		<< "  -V	print the software version info" << std::endl;
}

/**
 * readOptions
 * @brief Read the options from the command line.
 */
void ExtractCyclesApp::readOptions (int argc, char* argv[])
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
		case 'd':
		{
			mstrOutputDirectory = optarg;
			break;
		}
		case 'i':
		{
			mstrInputDirectory = optarg;
			break;
		}
		case 'h':
			usage ();
			help ();
			exit (EXIT_SUCCESS);
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	while (optind < argc)
	{
		mCyclesFiles.push_back((std::string) argv[optind]);
		++ optind;
	}

	if(mstrInputDirectory.empty())
	{
		usage ();
		exit (EXIT_FAILURE);
	} else if(0 == mCyclesFiles.size())
	{
		usage ();
		exit (EXIT_FAILURE);
	}

}

/**
 * readCyclesFile
 * @brief Read the cycles from a file.
 */
void ExtractCyclesApp::readCyclesFile()
{
	std::list<std::string>::const_iterator it;
	for(it = mCyclesFiles.begin(); it != mCyclesFiles.end(); ++ it)
	{
		annotate::CycleInfoFile infile;
		infile.read(it->c_str());
		std::set<annotate::CycleInfo> cycles = infile.cycles();
		mCycles.insert(cycles.begin(), cycles.end());
	}
}

/**
 * getModels
 * @brief Get a list of models from the cycles
 * @return List of models from which the cycles are originating
 */
std::list<annotate::ModelInfo> ExtractCyclesApp::getModels() const
{
	std::list<annotate::ModelInfo> models;
	std::set<annotate::CycleInfo>::const_iterator it;
	const annotate::ModelInfo* pInfo = NULL;
	for(it = mCycles.begin(); it != mCycles.end(); ++ it)
	{
		if(it != mCycles.begin() && (*pInfo != it->getModelInfo()))
		{
			models.push_back(*pInfo);
		}
		pInfo = &(it->getModelInfo());
	}
	if(NULL != pInfo)
	{
		models.push_back(*pInfo);
	}
	return models;
}

mccore::Molecule* ExtractCyclesApp::loadMoleculeFile (const std::string &filename) const
{
	mccore::Molecule *pMolecule = NULL;
	mccore::ResidueFM rFM;
	annotate::AnnotateModelFM aFM (ResIdSet(), 0/*environment*/, &rFM);

	mccore::izfPdbstream in;

	in.open (filename.c_str ());
	if (in.fail ())
	{
		mccore::gErr (0) << PACKAGE << ": cannot open pdb file '" << filename << "'." << std::endl;
		return 0;
	}else
	{
		pMolecule = new Molecule (&aFM);
		in >> *pMolecule;
		in.close ();
	}
	return pMolecule;
}

mccore::Molecule::const_iterator ExtractCyclesApp::getModel(
	mccore::Molecule& aMolecule,
	unsigned int auiModel) const
{
	unsigned int uiModel = 1;
	mccore::Molecule::const_iterator molIt;
	for (molIt = aMolecule.begin(); molIt != aMolecule.end (); ++molIt)
	{

		if (uiModel != auiModel)
		{
			++uiModel;
		}
		else
		{
			break;
		}
	}
	return molIt;
}

std::map<annotate::CycleInfo, mccore::GraphModel> ExtractCyclesApp::extract() const
{
	std::map<annotate::CycleInfo, mccore::GraphModel> results;
	std::list<annotate::ModelInfo> models = getModels();
	std::list<annotate::ModelInfo>::const_iterator itModel;
	for(itModel = models.begin(); itModel != models.end(); ++ itModel)
	{
		std::ostringstream oss;
		mccore::Molecule* pMolecule = NULL;
		oss << mstrInputDirectory << "/" << itModel->getPDBFile() << ".pdb.gz";
		pMolecule = loadMoleculeFile(oss.str());
		mccore::Molecule::const_iterator itMol = getModel(*pMolecule, itModel->getModel());
		GetModelRangeFunctor<annotate::CycleInfo> rangeFunctor;
		GetModelRangeFunctor<annotate::CycleInfo>::const_range modelRange;
		modelRange = rangeFunctor(mCycles, *itModel);
		std::map<annotate::CycleInfo, mccore::GraphModel> cyclesModels;
		mccore::GraphModel test = *itMol;
		cyclesModels = getCycleModels(*itMol, modelRange);
		results.insert(cyclesModels.begin(), cyclesModels.end());
		delete pMolecule;
	}
	return results;
}

std::map<annotate::CycleInfo, mccore::GraphModel> ExtractCyclesApp::getCycleModels(
	const mccore::GraphModel& aModel,
	GetModelRangeFunctor<annotate::CycleInfo>::const_range& aRange) const
{
	std::map<annotate::CycleInfo, mccore::GraphModel> models;
	std::set<annotate::CycleInfo>::const_iterator it;
	for(it = aRange.first; it != aRange.second; ++ it)
	{
		mccore::GraphModel model = getCycleModel(aModel, *it);
		models.insert(std::pair<annotate::CycleInfo, mccore::GraphModel>(*it, model));
	}
	return models;
}

mccore::GraphModel ExtractCyclesApp::getCycleModel(
	const mccore::GraphModel& aModel,
	const annotate::CycleInfo& aCycle) const
{
	mccore::GraphModel model;
	std::vector<std::string> resids = aCycle.getResIds();
	std::vector<std::string>::const_iterator it;
	for(it = resids.begin(); it != resids.end(); ++ it)
	{
		mccore::GraphModel::const_iterator itRes;
		itRes = aModel.find(mccore::ResId(it->c_str()));
		model.insert(*itRes);
	}
	return model;
}

void ExtractCyclesApp::writePDB(
	const std::map<annotate::CycleInfo, mccore::GraphModel>& aCycles) const
{
	std::map<annotate::CycleInfo, mccore::GraphModel>::const_iterator it;
	for(it = aCycles.begin(); it != aCycles.end(); ++ it)
	{
		std::ostringstream oss;
		ozfPdbstream outfile;
		oss << mstrOutputDirectory << "/";
		oss << it->first.getModelInfo().getPDBFile();
		oss << "." << it->first.getModelInfo().getModel();
		oss << "." << it->first.getProfile().toString();
		oss << "." << it->first.resIdsString();
		oss << ".pdb.gz";
		outfile.open(oss.str().c_str());
		outfile << it->second;
		outfile.close();
	}
}
