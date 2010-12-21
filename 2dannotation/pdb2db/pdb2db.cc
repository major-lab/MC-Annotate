//                              -*- Mode: C++ -*-
// pdb2db.cc
// Copyright © 2001-10 Laboratoire d'ingénierie des ARN.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Thu Jul 8 10:00:00 2010


#include "pdb2db.h"

#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <sstream>

#include "mccore/Binstream.h"
#include "mccore/Exception.h"
#include "mccore/Messagestream.h"
#include "mccore/ModelFactoryMethod.h"
#include "mccore/Pdbstream.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"
#include "mccore/ResidueFactoryMethod.h"
#include "mccore/ResIdSet.h"

// libmcannotate
#include "AlgorithmExtra.h"
#include "AnnotateModel.h"
#include "AnnotationChains.h"
#include "AnnotationInteractions.h"


#include "AnnotationStemsLoose.h"
#include "ModelDotBracketAnnotator.h"


PDB2DotBracket::PDB2DotBracket(
	const PDB2DotBracketParams& aParams,
	const std::vector<std::string>& aFiles)
{
	mParams = aParams;

	std::vector<std::string>::const_iterator itFile;
	for(itFile = aFiles.begin(); itFile != aFiles.end(); ++ itFile)
	{
		processFile(*itFile);
	}
}

PDB2DotBracket::PDB2DotBracket(
	const PDB2DotBracketParams& aParams,
	const std::string& astrFile,
	std::ostream& aFile)
{
	mParams = aParams;
	processStream(astrFile, aFile);
}

void PDB2DotBracket::processFile(const std::string& astrFile) const
{
	mccore::Molecule *molecule;
	mccore::Molecule::iterator molIt;
	unsigned int uiModelNumber = mParams.muiModelNumber;

	molecule = loadFile (astrFile);
	if (0 != molecule)
	{
		unsigned int uiCurrentModel = 1;
		for (molIt = molecule->begin (); molecule->end () != molIt; ++molIt)
		{
			if (0 != uiModelNumber)
			{
				--uiModelNumber;
			}
			else
			{
				ModelDotBracketAnnotator dbAnnotator(
					mParams.mbCompleteGaps,
					mParams.muiCombinedLayers,
					mParams.muiSplitLayers,
					mParams.muiMaxPerfectSearch);
				annotate::AnnotateModel &am = (annotate::AnnotateModel&) *molIt;
				std::string strPDBPrefix = getFilePrefix(astrFile);
				am.id(uiCurrentModel);
				am.name(strPDBPrefix);
				dbAnnotator.model(am);
				dbAnnotator.process();
				if (mParams.mbOneModel)
				{
					break;
				}
			}
			uiCurrentModel ++;
		}
		delete molecule;
	}
}

void PDB2DotBracket::processStream(
	const std::string& astrFile,
	const std::ostream& aFile) const
{
	mccore::Molecule *molecule;
	mccore::Molecule::iterator molIt;
	unsigned int uiModelNumber = mParams.muiModelNumber;

	molecule = loadStream(aFile);
	if (0 != molecule)
	{
		unsigned int uiCurrentModel = 1;
		for (molIt = molecule->begin (); molecule->end () != molIt; ++molIt)
		{
			if (0 != uiModelNumber)
			{
				--uiModelNumber;
			}
			else
			{
				ModelDotBracketAnnotator dbAnnotator(
					mParams.mbCompleteGaps,
					mParams.muiCombinedLayers,
					mParams.muiSplitLayers,
					mParams.muiMaxPerfectSearch);
				annotate::AnnotateModel &am = (annotate::AnnotateModel&) *molIt;
				std::string strPDBPrefix = getFilePrefix(astrFile);
				am.id(uiCurrentModel);
				am.name(strPDBPrefix);
				dbAnnotator.model(am);
				dbAnnotator.process();
				if (mParams.mbOneModel)
				{
					break;
				}
			}
			uiCurrentModel ++;
		}
		delete molecule;
	}
}

mccore::Molecule* PDB2DotBracket::loadStream(const std::ostream& aPDBStream) const
{
	Molecule *molecule = 0;
	ResidueFM rFM;
	ResIdSet residueSelection;
	annotate::AnnotateModelFM aFM (residueSelection, 0, &rFM);
	// TODO : Check how to make this work for compressed PDB
	iPdbstream in(aPDBStream.rdbuf());
	molecule = new Molecule (&aFM);
	in >> *molecule;
	return molecule;
}

mccore::Molecule* PDB2DotBracket::loadFile (const string &filename) const
{
	Molecule *molecule;
	ResidueFM rFM;
	ResIdSet residueSelection;
	annotate::AnnotateModelFM aFM (residueSelection, 0, &rFM);

	molecule = 0;
	if (mParams.mbBinary)
	{
		izfBinstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			mccore::gErr (0) << PACKAGE << ": cannot open binary file '" << filename << "'." << endl;
			return 0;
		}
		molecule = new Molecule (&aFM);
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
		izfPdbstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			mccore::gErr (0) << PACKAGE << ": cannot open pdb file '" << filename << "'." << endl;
			return 0;
		}
		molecule = new Molecule (&aFM);
		in >> *molecule;
		in.close ();
#ifdef HAVE_LIBRNAMLC__
		}
#endif
	}
	return molecule;
}

std::string PDB2DotBracket::getFilePrefix(const std::string& aFileName) const
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
