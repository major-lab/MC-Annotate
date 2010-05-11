/*
 * InteractionInfoFile.cc
 *
 *  Created on: Mar 10, 2010
 *      Author: blanchmf
 */

#include "InteractionInfoFile.h"

#include "StringUtil.h"

#include <cassert>
#include <cstdlib>
#include <fstream>

namespace annotate {

InteractionInfoFile::InteractionInfoFile()
{}

InteractionInfoFile::~InteractionInfoFile()
{
	mInteractions.clear();
}

void InteractionInfoFile::read(const std::string& astrFileName)
{
	std::ifstream infile;
	mInteractions.clear();

	infile.open (astrFileName.c_str(), std::ifstream::in);
	std::string strLine;

	while (std::getline(infile, strLine).good())
	{
		annotate::cleanString(strLine, ' ');

		std::vector<std::string> fields = splitStringFields(strLine, ":");

		// Read the File
		std::string strPDBFile = fields[0];

		// Read the model
		unsigned int uiModel = atol(fields[1].c_str());

		// Read the residues
		std::vector<std::string> residues = splitStringFields(fields[2], "-");

		// Read ( and drop ) the complete pairing info
		// std::string strFaces = fields[3];

		// Read the generalized pairing info
		std::vector<std::string> generalFaces = splitStringFields(fields[4], "/");
		InteractionInfo::enFace eFace1 = getFace(generalFaces[0]);
		InteractionInfo::enFace eFace2 = getFace(generalFaces[1]);

		//  Read the relative orientation
		std::string strOrientation = fields[5];
		annotate::BasePair::enOrientation eOrientation = getOrientation(strOrientation);

		// Create the instance of interaction info
		annotate::InteractionInfo info(
			strPDBFile,
			uiModel,
			residues[0],
			residues[1],
			eFace1,
			eFace2,
			eOrientation);
		mInteractions.insert(info);
	}
	infile.close();
}

BasePair::enOrientation InteractionInfoFile::getOrientation(const std::string& astrOrientation) const
{
	annotate::BasePair::enOrientation eOrientation = annotate::BasePair::eUnknown;
	if(astrOrientation == "Cis")
	{
		eOrientation = annotate::BasePair::eCis;
	}else if(astrOrientation == "Trans")
	{
		eOrientation = annotate::BasePair::eTrans;
	}else if(astrOrientation == "Unknown")
	{
		eOrientation = annotate::BasePair::eUnknown;
	}else
	{
		assert(false); // There is a problem with the file
	}
	return eOrientation;
}

InteractionInfo::enFace InteractionInfoFile::getFace(
	const std::string& astrGeneralFaces) const
{
	InteractionInfo::enFace eFace = InteractionInfo::eWatson;
	if(astrGeneralFaces == "W")
	{
		eFace = InteractionInfo::eWatson;
	}else if(astrGeneralFaces == "H")
	{
		eFace = InteractionInfo::eHoogsteen;
	}else if(astrGeneralFaces == "S")
	{
		eFace = InteractionInfo::eSugar;
	}else if(astrGeneralFaces == "Ribose")
	{
		eFace = InteractionInfo::eRibose;
	}
	else if(astrGeneralFaces == "Phosphate")
	{
		eFace = InteractionInfo::ePhosphate;
	}else
	{
		std::cerr << "Couldn't convert face string : " << astrGeneralFaces << std::endl;
		assert(false); // There is a problem with the file
	}
	return eFace;
}

};
