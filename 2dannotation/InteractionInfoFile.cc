/*
 * InteractionInfoFile.cc
 *
 *  Created on: Mar 10, 2010
 *      Author: blanchmf
 */

#include "InteractionInfoFile.h"

#include "StringUtil.h"

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

		std::size_t sep = strLine.rfind('-');
		std::string strRes2 = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);

		sep = strLine.rfind(':');
		std::string strRes1 = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);

		sep = strLine.rfind(':');
		std::string strModel = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);
		unsigned int uiModel = atol(strModel.c_str());

		std::string strPDBFile = strLine;

		annotate::InteractionInfo info(strPDBFile, uiModel, strRes1, strRes2);
		mInteractions.insert(info);
	}
	infile.close();
}

};
