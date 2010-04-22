/*
 * RNAGroupFile.cc
 *
 *  Created on: Apr 22, 2010
 *      Author: blanchmf
 */

#include "RNAGroupFile.h"

#include "StringUtil.h"

#include <cassert>
#include <fstream>

namespace annotate {

RNAGroupFile::RNAGroupFile()
{}

RNAGroupFile::~RNAGroupFile()
{
	mGroups.clear();
}

void RNAGroupFile::read(const char* aszFileName)
{
	std::ifstream infile;
	infile.open(aszFileName, std::ios_base::in);
	mGroups.clear();
	if(infile.good())
	{
		std::string strLine;
		while(std::getline(infile, strLine).good())
		{
			readLine(strLine);
			}
	}
	else
	{
		// TODO : This should be an exception
		std::cout << "Error opening file " << aszFileName << std::endl;
	}
	infile.close();
}

/**
 * readLine
 * @brief Read a line from an RNAGroup file.
  */
void RNAGroupFile::readLine(const std::string& astrLine) const
{
	std::string strLine = astrLine;
	cleanString(strLine, ' ');
	std::vector<std::string> fields = splitStringFields(strLine, ":");
	unsigned int uiGroup = atol(fields[0].c_str());
	int iOffset = atoi(fields[1].c_str());

	// Read the files
	for(unsigned int i = 0; i < fields.size(); ++ i)
	{
		std::vector<std::string> fileFields = splitStringFields(fields[i], "_");
		stGroupFileEntry entry;
		entry.strName = fileFields[0];
		entry.uiModel = atol(fields[1].c_str());
		entry.cChain = fileFields[0][0];
		entry.iOffset = iOffset;
		mGroups[uiGroup].push_back(entry);
	}
	return returnVal;
}

};

