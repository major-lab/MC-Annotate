/*
 * CycleInfoFile.cc
 *
 *  Created on: Feb 1, 2010
 *      Author: blanchmf
 */

#include "CycleInfoFile.h"

#include "StringUtil.h"

#include <fstream>

namespace annotate {

CycleInfoFile::CycleInfoFile()
{}

CycleInfoFile::~CycleInfoFile()
{
	mCycles.clear();
}

void CycleInfoFile::read(const char* aszFileName)
{
	std::ifstream infile;
	infile.open(aszFileName, std::ios_base::in);
	mCycles.clear();
	if(infile.good())
	{
		std::string strLine;
		while(std::getline(infile, strLine).good())
		{
			CycleInfo cInfo = readLine(strLine);
			mCycles.insert(cInfo);
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
 * @brief Read a line from a CycleInfo file.
 * @return CycleInfo corresponding to that line.
 */
CycleInfo CycleInfoFile::readLine(const std::string& astrLine) const
{
	std::string strLine = astrLine;
	cleanString(strLine, ' ');
	std::vector<std::string> fields = splitStringFields(strLine, ":");
	std::string strPDBFile = fields[0];
	unsigned int uiModel = atol(fields[1].c_str());
	std::string strFileProfile = fields[2];
	std::string strProfile = fields[3];
	std::string strResIds = fields[4];
	std::string strSeq = fields[5];

	CycleProfile prof(strProfile);
	CycleProfile fileProf(strFileProfile);
	CycleInfo::residue_profile resProfile;
	resProfile = getStrandResidues(strResIds, fileProf);

	std::vector<std::string> residues = splitStringFields(strSeq, "-");
	return CycleInfo(strPDBFile, uiModel, fileProf, prof, resProfile, residues);
}

std::vector<std::vector<std::string> > CycleInfoFile::getStrandResidues(
	const std::string& aResidues,
	const annotate::CycleProfile& aProfile) const
{
	std::list<std::string> residues = getResidues(aResidues);
	std::vector<std::vector<std::string> > strandResidues;

	std::list<std::string>::const_iterator itRes = residues.begin();
	std::list<unsigned int>::const_iterator itProf;
	for(itProf = aProfile.strandProfile().begin();
		itProf != aProfile.strandProfile().end();
		++ itProf)
	{
		std::vector<std::string> strand;
		unsigned int iRes = 0;
		for(iRes = 0; iRes < *itProf; ++ iRes)
		{
			strand.push_back(*itRes);
			++ itRes;
		}
		strandResidues.push_back(strand);
		strand.clear();
	}
	return strandResidues;
}

std::list<std::string> CycleInfoFile::getResidues(
	const std::string& aResidues) const
{
	std::list<std::string> residues;
	std::string strResidues = aResidues;

	annotate::cleanString(strResidues, ' ');
	while(0 < strResidues.size())
	{
		std::string strRes;
		std::size_t sep = strResidues.rfind('-');
		if(sep != std::string::npos)
		{
			strRes = strResidues.substr(sep + 1, strResidues.size() - (sep + 1));
			strResidues.erase(sep);
		}
		else
		{
			strRes = strResidues;
			strResidues.clear();
		}
		residues.push_front(strRes);
	}
	return residues;
}

};
