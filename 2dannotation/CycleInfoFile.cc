/*
 * CycleInfoFile.cc
 *
 *  Created on: Feb 1, 2010
 *      Author: blanchmf
 */

#include "CycleInfoFile.h"

#include "StringUtil.h"

#include <cassert>
#include <fstream>

namespace annotate {

CycleInfoFile::CycleInfoFile()
{}

CycleInfoFile::~CycleInfoFile()
{
	mCycles.clear();
}

void CycleInfoFile::read(const std::string& astrFileName, bool abFileProfile)
{
	std::ifstream infile;
	infile.open(astrFileName.c_str(), std::ios_base::in);
	mCycles.clear();
	if(infile.good())
	{
		std::string strLine;
		while(std::getline(infile, strLine).good())
		{
			CycleInfo cInfo = readLine(strLine, abFileProfile);
			mCycles.insert(cInfo);
		}
	}
	else
	{
		// TODO : This should be an exception
		std::cout << "Error opening file " << astrFileName << std::endl;
	}
	infile.close();
}

/**
 * readLine
 * @brief Read a line from a CycleInfo file.
 * @return CycleInfo corresponding to that line.
 */
CycleInfo CycleInfoFile::readLine(const std::string& astrLine, bool abFileProfile) const
{
	std::string strLine = astrLine;
	cleanString(strLine, ' ');
	std::vector<std::string> fields = splitStringFields(strLine, ":");
	std::string strPDBFile = fields[0];
	unsigned int uiModel = atol(fields[1].c_str());
	std::string strProfile = fields[2];
	std::string strFileProfile;
	std::string strResIds;
	std::string strSeq;
	if(abFileProfile)
	{
		assert(6 == fields.size());
		strFileProfile = fields[3];
		strResIds = fields[4];
		strSeq = fields[5];
	}else
	{
		assert(5 == fields.size());
		strFileProfile = strProfile;
		strResIds = fields[3];
		strSeq = fields[4];
	}

	CycleProfile prof(strProfile);
	CycleProfile fileProf(strFileProfile);
	CycleInfo::residue_profile resProfile;
	resProfile = getStrandResidues(strResIds, fileProf);

	std::vector<std::string> residues = splitStringFields(strSeq, "-");
	CycleInfo returnVal(strPDBFile, uiModel, fileProf, prof, resProfile, residues);
	if(strProfile != strFileProfile)
	{
		returnVal = CycleInfo::flipStrand(returnVal);
	}
	return returnVal;
}

std::vector<std::vector<mccore::ResId> > CycleInfoFile::getStrandResidues(
	const std::string& aResidues,
	const annotate::CycleProfile& aProfile) const
{
	CycleInfo::residue_strand residues = getResidues(aResidues);
	CycleInfo::residue_profile strandResidues;

	CycleInfo::residue_strand::const_iterator itRes = residues.begin();
	std::list<unsigned int>::const_iterator itProf;
	for(itProf = aProfile.strandProfile().begin();
		itProf != aProfile.strandProfile().end();
		++ itProf)
	{
		CycleInfo::residue_strand strand;
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

CycleInfo::residue_strand CycleInfoFile::getResidues(
	const std::string& aResidues) const
{
	CycleInfo::residue_strand residues;
	std::string strResidues = aResidues;

	annotate::cleanString(strResidues, ' ');
	std::vector<std::string> residuesString = splitStringFields(strResidues, "-");
	std::vector<std::string>::const_iterator it = residuesString.begin();
	for(it = residuesString.begin(); it != residuesString.end(); ++it)
	{
		mccore::ResId res(it->c_str());
		residues.push_back(res);
	}
	return residues;
}

};
