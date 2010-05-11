/*
 * CycleInfoFile.h
 *
 *  Created on: Feb 1, 2010
 *      Author: blanchmf
 */

#ifndef CYCLEINFOFILE_H_
#define CYCLEINFOFILE_H_

#include "CycleInfo.h"
#include <set>

namespace annotate
{

class CycleInfoFile
{
public:
	// LIFECYLE ----------------------------------------------------------------
	CycleInfoFile();
	~CycleInfoFile();

	// ACCESSOR ----------------------------------------------------------------
	const std::set<CycleInfo>& cycles() const {return mCycles;}

	// METHODS -----------------------------------------------------------------
	void read(const std::string& astrFilename, bool abFileProfile = true);


private:
	std::set<CycleInfo> mCycles;

	CycleInfo readLine(const std::string& astrLine, bool abFileProfile) const;
	CycleInfo::residue_profile getStrandResidues(
		const std::string& aResidues,
		const CycleProfile& aProfile) const;
	CycleInfo::residue_strand getResidues(const std::string& aResidues) const;
};

};

#endif /* CYCLEINFOFILE_H_ */
