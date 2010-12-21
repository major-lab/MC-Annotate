/*
 * pdb2db.h
 *
 *  Created on: June 30, 2010
 *      Author: blanchmf
 */

#ifndef _pdb2db_H_
#define _pdb2db_H_

#include <string>
#include <vector>

#include "mccore/Molecule.h"

class PDB2DotBracketParams
{
public:
	PDB2DotBracketParams()
	{
		mbOneModel = false;
		muiModelNumber = 0;
		mbBinary = false;
		muiCombinedLayers = 2;
		muiSplitLayers = 0;
		muiMaxPerfectSearch = 24;
		mbCompleteGaps = false;
	}

	// ATTRIBUTES --------------------------------------------------------------
	bool mbOneModel;
	unsigned int muiModelNumber;  // 1 based vector identifier, 0 means all
	unsigned int muiSplitLayers;
	unsigned int muiCombinedLayers;
	bool mbBinary;
	unsigned int muiMaxPerfectSearch; // Maximum number of stems for exhaustive search
	bool mbCompleteGaps;

	std::vector<std::string> mFiles;
};


class PDB2DotBracket
{
public:
	// LIFECYCLE -------------------------------------------
	PDB2DotBracket(
		const PDB2DotBracketParams& aParams,
		const std::vector<std::string>& aFiles);
	PDB2DotBracket(
		const PDB2DotBracketParams& aParams,
		const std::string& astrFile,
		std::ostream& aFile);
private:
	float mfAppVersion;
	std::string mstrAppName;

	PDB2DotBracketParams mParams;
	std::vector<std::string> mFiles;

	// void read_options (int argc, char* argv[]);
	mccore::Molecule* loadFile (const string &filename) const;
	mccore::Molecule* loadStream(const std::ostream& aPDBStream) const;

	void processStream(
		const std::string& astrFile,
		const std::ostream& aFile) const;
	void processFile(const std::string& astrFile) const;

	std::string getFilePrefix(const std::string& aFileName) const;
};

#endif // _pdb2db_H_
