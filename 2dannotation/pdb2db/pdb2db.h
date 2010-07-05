/*
 * pdb2db.h
 *
 *  Created on: June 30, 2010
 *      Author: blanchmf
 */

#ifndef _pdb2db_H_
#define _pdb2db_H_

#include <string>

#include "mccore/Molecule.h"

#include "AnnotationChains.h"
#include "AnnotationStems.h"

class PDB2DotBracket
{
public:
	PDB2DotBracket(int argc, char * argv []);

	// METHODS ---------------------------------------------
	void version () const;
	void usage () const;
	void help () const;
private:
	float mfAppVersion;
	std::string mstrAppName;
	std::string mstrPDBFile;

	void read_options (int argc, char* argv[]);
	mccore::Molecule* loadFile (const string &filename);

	std::string getFilePrefix(const std::string& aFileName) const;
	bool pseudoKnots(const std::vector<annotate::Stem>& usedStems, const annotate::Stem& otherStem) const;
	std::string toDotBracket(
		const annotate::AnnotationStems& aStems,
		const annotate::AnnotationChains& aChains) const;
};

#endif // _pdb2db_H_
