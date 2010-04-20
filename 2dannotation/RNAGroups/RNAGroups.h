/*
 * RNAGroups.h
 *
 *  Created on: Mar 27, 2010
 *      Author: blanchmf
 */

#include <string>
#include <list>
#include <map>

namespace mccore {
	class Molecule;
	class Residue;
	class GraphModel;
};

class RNAGroups
{
public:
	RNAGroups(int argc, char * argv []);

	// METHODS ---------------------------------------------
	void version () const;
	void usage () const;
	void help () const;
private:
	float mfAppVersion;
	std::string mstrAppName;
	std::string mstrRNASelectFile;
	std::string mstrMoleculePath;
	std::map<unsigned int, std::list<std::string> > mGroups;
	std::map<std::string, unsigned int> mFileToGroup;

	void read_options (int argc, char* argv[]);
	void readRNASelectOutput();
	std::list<std::string> readGroupContent(std::ifstream& aInputFile);
	mccore::Molecule loadMolecule(const std::string &filename, bool abBinary);
	std::list<mccore::GraphModel> getRNAModels(const mccore::Molecule& aMolecule) const;
	std::list<mccore::Residue> getResidues(const mccore::GraphModel& aModel) const;
	void validateGroup(unsigned int auiGroup);
	std::string toString() const;
};
