/*
 * molsepRNA.h
 *
 *  Created on: Apr 07, 2010
 *      Author: blanchmf
 */

#include <string>

namespace mccore {
	class Molecule;
	class GraphModel;
	class AbstractModel;
};

class molsepRNA
{
public:
	molsepRNA(int argc, char * argv []);

	// METHODS ---------------------------------------------
	void version () const;
	void usage () const;
	void help () const;
private:
	float mfAppVersion;
	std::string mstrAppName;
	std::string mstrOutputPath;
	bool mbSplitChains;

	void read_options (int argc, char* argv[]);
	mccore::Molecule loadMolecule(const std::string &filename, bool abBinary);
	std::string getFilePrefix(const std::string& astrFile) const;
	void writeModel(
		const std::string astrOriginalFile,
		unsigned int auiModel,
		const mccore::AbstractModel& aModel) const;
	void writeChains(
		const std::string astrOriginalFile,
		unsigned int auiModel,
		const mccore::AbstractModel& aModel) const;
};
