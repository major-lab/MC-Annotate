#ifndef _mcsc2resids_h_
#define _mcsc2resids_h_

#include <string>
#include <set>

namespace mccore
{
	class Molecule;
};

namespace annotate
{
	class Cycle;
};


// Prototypes
void version ();
void usage ();
void help ();
void read_options (int argc, char* argv[]);
mccore::Molecule* loadFile (const std::string &filename);
unsigned int getModelIndex(const std::string& aFileName);
std::string getFilePrefix(const std::string& aFileName);
std::string getPdbFileName(const std::string& aFileName);
annotate::Cycle mergeCycles(const std::set<annotate::Cycle>& aCycles);

#endif /*_mcsc2resids_h_*/
