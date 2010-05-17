#ifndef _mcnaip_h_
#define _mcnaip_h_

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
std::string getPdbFileName(const std::string& aFileName);

#endif /*_mcnaip_h_*/
