#ifndef _mcsc2resids_h_
#define _mcsc2resids_h_

#include <list>
#include <set>
#include <string>
#include <vector>


namespace mccore
{
	class Molecule;
	class ResId;
};

namespace annotate
{
	class AnnotateModel;
	class BaseInteraction;
	class BaseLink;
	class BasePair;
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
const annotate::BaseLink* findFirstBaseLink(
	std::list<const annotate::BaseInteraction*>& aInteractions);
const annotate::BasePair* findFirstBasePair(
	std::list<const annotate::BaseInteraction*>& aInteractions);
annotate::Cycle createCycleFromStrands(
	const annotate::AnnotateModel& aModel, 
	std::vector<std::vector<mccore::ResId> > aStrands);

#endif /*_mcsc2resids_h_*/
