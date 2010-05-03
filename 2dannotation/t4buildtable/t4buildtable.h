#ifndef _t4buildtable_h_
#define _t4buildtable_h_

#include <string>

class T4BuildTable
{
public:
	// LIFECYCLE
	T4BuildTable(int argc, char * argv []);

	// Console specific methods
	void version() const;
	void usage() const;
	void help() const;
private:
	void readOptions(int argc, char* argv[]);
	std::string mstrInteractionsFile;
	std::string mstrPDBGroupsFile;
	int miHypothesis;
};
#endif /*_t4buildtable_h_*/
