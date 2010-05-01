#ifndef _t2buildtable_h_
#define _t2buildtable_h_

#include <string>

class T2BuildTable
{
public:
	// LIFECYCLE
	T2BuildTable(int argc, char * argv []);

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
#endif /*_t2buildtable_h_*/
