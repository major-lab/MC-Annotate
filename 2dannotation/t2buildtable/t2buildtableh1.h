#ifndef _t2buildtableh1_h_
#define _t2buildtableh1_h_

#include "t2tablebuilder.h"

class T2BuildTableH1 : public T2TableBuilder
{
public:
	// LIFECYCLE
	T2BuildTableH1(
		const std::string& astrGroupFile,
		const std::string& astrInteractions);

private:
	typedef std::set<annotate::CycleInfo> cycle_set;
	std::vector<std::vector<std::vector<unsigned int> > > mInteractionsCount;
	virtual void computeStatistics();
	void initializeTables();
	void countInteractions();
	void computeFrequencies();
};

#endif /*_t2buildtableh1_h_*/
