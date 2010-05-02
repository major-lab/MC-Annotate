#ifndef _t3tablebuilderh2_h_
#define _t3tablebuilderh2_h_

#include "t3tablebuilder.h"

class T3TableBuilderH2 : public T3TableBuilder
{
public:
	// LIFECYCLE
	T3TableBuilderH2(
		const std::string& astrGroupFile,
		const std::string& astrInteractions);

private:
	typedef std::set<annotate::CycleInfo> cycle_set;
	std::vector< std::vector<std::vector<unsigned int> > > mInteractionsCount;
	virtual void computeStatistics();
	void initializeTables();
	void countInteractions();
	void computeFrequencies();
};

#endif /*_t3tablebuilderh2_h_*/
