#ifndef _t2tablebuilderh2_h_
#define _t2tablebuilderh2_h_

#include "t2tablebuilder.h"

class T2TableBuilderH2 : public T2TableBuilder
{
public:
	// LIFECYCLE
	T2TableBuilderH2(
		const std::string& astrGroupFile,
		const std::string& astrInteractions);

private:
	typedef std::set<annotate::CycleInfo> cycle_set;
	std::vector<std::vector<unsigned int > > mInteractionsCount;
	virtual void computeStatistics();
	void initializeTables();
	void countInteractions();
	void computeFrequencies();
};

#endif /*_t2tablebuilderh2_h_*/
