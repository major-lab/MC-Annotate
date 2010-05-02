#ifndef _t3tablebuilderh1_h_
#define _t3tablebuilderh1_h_

#include "t3tablebuilder.h"

class T3BuildTableH1 : public T3TableBuilder
{
public:
	// LIFECYCLE
	T3BuildTableH1(
		const std::string& astrGroupFile,
		const std::string& astrInteractions);

private:
	typedef std::pair<unsigned int, unsigned int> interaction_coordinate;
	typedef std::set<annotate::CycleInfo> cycle_set;
	std::vector<std::vector< std::vector<std::vector<std::vector<unsigned int> > > > > mInteractionsCount;
	virtual void computeStatistics();
	void initializeTables();
	void countInteractions();
	void computeFrequencies();
	void computeFrequencies(
		unsigned int auiProfile1,
		unsigned int auiProfile2,
		unsigned int auiTypeId );
};

#endif /*_t3tablebuilderh1_h_*/
