#ifndef _t4tablebuilderh1_h_
#define _t4tablebuilderh1_h_

#include "t4tablebuilder.h"

class T4BuildTableH1 : public T4TableBuilder
{
public:
	// LIFECYCLE
	T4BuildTableH1(
		const std::string& astrGroupFile,
		const std::string& astrInteractions);

private:
	std::vector<std::vector< std::vector<unsigned int> > > mInteractionsCount;
	virtual void computeStatistics();
	void initializeTables();
	void countInteractions();
	void computeFrequencies();
	unsigned int getNucleotideId(
		const mccore::ResId& aResId,
		const annotate::CycleInfo& aCycle) const;
};

#endif /*_t4tablebuilderh1_h_*/
