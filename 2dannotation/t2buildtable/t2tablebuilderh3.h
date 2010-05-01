#ifndef _t2tablebuilderh3_h_
#define _t2tablebuilderh3_h_

#include "t2tablebuilder.h"

class T2TableBuilderH3 : public T2TableBuilder
{
public:
	// LIFECYCLE
	T2TableBuilderH3(
		const std::string& astrGroupFile,
		const std::string& astrInteractions);

private:
	typedef std::pair<face, annotate::BasePair::enOrientation> oriented_face;
	typedef std::set<annotate::CycleInfo> cycle_set;
	std::vector<std::vector<unsigned int > > mInteractionsCount;
	std::map<oriented_face, unsigned int> mOrientedFaceMap;
	virtual void computeStatistics();
	void initializeTables();
	void countInteractions();
	void computeFrequencies();

	std::pair<unsigned int, unsigned int> getOrientedFacesIds(
		const annotate::InteractionInfo& aInteraction) const;
	std::pair<unsigned int, unsigned int> internalGetOrientedFacesIds(
		const oriented_face& aFace1, const oriented_face& aFace2) const;
	std::pair<unsigned int, unsigned int> getFacesIdFromTypeId(
		unsigned int auiTypeId) const;
};

#endif /*_t2tablebuilderh3_h_*/
