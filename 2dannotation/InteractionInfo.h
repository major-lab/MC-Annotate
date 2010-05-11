#ifndef _InteractionInfo_H_
#define _InteractionInfo_H_

#include <string>

#include "BasePair.h"

#include "ModelInfo.h"

namespace annotate {

class InteractionInfo
{
public:
	enum enFace
	{
		eWatson,
		eHoogsteen,
		eSugar,
		eRibose,
		ePhosphate
	};

	// LIFECYCLE ---------------------------------------------------------------
	InteractionInfo(
		const std::string& astrFile,
		unsigned int auiModel,
		const mccore::ResId& aResId1,
		const mccore::ResId& aResId2,
		const enFace& aFace1,
		const enFace& aFace2,
		const annotate::BasePair::enOrientation& aeOrientation);
	InteractionInfo(
		const std::string& astrFile,
		unsigned int auiModel,
		const std::string& astrRes1,
		const std::string& astrRes2,
		const enFace& aFace1,
		const enFace& aFace2,
		const annotate::BasePair::enOrientation& aeOrientation);
	~InteractionInfo() {}

	// ACCESS ------------------------------------------------------------------
	const ModelInfo& getModelInfo() const {return mModelInfo;}
	const mccore::ResId& resId1() const {return mResId1;}
	const mccore::ResId& resId2() const {return mResId2;}
	const enFace& face1() const {return meFace1;}
	const enFace& face2() const {return meFace2;}
	const annotate::BasePair::enOrientation& orientation() const {return meOrientation;}

	// METHODS -----------------------------------------------------------------
	std::set<char> getChains() const;
	void setChainAndOffset(char acChain, int aiOffset);

	// OPERATOR ----------------------------------------------------------------
	bool operator <(const InteractionInfo& aRight) const;
	bool operator ==(const InteractionInfo& aRight) const;

private:
	ModelInfo		mModelInfo;
	mccore::ResId	mResId1;
	mccore::ResId	mResId2;

	enFace meFace1;
	enFace meFace2;
	annotate::BasePair::enOrientation meOrientation;
};

}; // namespace annotate

#endif /*_InteractionInfo_H_*/
