#include "Interaction.h"

// libmcannotate
#include "AlgorithmExtra.h"
#include "CycleInfo.h"

#include <cassert>
#include <sstream>

namespace annotate {

// TODO : Remove this (currently required to allow empty map initialization)
CycleInfo::CycleInfo()
: 	mModelInfo("", 0),
	mProfile("2_2"),
	mFileProfile("2_2")
{
}

CycleInfo::CycleInfo(
	const std::string& aFile,
	unsigned int auiModel,
	const annotate::CycleProfile& aFileProfile,
	const annotate::CycleProfile& aProfile,
	const residue_profile& aResIds,
	const std::vector<std::string>& aResidues)
: 	mModelInfo(aFile, auiModel),
	mProfile(aProfile),
	mFileProfile(aFileProfile)
{
	unsigned int uiNbRes = 0;
	for(unsigned int uiStrand = 0; uiStrand < aResIds.size(); ++ uiStrand)
	{
		uiNbRes += aResIds[uiStrand].size();
	}
	mResidueIds.resize(uiNbRes);
	unsigned int uiIndex = 0;
	for(unsigned int uiStrand = 0; uiStrand < aResIds.size(); ++ uiStrand)
	{
		for(unsigned int uiRes = 0; uiRes < aResIds[uiStrand].size(); ++ uiRes)
		{
			mResidueIds[uiIndex] = aResIds[uiStrand][uiRes];
			mIdToResMap.insert(std::pair<mccore::ResId, std::string>(aResIds[uiStrand][uiRes], aResidues[uiIndex]));
			uiIndex ++;
		}
	}
	assert(mIdToResMap.size() == mResidueIds.size());
}

CycleInfo::CycleInfo(
		const std::string& aFile,
		unsigned int auiModel,
		const annotate::CycleProfile& aFileProfile,
		const annotate::CycleProfile& aProfile,
		const std::vector<mccore::ResId>& aResIds,
		const std::vector<std::string>& aResidues)
	: 	mModelInfo(aFile, auiModel),
		mProfile(aProfile),
		mFileProfile(aFileProfile)
{
	mResidueIds = aResIds;
	for(unsigned int uiRes = 0; uiRes < mResidueIds.size(); ++ uiRes)
	{
		mIdToResMap.insert(std::pair<mccore::ResId, std::string>(mResidueIds[uiRes], aResidues[uiRes]));
	}
	assert(mIdToResMap.size() == mResidueIds.size());
}

bool CycleInfo::operator <(const CycleInfo& aRight) const
{
	bool bLess = false;

	if(mModelInfo < aRight.mModelInfo)
	{
		bLess = true;
	}
	else if(mModelInfo == aRight.mModelInfo)
	{
		if(mProfile < aRight.mProfile)
		{
			bLess = true;
		}else if(mProfile == aRight.mProfile)
		{
			for(unsigned int uiIndex = 0; uiIndex < mResidueIds.size(); ++ uiIndex)
			{
				if(mResidueIds[uiIndex] < aRight.mResidueIds[uiIndex])
				{
					bLess = true;
					break;
				}else if(aRight.mResidueIds[uiIndex] < mResidueIds[uiIndex])
				{
					// Greather than
					break;
				}
			}
		}
	}
	return bLess;
}

bool CycleInfo::operator ==(const CycleInfo& aRight) const
{
	bool bEqual = true;

	if(mModelInfo == aRight.mModelInfo)
	{
		if(mProfile == aRight.mProfile)
		{
			for(unsigned int uiIndex = 0; uiIndex < mResidueIds.size(); ++ uiIndex)
			{
				if(!(mResidueIds[uiIndex] == aRight.mResidueIds[uiIndex]))
				{
					bEqual = false;
					break;
				}
			}
		}else
		{
			bEqual = false;
		}
	}
	return bEqual;
}

// TODO : Check if this can be removed
std::pair<int, int> CycleInfo::findResId(const mccore::ResId& astrResId) const
{
	std::pair<int, int> coord(-1, -1);
	for(int iRow = 0; iRow < (int)mProfile.strandProfile().size(); ++ iRow)
	{
		std::vector<mccore::ResId> resids = getStrandResIds(iRow);
		for(int iCol = 0; iCol < (int)resids.size(); ++ iCol)
		{
			if(resids[iCol] == astrResId)
			{
				coord.first = iRow;
				coord.second = iCol;
				break;
			}
		}
	}
	return coord;
}

bool CycleInfo::contains(const InteractionInfo& aInteraction) const
{
	bool bContains = false;
	std::vector<mccore::ResId> resIds = getResIds();
	if(2 < resIds.size())
	{
		unsigned int i;
		for(i = 0; i < resIds.size(); ++ i)
		{
			if(resIds[i] == aInteraction.resId1())
			{
				break;
			}
		}
		if(i < resIds.size())
		{
			mccore::ResId res2 = aInteraction.resId2();
			if(0 == i)
			{
				if(resIds[i + 1] == res2)
				{
					bContains = true;
				}
				else if(resIds.back() == res2)
				{
					bContains = true;
				}
			}
			else if(resIds.size() - 1 == i)
			{
				if(resIds.front() == res2)
				{
					bContains = true;
				}
				else if(resIds[i - 1] == res2)
				{
					bContains = true;
				}
			}
			else
			{
				if(resIds[i + 1] == res2 || resIds[i - 1] == res2)
				{
					bContains = true;
				}
			}
		}
	}
	return bContains;
}

bool CycleInfo::contains(const mccore::ResId& aResId) const
{
	bool bContains = false;
	std::vector<mccore::ResId> resIds = getResIds();
	unsigned int i;
	for(i = 0; i < resIds.size() && !bContains; ++ i)
	{
		if(resIds[i] == aResId)
		{
			bContains = true;
		}
	}
	return bContains;
}

std::set<Interaction> CycleInfo::getStrandInteractions(
	unsigned int auiStrand) const
{
	std::set<Interaction> interactions;

	assert(auiStrand < mProfile.strandProfile().size());

	std::vector<mccore::ResId>::const_iterator it;
	std::vector<mccore::ResId>::const_iterator itPrev;
	std::vector<mccore::ResId> resids = getStrandResIds(auiStrand);
	for(it = resids.begin(); it != resids.end(); ++ it)
	{
		if(it != resids.begin())
		{
			interactions.insert(Interaction(*itPrev, *it));
		}
		itPrev = it;
	}
	return interactions;
}

bool CycleInfo::shareInteraction(
	const std::set<Interaction>& aInteractions) const
{
	bool bShares = false;
	for(unsigned int uiStrand = 0;
		uiStrand < mProfile.strandProfile().size() && !bShares;
		++ uiStrand)
	{
		std::set<Interaction> interactions = getStrandInteractions(uiStrand);
		bShares = annotate::set_intersects(
			interactions.begin(), interactions.end(),
			aInteractions.begin(), aInteractions.end());
	}
	return bShares;
}

bool CycleInfo::isSubCycleOf(const CycleInfo& aCycleInfo) const
{
	bool bSubCycle = false;
	unsigned int uiNbStrands = mProfile.strandProfile().size();
	unsigned int uiRightNbStrands = aCycleInfo.mProfile.strandProfile().size();

	if(1 == uiNbStrands && 1 == uiRightNbStrands)
	{
		// Loop
		bSubCycle = isSubLoopOf(aCycleInfo);
	}else if(2 == uiNbStrands && 2 == uiRightNbStrands)
	{
		// 2 Strand
		bSubCycle = isSub2StrandsCycle(aCycleInfo);
	}
	else
	{
		assert(1 == uiNbStrands || 2 == uiNbStrands);
		assert(1 == uiRightNbStrands || 2 == uiRightNbStrands);

		if(2 == uiNbStrands)
		{
			bSubCycle = is2StrandsSubCycleOfLoop(aCycleInfo);
		}else
		{
			bSubCycle = isLoopSubCycleOf2Strands(aCycleInfo);
		}

	}
	return bSubCycle;
}

bool CycleInfo::isSubLoopOf(const CycleInfo& aCycleInfo) const
{
	bool bSubLoop = false;
	assert(mProfile.strandProfile().size() == 1);
	assert(aCycleInfo.mProfile.strandProfile().size() == 1);

	if(mProfile.strandProfile()[0] < aCycleInfo.mProfile.strandProfile()[0])
	{
		std::set<Interaction> interactionsLeft;
		std::set<Interaction> interactionsRight;
		interactionsLeft = getStrandInteractions(0);
		interactionsRight = aCycleInfo.getStrandInteractions(0);

		std::set<Interaction> intersection;
		intersection = annotate::SetIntersection(
			interactionsLeft,
			interactionsRight);
		bSubLoop = (intersection.size() == interactionsLeft.size());
	}

	return bSubLoop;
}

bool CycleInfo::isSub2StrandsCycle(const CycleInfo& aCycleInfo) const
{
	bool bSubLoop = false;
	unsigned int uiNbStrands = mProfile.strandProfile().size();
	unsigned int uiRightNbStrands = aCycleInfo.mProfile.strandProfile().size();
	assert(uiNbStrands == 2);
	assert(uiRightNbStrands == 2);

	unsigned int iStrand;
	for(iStrand = 0; iStrand < uiNbStrands; ++ iStrand)
	{
		if(mProfile.strandProfile()[iStrand] <= aCycleInfo.mProfile.strandProfile()[iStrand])
		{
			std::set<Interaction> interactionsLeft;
			std::set<Interaction> interactionsRight;
			interactionsLeft = getStrandInteractions(iStrand);
			interactionsRight = aCycleInfo.getStrandInteractions(iStrand);

			std::set<Interaction> intersection;
			intersection = annotate::SetIntersection(
				interactionsLeft,
				interactionsRight);
			if(intersection.size() != interactionsLeft.size())
			{
				break;
			}
		}
		else
		{
			break;
		}
	}

	// If we got through, than it is a sub cycle
	bSubLoop = (iStrand == uiNbStrands);

	return bSubLoop;
}

bool CycleInfo::is2StrandsSubCycleOfLoop(const CycleInfo& aCycleInfo) const
{
	bool bSubLoop = false;
	unsigned int uiNbStrands = mProfile.strandProfile().size();
	unsigned int uiRightNbStrands = aCycleInfo.mProfile.strandProfile().size();
	assert(uiNbStrands == 2);
	assert(uiRightNbStrands == 1);


	std::set<Interaction> interactionsRight;
	interactionsRight = aCycleInfo.getStrandInteractions(0);

	std::set<Interaction> interactionsLeft;
	interactionsLeft = getStrandInteractions(0);
	std::set<Interaction> interactionStrand2 = getStrandInteractions(1);
	interactionsLeft.insert(interactionStrand2.begin(), interactionStrand2.end());

	if(interactionsLeft.size() < interactionsRight.size())
	{
		std::set<Interaction> intersection;
		intersection = annotate::SetIntersection(
			interactionsLeft,
			interactionsRight);
		bSubLoop = (intersection.size() == interactionsLeft.size());
	}

	return bSubLoop;
}

bool CycleInfo::isLoopSubCycleOf2Strands(const CycleInfo& aCycleInfo) const
{
	bool bSubLoop = false;
	unsigned int uiNbStrands = mProfile.strandProfile().size();
	unsigned int uiRightNbStrands = aCycleInfo.mProfile.strandProfile().size();
	assert(uiNbStrands == 1);
	assert(uiRightNbStrands == 2);

	std::set<Interaction> interactionsRight;
	interactionsRight = aCycleInfo.getStrandInteractions(0);
	std::set<Interaction> interactionStrand2 = aCycleInfo.getStrandInteractions(1);
	interactionsRight.insert(interactionStrand2.begin(), interactionStrand2.end());

	std::set<Interaction> interactionsLeft;
	interactionsLeft = getStrandInteractions(0);

	if(interactionsLeft.size() < interactionsRight.size())
	{
		std::set<Interaction> intersection;
		intersection = annotate::SetIntersection(
			interactionsLeft,
			interactionsRight);
		bSubLoop = (intersection.size() == interactionsLeft.size());
	}

	return bSubLoop;
}

bool CycleInfo::strandCoversInteractions(
	unsigned int uiStrand,
	const std::set<Interaction>& aInteractions) const
{
	bool bCovers;
	std::set<Interaction> strandInteractions = getStrandInteractions(uiStrand);

	std::set<Interaction> intersection = annotate::SetIntersection(
		aInteractions,
		strandInteractions);

	bCovers = (aInteractions.size() == intersection.size());
	return bCovers;
}

bool CycleInfo::hasStrandCoveringInteractions(
	const std::set<Interaction>& aInteractions) const
{
	bool bMatch = false;
	unsigned int uiStrand = 0;
	for(uiStrand = 0; uiStrand < mProfile.strandProfile().size() && !bMatch; ++ uiStrand)
	{
		bMatch = strandCoversInteractions(uiStrand, aInteractions);
	}
	return bMatch;
}

const std::vector<std::string> CycleInfo::getSequence() const
{
	std::vector<std::string> sequence;
	for(unsigned int i = 0; i < mResidueIds.size(); ++ i)
	{
		std::map<mccore::ResId, std::string>::const_iterator it;
		it = mIdToResMap.find(mResidueIds[i]);
		assert(it != mIdToResMap.end());
		sequence.push_back(it->second);
	}
	return sequence;
}

std::string CycleInfo::getNucleotideString(const mccore::ResId& aResId) const
{
	std::map<mccore::ResId, std::string>::const_iterator it;
	it = mIdToResMap.find(aResId);
	assert(it != mIdToResMap.end());
	return it->second;
}

std::string CycleInfo::toString(const std::string astrSeparator) const
{
	std::ostringstream oss;

	oss << mProfile.toString();
	oss << " " << astrSeparator << " ";
	oss << groupedResIdsString();
	oss << " " << astrSeparator << " ";
	oss << residuesString();
	return oss.str();
}

std::string CycleInfo::residuesString(const std::string& astrSeparator) const
{
	std::ostringstream oss;

	std::vector<std::string> sequence = getSequence();
	std::vector<std::string>::const_iterator it;
	for(it = sequence.begin(); it != sequence.end(); ++ it)
	{
		if(it != sequence.begin())
		{
			oss << astrSeparator;
		}
		oss << *it;
	}
	return oss.str();
}

std::string CycleInfo::resIdsString(const std::string& astrSeparator) const
{
	std::ostringstream oss;

	for(unsigned int i = 0; i < mResidueIds.size(); ++ i)
	{
		if(0 < i)
		{
			oss << astrSeparator;
		}
		oss << mResidueIds[i];
	}
	return oss.str();
}

std::string CycleInfo::groupedResIdsString() const
{
	std::ostringstream oss;

	unsigned int uiNbStrands = mProfile.strandProfile().size();
	CycleInfo::residue_profile::const_iterator itResStrand;
	for(unsigned int uiStrand = 0; uiStrand < uiNbStrands; ++ uiStrand)
	{
		if(0 != uiStrand)
		{
			oss << ",";
		}
		oss << "{";
		std::vector<mccore::ResId> res = getStrandResIds(uiStrand);
		for(unsigned int uiRes = 0; uiRes < res.size(); ++ uiRes)
		{
			if(0 != uiRes)
			{
				oss << "-";
			}
			oss << res[uiRes];
		}
		oss << "}";
	}

	return oss.str();
}

const std::vector<std::string> CycleInfo::getFlipStrandSequence() const
{
	assert(2 == mProfile.strandProfile().size());
	std::vector<std::string> sequence;
	residue_strand::const_iterator itRes;
	std::vector<mccore::ResId> resIds1 = getStrandResIds(0);
	std::vector<mccore::ResId> resIds2 = getStrandResIds(1);
	for(itRes = resIds2.begin(); itRes != resIds2.end(); ++ itRes)
	{
		std::map<mccore::ResId, std::string>::const_iterator it;
		it = mIdToResMap.find(*itRes);
		assert(it != mIdToResMap.end());
		sequence.push_back(it->second);
	}
	for(itRes = resIds1.begin(); itRes != resIds1.end(); ++ itRes)
	{
		std::map<mccore::ResId, std::string>::const_iterator it;
		it = mIdToResMap.find(*itRes);
		assert(it != mIdToResMap.end());
		sequence.push_back(it->second);
	}
	return sequence;
}

CycleInfo CycleInfo::flipStrand(const CycleInfo& aCycleInfo)
{
	CycleInfo newCycle;

	if(aCycleInfo.getNbStrands() == 2)
	{
		std::vector<std::vector<mccore::ResId> > flipResidueProfile;
		std::vector<std::string> flipResidues;

		flipResidueProfile.resize(2);
		flipResidueProfile[0] = aCycleInfo.getStrandResIds(1);
		flipResidueProfile[1] = aCycleInfo.getStrandResIds(0);

		flipResidues = aCycleInfo.getFlipStrandSequence();

		newCycle = CycleInfo(
			aCycleInfo.getPDBFile(),
			aCycleInfo.getModel(),
			aCycleInfo.mFileProfile,
			aCycleInfo.mProfile,
			flipResidueProfile,
			flipResidues);
	}
	else
	{
		newCycle = aCycleInfo;
	}

	return newCycle;
}

std::set<char> CycleInfo::getChains() const
{
	std::set<char> chains;
	for(unsigned int i = 0; i < mResidueIds.size(); ++ i)
	{
		chains.insert(mResidueIds[i].getChainId());
	}
	return chains;
}

void CycleInfo::setChainAndOffset(char acChain, int aiOffset)
{
	for(unsigned int i = 0; i < mResidueIds.size(); ++ i)
	{
		mccore::ResId& res = mResidueIds[i];
		res.setChainId(acChain);
		res.setResNo(res.getResNo() + aiOffset);
	}
}

std::pair<unsigned int, unsigned int> CycleInfo::getStrandRange(unsigned int auiStrand) const
{
	std::pair<unsigned int, unsigned int> range(0,0);
	unsigned int uiNbStrand = mProfile.strandProfile().size();
	assert(auiStrand < uiNbStrand);

	unsigned int uiStartIndex = 0;

	for(unsigned int i = 0; i < uiNbStrand; ++ i)
	{
		if(i != auiStrand)
		{
			uiStartIndex += mProfile.strandProfile()[i];
		}else
		{
			break;
		}
	}
	range.first = uiStartIndex;
	range.second = uiStartIndex + mProfile.strandProfile()[auiStrand];
	return range;
}

std::vector<mccore::ResId> CycleInfo::getStrandResIds(unsigned int auiStrand) const
{
	std::vector<mccore::ResId> resids;
	std::pair<unsigned int, unsigned int> range = getStrandRange(auiStrand);
	resids.resize(range.second - range.first);
	for(unsigned int i = range.first; i < range.second; ++ i)
	{
		resids[i - range.first] = mResidueIds[i];
	}
	return resids;
}

}; // namespace annotate
