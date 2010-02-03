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
	mResIds = aResIds;

	std::vector<std::string>::const_iterator it = aResidues.begin();
	residue_profile::const_iterator itStrand;
	for(itStrand = mResIds.begin(); itStrand != mResIds.end(); ++ itStrand)
	{
		residue_strand::const_iterator itRes;
		for(itRes = itStrand->begin(); itRes != itStrand->end(); ++ itRes)
		{
			mIdToResMap.insert(std::pair<std::string, std::string>(*itRes, *it));
			++ it;
		}
	}
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
		if(mResIds.size() < aRight.mResIds.size())
		{
			bLess = true;
		}
		else if(mResIds.size() == aRight.mResIds.size())
		{
			residue_profile::const_iterator itLeft = mResIds.begin();
			residue_profile::const_iterator itRight = aRight.mResIds.begin();
			for(;
				itLeft != mResIds.end() && itRight != aRight.mResIds.end();
				++ itLeft, ++ itRight)
			{
				int iCompare = compareStrand(*itLeft, *itRight);
				if(iCompare < 0)
				{
					bLess = true;
					break;
				}else if(iCompare > 0)
				{
					break;
				}
			}
		}
	}
	return bLess;
}

int CycleInfo::compareStrand(
	const residue_strand& aLeft,
	const residue_strand& aRight) const
{
	int iCompare = 0;

	if(aLeft.size() == aRight.size())
	{
		residue_strand::const_iterator leftIt = aLeft.begin();
		residue_strand::const_iterator rightIt = aRight.begin();
		for(;
			leftIt != aLeft.end() && rightIt != aRight.end();
			++leftIt, ++rightIt)
		{
			iCompare = leftIt->compare(*rightIt);
			if(iCompare != 0)
			{
				break;
			}
		}
	}
	else if(aLeft.size() < aRight.size())
	{
		iCompare = -1;
	}
	else
	{
		iCompare = 1;
	}
	return iCompare;
}

std::vector<std::string> CycleInfo::getResIds() const
{
	std::vector<std::string> residues;
	residue_profile::const_iterator itStrand;
	for(itStrand = mResIds.begin(); itStrand != mResIds.end(); ++ itStrand)
	{
		residue_strand::const_iterator itRes;
		for(itRes = itStrand->begin(); itRes != itStrand->end(); ++ itRes)
		{
			residues.push_back(*itRes);
		}
	}
	return residues;
}

std::pair<int, int> CycleInfo::findResId(const std::string& astrResId) const
{
	std::pair<int, int> coord(-1, -1);
	for(int iRow = 0; iRow < (int)mResIds.size(); ++ iRow)
	{
		for(int iCol = 0; iCol < (int)mResIds[iRow].size(); ++ iCol)
		{
			if(mResIds[iRow][iCol] == astrResId)
			{
				coord.first = iRow;
				coord.second = iCol;
			}
		}
	}
	return coord;
}

bool CycleInfo::contains(const InteractionInfo& aInteraction) const
{
	bool bContains = false;
	std::vector<std::string> resIds = getResIds();
	if(2 < resIds.size())
	{
		unsigned int i;
		for(i = 0; i < resIds.size(); ++ i)
		{
			if(resIds[i] == aInteraction.getRes1())
			{
				break;
			}
		}
		if(i < resIds.size())
		{
			std::string res2 = aInteraction.getRes2();
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

bool CycleInfo::contains(const std::string& aResId) const
{
	bool bContains = false;
	std::vector<std::string> resIds = getResIds();
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

	assert(auiStrand < mResIds.size());

	std::vector<std::string>::const_iterator it;
	std::vector<std::string>::const_iterator itPrev;
	for(it = mResIds[auiStrand].begin();
		it != mResIds[auiStrand].end();
		++ it)
	{
		if(it != mResIds[auiStrand].begin())
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
		uiStrand < mResIds.size() && !bShares;
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

	if(1 == mResIds.size() && 1 == aCycleInfo.mResIds.size())
	{
		// Loop
		bSubCycle = isSubLoopOf(aCycleInfo);
	}else if(2 == mResIds.size() && 2 == aCycleInfo.mResIds.size())
	{
		// 2 Strand
		bSubCycle = isSub2StrandsCycle(aCycleInfo);
	}
	else
	{
		assert(1 == mResIds.size() || 2 == mResIds.size());
		assert(1 == aCycleInfo.mResIds.size() || 2 == aCycleInfo.mResIds.size());

		if(2 == mResIds.size())
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
	assert(mResIds.size() == 1);
	assert(aCycleInfo.mResIds.size() == 1);

	if(mResIds[0].size() < aCycleInfo.mResIds[0].size())
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
	assert(mResIds.size() == 2);
	assert(aCycleInfo.mResIds.size() == 2);

	unsigned int iStrand;
	for(iStrand = 0; iStrand < mResIds.size(); ++ iStrand)
	{
		if(mResIds[iStrand].size() <= aCycleInfo.mResIds[iStrand].size())
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
	bSubLoop = (iStrand == mResIds.size());

	return bSubLoop;
}

bool CycleInfo::is2StrandsSubCycleOfLoop(const CycleInfo& aCycleInfo) const
{
	bool bSubLoop = false;
	assert(mResIds.size() == 2);
	assert(aCycleInfo.mResIds.size() == 1);


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
	assert(mResIds.size() == 1);
	assert(aCycleInfo.mResIds.size() == 2);

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
	for(uiStrand = 0; uiStrand < mResIds.size() && !bMatch; ++ uiStrand)
	{
		bMatch = strandCoversInteractions(uiStrand, aInteractions);
	}
	return bMatch;
}

const std::vector<std::string> CycleInfo::getSequence() const
{
	std::vector<std::string> sequence;
	residue_profile::const_iterator itStrand;
	for(itStrand = mResIds.begin(); itStrand != mResIds.end(); ++ itStrand)
	{
		residue_strand::const_iterator itRes;
		for(itRes = itStrand->begin(); itRes != itStrand->end(); ++ itRes)
		{
			std::map<std::string, std::string>::const_iterator it;
			it = mIdToResMap.find(*itRes);
			assert(it != mIdToResMap.end());
			sequence.push_back(it->second);
		}
	}
	return sequence;
}

std::string CycleInfo::toString(const std::string astrSeparator) const
{
	std::ostringstream oss;

	oss << mProfile.toString();
	oss << "\t" << astrSeparator << " ";
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

	CycleInfo::residue_profile::const_iterator itResStrand;
	for(itResStrand = mResIds.begin();
		itResStrand != mResIds.end();
		++ itResStrand)
	{
		CycleInfo::residue_strand::const_iterator itRes;
		for(itRes = itResStrand->begin();
			itRes != itResStrand->end();
			++ itRes)
		{
			if(!(itResStrand == mResIds.begin() && itRes == itResStrand->begin()))
			{
				oss << "-";
			}
			oss << *itRes;
		}
	}

	return oss.str();
}

std::string CycleInfo::groupedResIdsString() const
{
	std::ostringstream oss;

	CycleInfo::residue_profile::const_iterator itResStrand;
	for(itResStrand = mResIds.begin();
		itResStrand != mResIds.end();
		++ itResStrand)
	{
		if(itResStrand != mResIds.begin())
		{
			oss << ",";
		}
		oss << "{";
		CycleInfo::residue_strand::const_iterator itRes;
		for(itRes = itResStrand->begin();
			itRes != itResStrand->end();
			++ itRes)
		{
			if(itRes != itResStrand->begin())
			{
				oss << "-";
			}
			oss << *itRes;
		}
		oss << "}";
	}

	return oss.str();
}

}; // namespace annotate
