#include "CycleInfo.h"

#include "Interaction.h"

#include "../AlgorithmExtra.h"

#include <cassert>

std::vector<unsigned int> CycleInfo::getProfile() const
{
	std::vector<unsigned int> profile;
	std::vector<std::vector<std::string> >::const_iterator it;
	for(it = mResidues.begin();
		it != mResidues.end();
		++ it)
	{
		profile.push_back(it->size());
	}
	return profile;
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
		if(mResidues.size() < aRight.mResidues.size())
		{
			bLess = true;
		}
		else if(mResidues.size() == aRight.mResidues.size())
		{
			residue_profile::const_iterator itLeft = mResidues.begin();
			residue_profile::const_iterator itRight = aRight.mResidues.begin();
			for(;
				itLeft != mResidues.end() && itRight != aRight.mResidues.end(); 
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

std::vector<std::string> CycleInfo::getResidues() const
{
	std::vector<std::string> residues;
	residue_profile::const_iterator itStrand;
	for(itStrand = mResidues.begin(); itStrand != mResidues.end(); ++ itStrand)
	{
		residue_strand::const_iterator itRes;
		for(itRes = itStrand->begin(); itRes != itStrand->end(); ++ itRes)
		{
			residues.push_back(*itRes);
		}
	}
	return residues;	
}

std::pair<int, int> CycleInfo::findResidue(const std::string& astrResidue) const
{
	std::pair<int, int> coord(-1, -1);
	for(int iRow = 0; iRow < (int)mResidues.size(); ++ iRow)
	{
		for(int iCol = 0; iCol < (int)mResidues[iRow].size(); ++ iCol)
		{
			if(mResidues[iRow][iCol] == astrResidue)
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
	std::vector<std::string> residues = getResidues();
	if(2 < residues.size())
	{
		unsigned int i;
		for(i = 0; i < residues.size(); ++ i)
		{
			if(residues[i] == aInteraction.getRes1())
			{
				break;
			}
		}
		if(i < residues.size())
		{
			std::string res2 = aInteraction.getRes2();
			if(0 == i)
			{
				if(residues[i + 1] == res2)
				{
					bContains = true;
				}
				else if(residues.back() == res2)
				{
					bContains = true;
				}
			}
			else if(residues.size() - 1 == i)
			{
				if(residues.front() == res2)
				{
					bContains = true;
				}
				else if(residues[i - 1] == res2)
				{
					bContains = true;
				}
			}
			else
			{
				if(residues[i + 1] == res2 || residues[i - 1] == res2)
				{
					bContains = true;	
				}
			}
		}
	}
	return bContains;	
}

std::set<Interaction> CycleInfo::getStrandInteractions(
	unsigned int auiStrand) const
{
	std::set<Interaction> interactions;
	
	assert(auiStrand < mResidues.size());
	
	std::vector<std::string>::const_iterator it;
	std::vector<std::string>::const_iterator itPrev;
	for(it = mResidues[auiStrand].begin(); 
		it != mResidues[auiStrand].end(); 
		++ it)
	{
		if(it != mResidues[auiStrand].begin())
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
		uiStrand < mResidues.size() && !bShares; 
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
	if(1 == mResidues.size())
	{
		// Loop
		bSubCycle = isSubLoopOf(aCycleInfo);
	}else if(2 == mResidues.size())
	{
		// 2 Strand
		bSubCycle = isSub2StrandsCycle(aCycleInfo);
	}
	return bSubCycle;
}

bool CycleInfo::isSubLoopOf(const CycleInfo& aCycleInfo) const
{
	bool bSubLoop = false;
	assert(mResidues.size() == 1);
	assert(aCycleInfo.mResidues.size() == 1);
	
	if(mResidues[0].size() < aCycleInfo.mResidues[0].size())
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
	assert(mResidues.size() == 2);
	assert(aCycleInfo.mResidues.size() == 2);
	
	unsigned int iStrand;
	for(iStrand = 0; iStrand < mResidues.size(); ++ iStrand)
	{
		if(mResidues[iStrand].size() <= aCycleInfo.mResidues[iStrand].size())
		{
			std::set<Interaction> interactionsLeft;
			std::set<Interaction> interactionsRight;
			interactionsLeft = getStrandInteractions(iStrand);
			interactionsRight = aCycleInfo.getStrandInteractions(iStrand);
			
			std::set<Interaction> intersection;
			intersection = annotate::SetIntersection(
				interactionsLeft, 
				interactionsRight);
			bSubLoop = (intersection.size() == interactionsLeft.size());
		}
		else
		{
			break;
		}
	}
	
	// If we got through, than it is a sub cycle
	bSubLoop = (iStrand == mResidues.size());
	
	return bSubLoop;
}

