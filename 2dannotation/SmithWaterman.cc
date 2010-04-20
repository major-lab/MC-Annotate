/*
 * SmithWaterman.cc
 *
 *  Created on: Apr 8, 2010
 *      Author: blanchmf
 */

#include "SmithWaterman.h"

#include <mccore/Residue.h>

#include <cassert>

namespace annotate {

SmithWaterman::SmithWaterman()
{
	mfGapScore = -4.0f;
	mfMatchScore = 10.0f;
	mfMismatchScore = -7.0f;
}

void SmithWaterman::align(
	const std::vector<mccore::Residue>& aSeq1,
	const std::vector<mccore::Residue>& aSeq2)
{
	mSequence1 = aSeq1;
	mSequence2 = aSeq2;

	mTable.resize(aSeq2.size() + 1);
	dynamic_table::iterator it;
	for(it = mTable.begin(); it != mTable.end(); ++ it)
	{
		it->resize(aSeq1.size() + 1);
		std::vector<stCell>::iterator it2;
		for(it2 = it->begin(); it2 != it->end(); ++ it2)
		{
			it2->fScore = 0.0f;
		}
	}

	// Fill the table
	fillTable();
}

void SmithWaterman::fillTable()
{
	float fBestScore = 0.0f;
	int iSeq1;
	int iSeq2;
	int iSeq1Size = (int)mSequence1.size();
	int iSeq2Size = (int)mSequence2.size();

	for(iSeq1 = 0; iSeq1 < iSeq1Size; ++ iSeq1)
	{
		for(iSeq2 = 0; iSeq2 < iSeq2Size; ++ iSeq2)
		{
			float fScore = fillCell(iSeq1, iSeq2);
			if(fBestScore < fScore)
			{
				fBestScore = fScore;
				mBestCells.clear();
			}
			if(fBestScore == fScore)
			{
				mBestCells.push_back(std::pair<int, int>(iSeq1, iSeq2));
			}
		}
	}
}

float SmithWaterman::fillCell(int aiSeq1, int aiSeq2)
{
	assert(0 < aiSeq1);
	assert(0 < aiSeq2);

	float fBestScore = 0.0f;
	float fMatchScore = 0.0f;

	// Gap in sequence 1
	float fGapSeq1Score = 0.0f;
	if(0 < aiSeq1)
	{
		fGapSeq1Score = mTable[aiSeq1 - 1][aiSeq2].fScore + mfGapScore;
	}

	// Gap in sequence 2
	float fGapSeq2Score = 0.0f;
	if(0 < aiSeq2)
	{
		fGapSeq2Score = mTable[aiSeq1][aiSeq2 - 1].fScore + mfGapScore;
	}

	// Match or mismatch
	if(match(mSequence1[aiSeq1], mSequence2[aiSeq2]))
	{
		fMatchScore = mfMatchScore;
	}else
	{
		fMatchScore += mfMismatchScore;
	}

	fBestScore = std::max(fGapSeq1Score, std::max(fGapSeq2Score, fMatchScore));
	if(0.0f < fBestScore)
	{
		if(fBestScore == fGapSeq1Score)
		{
			mTable[aiSeq1][aiSeq2].previous.push_back(std::pair<int, int>(aiSeq1 - 1, aiSeq2));
		}
		if(fBestScore == fGapSeq2Score)
		{
			mTable[aiSeq1][aiSeq2].previous.push_back(std::pair<int, int>(aiSeq1, aiSeq2 - 1));
		}
		if(fBestScore == fMatchScore)
		{
			mTable[aiSeq1][aiSeq2].previous.push_back(std::pair<int, int>(aiSeq1 - 1, aiSeq2 - 1));
		}
	}
	else
	{
		fBestScore = 0.0f;
	}
	mTable[aiSeq1][aiSeq2].fScore = fBestScore;
	return fBestScore;
}

bool SmithWaterman::match(
	const mccore::Residue& aRes1,
	const mccore::Residue& aRes2) const
{
	return aRes1.getType() == aRes2.getType();
}

}
