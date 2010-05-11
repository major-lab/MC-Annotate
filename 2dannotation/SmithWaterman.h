/*
 * SmithWaterman.h
 *
 *  Created on: Apr 8, 2010
 *      Author: blanchmf
 */

#ifndef _annotate_SmithWaterman_H_
#define _annotate_SmithWaterman_H_

#include <list>
#include <vector>

namespace mccore {
	class Residue;
};

namespace annotate {

class SmithWaterman
{
public:
	typedef std::list<std::pair<int, int> > alignment;

	SmithWaterman(
		const std::vector<mccore::Residue>& aSeq1,
		const std::vector<mccore::Residue>& aSeq2);
	void align();
	const std::vector<mccore::Residue>& sequence1() const {return mSequence1;}
	const std::vector<mccore::Residue>& sequence2() const {return mSequence2;}
	const alignment& bestAlignment() const {return mBestAlignment;}

private:
	struct stCell
	{
		float fScore;
		std::list<std::pair<int, int> > previous;
	};
	typedef std::vector<std::vector<stCell> > dynamic_table;
	std::vector<mccore::Residue> mSequence1;
	std::vector<mccore::Residue> mSequence2;
	dynamic_table mTable; // Dynamic programming table
	std::list<std::pair<int, int> > mBestCells;
	float mfGapScore;
	float mfMatchScore;
	float mfMismatchScore;
	alignment mBestAlignment;

	float fillCell(int aiSeq1, int aiSeq2);
	bool match(const mccore::Residue& aRes1, const mccore::Residue& aRes2) const;
	alignment getSingleBestAlignment() const;
};

}; // namespace annotate

#endif /* _annotate_SmithWaterman_H_ */
