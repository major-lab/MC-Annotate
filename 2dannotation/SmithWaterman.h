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
	SmithWaterman();
	void align(
		const std::vector<mccore::Residue>& aSeq1,
		const std::vector<mccore::Residue>& aSeq2);
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

	void fillTable();
	float fillCell(int aiSeq1, int aiSeq2);
	bool match(const mccore::Residue& aRes1, const mccore::Residue& aRes2) const;
};

}; // namespace annotate

#endif /* _annotate_SmithWaterman_H_ */
