#ifndef CycleInfo_H_
#define CycleInfo_H_

#include <string>
#include <vector>

#include "InteractionInfo.h"
#include "ModelInfo.h"

class CycleInfo
{
public:
	typedef std::vector<std::string> residue_strand; 
	typedef std::vector<residue_strand> residue_profile;

	// LIFECYLE ----------------------------------------------------------------
	CycleInfo(
		const std::string& aFile, unsigned int auiModel, 
		const std::vector<std::vector<std::string> >& aResidues) 
	: mModelInfo(aFile, auiModel)
	{
		mResidues = aResidues;
	}
	~CycleInfo() {}
	
	// ACCESS ------------------------------------------------------------------
	const ModelInfo& getModelInfo() const {return mModelInfo;}
	const std::string& getPDBFile() const {return mModelInfo.getPDBFile();}
	unsigned int getModel() const {return mModelInfo.getModel();}
	const std::vector<std::vector<std::string> >& getStrandResidues() const 
	{return mResidues;}
	
	// METHOD ------------------------------------------------------------------
	std::vector<unsigned int> getProfile() const;
	bool contains(const InteractionInfo& aInteraction) const;
	std::vector<std::string> getResidues() const;
	
	// OPERATOR ----------------------------------------------------------------
	bool operator <(const CycleInfo& aRight) const;
private:
	ModelInfo mModelInfo;
	residue_profile mResidues;
	
	int compareStrand(
		const residue_strand& aLeft, 
		const residue_strand& aRight) const;
	std::pair<int, int> findResidue(const std::string& astrResidue) const;
};

#endif /*_CycleInfo_H_*/
