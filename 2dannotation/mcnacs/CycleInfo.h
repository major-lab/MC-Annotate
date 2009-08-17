#ifndef CycleInfo_H_
#define CycleInfo_H_

#include <set>
#include <string>
#include <vector>

#include "InteractionInfo.h"
#include "ModelInfo.h"

class Interaction;

class CycleInfo
{
public:
	typedef std::pair<std::string, std::string> interaction;
	typedef std::vector<std::string> residue_strand; 
	typedef std::vector<residue_strand> residue_profile;

	// LIFECYLE ----------------------------------------------------------------
	CycleInfo(
		const std::string& aFile, unsigned int auiModel, 
		const residue_profile& aResidues) 
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
	unsigned int getNbStrands() const {return mResidues.size();}
	std::set<Interaction> getStrandInteractions(unsigned int auiStrand) const;
	bool shareInteraction(const std::set<Interaction>& aInteraction) const;
	bool isSubCycleOf(const CycleInfo& aCycleInfo) const;
	bool isSubLoopOf(const CycleInfo& aCycleInfo) const;
	bool isSub2StrandsCycle(const CycleInfo& aCycleInfo) const;
	bool is2StrandsSubCycleOfLoop(const CycleInfo& aCycleInfo) const;
	bool isLoopSubCycleOf2Strands(const CycleInfo& aCycleInfo) const;
	
	/**
	 * @brief Verify if the given interactions are all covered by a single 
	 * strand of this cycle.
	 * @return true if there is a perfect match, false otherwise.
	 */
	bool hasStrandCoveringInteractions(
		const std::set<Interaction>& aInteractions) const;
		
	bool strandCoversInteractions(
		unsigned int uiStrand, 
		const std::set<Interaction>& aInteractions) const;
	
	// OPERATOR ----------------------------------------------------------------
	bool operator <(const CycleInfo& aRight) const;
protected:
	ModelInfo mModelInfo;
	residue_profile mResidues;
private:
	
	int compareStrand(
		const residue_strand& aLeft, 
		const residue_strand& aRight) const;
	std::pair<int, int> findResidue(const std::string& astrResidue) const;
};

#endif /*_CycleInfo_H_*/
