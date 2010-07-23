#ifndef CycleInfo_H_
#define CycleInfo_H_

#include "InteractionInfo.h"
#include "ModelInfo.h"

// libmcannotate
#include "CycleProfile.h"

#include <set>
#include <string>
#include <vector>

namespace annotate {

class Interaction;

class CycleInfo
{
public:
	typedef std::pair<std::string, std::string> interaction;
	typedef std::vector<mccore::ResId> residue_strand;
	typedef std::vector<residue_strand> residue_profile;

	// LIFECYLE ----------------------------------------------------------------
	CycleInfo();
	CycleInfo(
		const std::string& aFile,
		unsigned int auiModel,
		const annotate::CycleProfile& aFileProfile,
		const annotate::CycleProfile& aProfile,
		const residue_profile& aResIds,
		const std::vector<std::string>& aResidues);
	CycleInfo(
		const std::string& aFile,
		unsigned int auiModel,
		const annotate::CycleProfile& aFileProfile,
		const annotate::CycleProfile& aProfile,
		const std::vector<mccore::ResId>& aResIds,
		const std::vector<std::string>& aResidues);
	~CycleInfo() {}

	// ACCESS ------------------------------------------------------------------
	const ModelInfo& getModelInfo() const {return mModelInfo;}
	const std::string& getPDBFile() const {return mModelInfo.getPDBFile();}
	unsigned int getModel() const {return mModelInfo.getModel();}
	const annotate::CycleProfile& getProfile() const {return mProfile;}
	const annotate::CycleProfile& getFileProfile() const {return mFileProfile;}
	const std::vector<mccore::ResId>& getResIds() const {return mResidueIds;}
	unsigned int getNbStrands() const {return mProfile.strandProfile().size();}

	// METHOD ------------------------------------------------------------------
	const std::vector<std::string> getSequence() const;
	std::string getNucleotideString(const mccore::ResId& aResId) const;
	bool contains(const InteractionInfo& aInteraction) const;
	bool contains(const mccore::ResId& aResId) const;

	std::set<Interaction> getStrandInteractions(unsigned int auiStrand) const;
	bool shareInteraction(const std::set<Interaction>& aInteraction) const;
	bool isSubCycleOf(const CycleInfo& aCycleInfo) const;
	bool isSubLoopOf(const CycleInfo& aCycleInfo) const;
	bool isSub2StrandsCycle(const CycleInfo& aCycleInfo) const;
	bool is2StrandsSubCycleOfLoop(const CycleInfo& aCycleInfo) const;
	bool isLoopSubCycleOf2Strands(const CycleInfo& aCycleInfo) const;


	static CycleInfo flipStrand(const CycleInfo& aCycleInfo);
	const std::vector<std::string> getFlipStrandSequence() const;

	// Get the chains forming the cycle
	std::set<char> getChains() const;
	void setChainAndOffset(char acChain, int aiOffset);

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

	residue_profile getStrandsResIds() const;
	std::vector<mccore::ResId> getStrandResIds(unsigned int auiStrand) const;

	// OPERATOR ----------------------------------------------------------------
	bool operator <(const CycleInfo& aRight) const;
	bool operator ==(const CycleInfo& aRight) const;

	std::string toString(const std::string astrSeparator = ":") const;
	std::string residuesString(const std::string& astrSeparator = "-") const;
	std::string resIdsString(const std::string& astrSeparator = "-") const;
	std::string groupedResIdsString() const;

protected:
	ModelInfo mModelInfo;
	std::vector<mccore::ResId> mResidueIds;
	annotate::CycleProfile mProfile; // Profile family
	annotate::CycleProfile mFileProfile;	// Profile in order of the file
private:
	std::map<mccore::ResId, std::string> mIdToResMap; // TODO : Is this really needed ?

	std::pair<int, int> findResId(const mccore::ResId& astrResId) const;
	std::pair<unsigned int, unsigned int> getStrandRange(unsigned int auiStrand) const;
};

}; // namespace annotate

#endif /*_CycleInfo_H_*/
