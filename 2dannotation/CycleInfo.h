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
	typedef std::vector<std::string> residue_strand;
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
	~CycleInfo() {}

	// ACCESS ------------------------------------------------------------------
	const ModelInfo& getModelInfo() const {return mModelInfo;}
	const std::string& getPDBFile() const {return mModelInfo.getPDBFile();}
	unsigned int getModel() const {return mModelInfo.getModel();}
	const std::vector<std::vector<std::string> >& getStrandResIds() const
	{return mResIds;}
	const annotate::CycleProfile& getProfile() const {return mProfile;}
	const annotate::CycleProfile& getFileProfile() const {return mFileProfile;}

	// METHOD ------------------------------------------------------------------
	const std::vector<std::string> getSequence() const;
	bool contains(const InteractionInfo& aInteraction) const;
	bool contains(const std::string& aResId) const;
	std::vector<std::string> getResIds() const;
	unsigned int getNbStrands() const {return mResIds.size();}
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

	std::string toString(const std::string astrSeparator = ":") const;
	std::string residuesString(const std::string& astrSep = "-") const;
	std::string resIdsString() const;
protected:
	ModelInfo mModelInfo;
	residue_profile mResIds;
private:

	annotate::CycleProfile mProfile; // Profile family
	annotate::CycleProfile mFileProfile;	// Profile in order of the file
	std::map<std::string, std::string> mIdToResMap;

	int compareStrand(
		const residue_strand& aLeft,
		const residue_strand& aRight) const;
	std::pair<int, int> findResId(const std::string& astrResId) const;
};

}; // namespace annotate

#endif /*_CycleInfo_H_*/
