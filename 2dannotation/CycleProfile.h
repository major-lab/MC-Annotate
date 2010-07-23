#ifndef _annotate_CycleProfile_H_
#define _annotate_CycleProfile_H_


#include "Cycle.h"

#include <string>
#include <vector>

namespace annotate
{
	class CycleProfile
	{
	public:
		CycleProfile(const std::string& astrProfile);
		~CycleProfile();

		// ACCESS --------------------------------------------------------------
		const std::vector<unsigned int>& strandProfile() const
		{
			return mStrandProfile;
		}
		const Cycle::enType& type() const {return meType;}

		// OPERATORS -----------------------------------------------------------
		bool operator <(const CycleProfile& aRight) const;
		bool operator ==(const CycleProfile& aRight) const;

		// METHODS -------------------------------------------------------------
		std::string toString() const;
		void fromString(const std::string& astrString);
		bool isSymmetric() const;
		static CycleProfile Rotate(const CycleProfile& aProfile);
	private:
		std::vector<unsigned int> mStrandProfile;
		Cycle::enType			meType;

		void clear();
	};

}; // annotate

#endif /*_annotate_CycleProfile_H_*/
