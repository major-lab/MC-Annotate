#ifndef _annotate_CycleProfile_H_
#define _annotate_CycleProfile_H_


#include "Cycle.h"

#include <list>
#include <string>

namespace annotate
{
	class CycleProfile
	{
	public:
		CycleProfile(const std::string& astrProfile);
		~CycleProfile();

		// ACCESS --------------------------------------------------------------
		const std::list<unsigned int>& strandProfile() const
		{
			return mStrandProfile;
		}

		const Cycle::enType& type() const {return meType;}

		// METHODS -------------------------------------------------------------
		std::string toString() const;
		void fromString(const std::string& astrString);
		bool isSymmetric() const;
	private:
		std::list<unsigned int> mStrandProfile;
		Cycle::enType			meType;

		void clear();
	};

}; // annotate

#endif /*_annotate_CycleProfile_H_*/
