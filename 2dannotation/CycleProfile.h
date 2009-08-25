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
		
		// METHODS -------------------------------------------------------------
		std::string toString() const;
	private:
		std::list<unsigned int> mStrandProfile;
		Cycle::enType			meType;
		
		void clear();
	};
	
}; // annotate

#endif /*_annotate_CycleProfile_H_*/
