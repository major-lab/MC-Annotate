#include "CycleProfile.h"

#include "StringUtil.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <sstream>

namespace annotate
{
	CycleProfile::CycleProfile(const std::string& astrString)
	{
		this->fromString(astrString);
	}

	CycleProfile::~CycleProfile()
	{
		clear();
	}

	void CycleProfile::clear()
	{
		mStrandProfile.clear();
	}

	std::string CycleProfile::toString() const
	{
		std::ostringstream oss;
		std::vector<unsigned int>::const_iterator it;
		for(it = mStrandProfile.begin(); it != mStrandProfile.end(); ++ it)
		{
			if(it != mStrandProfile.begin())
			{
				oss << "_";
			}
			oss << *it;
		}

		switch(meType)
		{
		case Cycle::eLOOSE:
			oss << "_L";
			break;
		case Cycle::e2STRANDS_PARALLEL:
			oss << "_p";
			break;
		case Cycle::eMULTIBRANCH:
			oss << "_M";
			break;
		default:
			break;
		}

		return oss.str();
	}

	void CycleProfile::fromString(const std::string& astrString)
	{
		// Reset the CycleProfile
		this->clear();

		// Read the cycle profile from the string
		std::vector<std::string> fields = splitStringFields(astrString, "_");
		std::vector<std::string>::const_iterator it;
		for(it = fields.begin(); it != fields.end(); ++ it)
		{
			if(std::isdigit((*it)[it->size() - 1]))
			{
				unsigned int uiSize = std::atol(it->c_str());
				mStrandProfile.push_back(uiSize);
			}
		}

		if(1 == mStrandProfile.size())
		{
			if(fields.back() == "L")
			{
				meType = Cycle::eLOOSE;
			}
			else if(fields.back() == "M")
			{
				meType = Cycle::eMULTIBRANCH;
			}
			else
			{
				meType = Cycle::eLOOP;
			}
		}
		else if(2 == mStrandProfile.size())
		{
			if(fields.back() == "p")
			{
				meType = Cycle::e2STRANDS_PARALLEL;
			}
			else
			{
				meType = Cycle::e2STRANDS_ANTIPARALLEL;
			}
		}
		else
		{
			assert(false);
		}
	}

	// A profile is symmetric if all the strands have the same length
	bool CycleProfile::isSymmetric() const
	{
		std::vector<unsigned int>::const_iterator it;
		bool bSymmetric = true;
		unsigned int uiLength = 0;

		for(it = mStrandProfile.begin(); it != mStrandProfile.end() && bSymmetric; ++ it)
		{
			if(it == mStrandProfile.begin())
			{
				uiLength = *it;
			}else if(uiLength != *it)
			{
				bSymmetric = false;
			}
		}
		return bSymmetric;
	}

	CycleProfile CycleProfile::Rotate(const CycleProfile& aProfile)
	{
		CycleProfile profile = aProfile;
		std::vector<unsigned int>::iterator it = profile.mStrandProfile.begin();
		it ++;
		std::rotate(profile.mStrandProfile.begin(),it, profile.mStrandProfile.end());
		return profile;
	}

	bool CycleProfile::operator <(const CycleProfile& aRight) const
	{
		bool bLess = false;

		if(mStrandProfile.size() < aRight.mStrandProfile.size())
		{
			bLess = true;
		}
		else if(mStrandProfile.size() == aRight.mStrandProfile.size())
		{
			std::vector<unsigned int>::const_iterator it1;
			std::vector<unsigned int>::const_iterator it2;
			for(it1 = mStrandProfile.begin(), it2 = aRight.mStrandProfile.begin();
				it1 != mStrandProfile.end();
				++ it1, ++ it2)
			{
				if(*it1 < *it2)
				{
					bLess = true;
					break;
				}else if(*it2 < *it1)
				{
					break;
				}
			}
		}
		return bLess;
	}

	bool CycleProfile::operator ==(const CycleProfile& aRight) const
	{
		bool bEqual = true;

		if(mStrandProfile.size() == aRight.mStrandProfile.size())
		{
			std::vector<unsigned int>::const_iterator it1;
			std::vector<unsigned int>::const_iterator it2;
			for(it1 = mStrandProfile.begin(), it2 = aRight.mStrandProfile.begin();
				it1 != mStrandProfile.end();
				++ it1, ++ it2)
			{
				if(*it1 != *it2)
				{
					bEqual = false;
					break;
				}
			}
		}
		else
		{
			bEqual = false;
		}
		return bEqual;
	}
};
