#include "CycleProfile.h"

#include "StringUtil.h"

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
		std::list<unsigned int>::const_iterator it;
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
};
