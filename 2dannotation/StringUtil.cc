#include "StringUtil.h"

namespace annotate
{
	std::vector<std::string> splitStringFields(
		const std::string& aString,
		const std::string& aDelimiter)
	{
		std::string strWork = aString;
		std::vector<std::string> stringFields;

		while(!strWork.empty())
		{
			std::string::size_type index;
			std::string strField;
			if(std::string::npos != (index = strWork.find(aDelimiter)))
		    {
		    	strField = strWork.substr(0, index);
		    	strWork.erase (0, index + aDelimiter.size());
		    }
		    else
		    {
		    	strField = strWork;
		    	strWork.clear();
		    }
		    stringFields.push_back(strField);
		}
		return stringFields;
	}

	void cleanString(std::string& aString, const char& aChar)
	{
		std::string::iterator it = aString.begin();
		while(it != aString.end())
		{
			if(*it == aChar)
			{
				it = aString.erase(it);
			}
			else
			{
				++ it;
			}
		}
	}
};
