#ifndef _annotate_StringUtil_H_
#define _annotate_StringUtil_H_

#include <vector>
#include <string>

namespace annotate
{
	std::vector<std::string> splitStringFields(
		const std::string& aString,
		const std::string& aDelimiter);

	void cleanString(std::string& aString, const char& aChar);
};

#endif /*_annotate_StringUtil_H_*/
