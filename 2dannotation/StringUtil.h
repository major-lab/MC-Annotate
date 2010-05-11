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

	// Same as clean string, but stop at first character not to be deleted
	void cleanStringStart(std::string& aString, const char& aChar);

	// Remove everything from the string after given characters
	std::string cutStringAfter(const std::string& aString, const std::string& astrCutPoint);
};

#endif /*_annotate_StringUtil_H_*/
