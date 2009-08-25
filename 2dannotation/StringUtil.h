#ifndef _annotate_StringUtil_H_
#define _annotate_StringUtil_H_

#include <list>
#include <string>

namespace annotate
{
	std::list<std::string> splitStringFields(
		const std::string& aString, 
		const std::string& aDelimiter);
};

#endif /*_annotate_StringUtil_H_*/
