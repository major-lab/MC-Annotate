/*
 * RNAGroupFile.h
 *
 *  Created on: Apr 22, 2010
 *      Author: blanchmf
 */

#ifndef _annotate_RNAGroupFile_H_
#define _annotate_RNAGroupFile_H_

#include <string>

namespace annotate
{

class RNAGroupFile
{
public:
	// LIFECYLE ----------------------------------------------------------------
	RNAGroupFile();
	~RNAGroupFile();

	// ACCESSOR ----------------------------------------------------------------

	// METHODS -----------------------------------------------------------------
	void read(const char* aszFilename);

private:
	struct stGroupFileEntry
	{
		std::string strName;
		unsigned int uiModel;
		char cChain;
		int iOffset;
	};
	std::map<unsigned int, std::list<stGroupFileEntry> > mGroups;
};

};

#endif /* _annotate_RNAGroupFile_H_ */
