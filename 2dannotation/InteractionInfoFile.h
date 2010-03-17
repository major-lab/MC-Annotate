/*
 * InteractionInfoFile.h
 *
 *  Created on: Mar 10, 2010
 *      Author: blanchmf
 */

#ifndef _annotate_InteractionInfoFile_H_
#define _annotate_InteractionInfoFile_H_

#include "InteractionInfo.h"
#include <set>

namespace annotate
{

class InteractionInfoFile
{
public:
	// LIFECYLE ----------------------------------------------------------------
	InteractionInfoFile();
	~InteractionInfoFile();

	// ACCESSOR ----------------------------------------------------------------
	const std::set<InteractionInfo>& interactions() const {return mInteractions;}

	// METHODS -----------------------------------------------------------------
	void read(const std::string& astrFileName);

private:
	std::set<InteractionInfo> mInteractions;
};

};

#endif /* _annotate_InteractionInfoFile_H_ */

