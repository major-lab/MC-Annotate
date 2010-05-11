/*
 * InteractionInfoFile.h
 *
 *  Created on: Mar 10, 2010
 *      Author: blanchmf
 */

#ifndef _annotate_InteractionInfoFile_H_
#define _annotate_InteractionInfoFile_H_

#include "InteractionInfo.h"

// Annotate lib
#include "BasePair.h"

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

	BasePair::enOrientation getOrientation(const std::string& astrOrient) const;
	InteractionInfo::enFace getFace(const std::string& astrGeneralFaces) const;
};

};

#endif /* _annotate_InteractionInfoFile_H_ */

