#ifndef _InteractionInfo_H_
#define _InteractionInfo_H_

#include <string>

#include "ModelInfo.h"

class InteractionInfo
{
public:
	// LIFECYCLE ---------------------------------------------------------------
	InteractionInfo(
		const std::string& astrFile, 
		unsigned int auiModel, 
		const std::string& astrRes1, 
		const std::string& astrRes2);
	~InteractionInfo() {}
	
	// ACCESS ------------------------------------------------------------------
	const ModelInfo& getModelInfo() const {return mModelInfo;}
	const std::string& getPDBFile() const {return mModelInfo.getPDBFile();}
	unsigned int getModel() const {return mModelInfo.getModel();}
	const std::string& getRes1() const {return mstrRes1;}
	const std::string& getRes2() const {return mstrRes2;}
	
private:
	ModelInfo		mModelInfo;
	std::string 	mstrRes1;
	std::string 	mstrRes2;
};

#endif /*_InteractionInfo_H_*/
