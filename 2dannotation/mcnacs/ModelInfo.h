#ifndef _ModelInfo_H_
#define _ModelInfo_H_

#include <string>

class ModelInfo
{
public:
	// LIFECYCLE ---------------------------------------------------------------
	ModelInfo(const std::string& astrFile, unsigned int auiModel);
	~ModelInfo() {}
	
	// ACCESS ------------------------------------------------------------------
	const std::string& getPDBFile() const {return mstrPDBFile;}
	unsigned int getModel() const {return muiModel;}
	
	// OPERATOR ----------------------------------------------------------------
	bool operator <(const ModelInfo& aRight) const;
	bool operator ==(const ModelInfo& aRight) const;
	bool operator !=(const ModelInfo& aRight) const;
private:
	std::string mstrPDBFile;
	unsigned int muiModel;
};

#endif /*_ModelInfo_H_*/
