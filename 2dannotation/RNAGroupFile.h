/*
 * RNAGroupFile.h
 *
 *  Created on: Apr 22, 2010
 *      Author: blanchmf
 */

#ifndef _annotate_RNAGroupFile_H_
#define _annotate_RNAGroupFile_H_

#include <list>
#include <map>
#include <string>

namespace annotate
{

// PROTOTYPES ------------------------------------------------------------------
class CycleInfo;
class InteractionInfo;

class RNAGroupFileEntry
{
public:
	// LIFECYLE ----------------------------------------------------------------
	RNAGroupFileEntry(
		const std::string astrPDB,
		unsigned int auiModel,
		char acChain,
		int aiOffset)
	{
		mstrName = astrPDB;
		muiModel = auiModel;
		mcChain = acChain;
		miOffset = aiOffset;
	}

	// ACCESSOR ----------------------------------------------------------------
	std::string name() const {return mstrName;}
	unsigned int model() const {return muiModel;}
	char chain() const {return mcChain;}
	int& offset() {return miOffset;}
	int offset() const {return miOffset;}

	// METHODS -----------------------------------------------------------------
	bool operator<(const RNAGroupFileEntry& aEntry) const
	{
		bool bLess = true;
		if(mstrName > aEntry.mstrName)
		{
			bLess = false;
		}else if(mstrName == aEntry.mstrName)
		{
			if(muiModel > aEntry.muiModel)
			{
				bLess = false;
			}else
			{
				bLess = (mcChain < aEntry.mcChain);
			}
		}
		return bLess;
	}
private:
	std::string mstrName;
	unsigned int muiModel;
	char mcChain;
	int miOffset;
};

class RNAGroupFile
{
public:
	// LIFECYLE ----------------------------------------------------------------
	RNAGroupFile();
	~RNAGroupFile();

	// ACCESSOR ----------------------------------------------------------------
	const std::map<unsigned int, std::list<RNAGroupFileEntry> >& groups() const
	{return mGroups;}

	// METHODS -----------------------------------------------------------------
	void read(const std::string& astrFilename);
	unsigned int getGroup(const RNAGroupFileEntry& aEntry) const;
	unsigned int getGroup(const CycleInfo& aCycle) const;
	unsigned int getGroup(const InteractionInfo& aCycle) const;

private:
	std::map<unsigned int, std::list<RNAGroupFileEntry> > mGroups;
	std::map<RNAGroupFileEntry, unsigned int> mEntryGroupMap;


	void readLine(const std::string& astrLine);
	void renumberOffsets();
};

};

#endif /* _annotate_RNAGroupFile_H_ */
