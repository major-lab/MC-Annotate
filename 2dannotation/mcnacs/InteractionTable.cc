/*
 * InteractionTable.cc
 *
 *  Created on: Jan 18, 2010
 *      Author: blanchmf
 */

#include "InteractionTable.h"

#include "CycleInfo.h"
#include "StringTable.h"

#include <sstream>

InteractionTable::InteractionTable()
{
	// Selection of the profile
	std::list<std::string> profileList = getProfileList();
	std::list<std::string>::const_iterator itProf = profileList.begin();
	unsigned int uiProf = 0;
	for(; itProf != profileList.end(); ++ itProf, ++ uiProf)
	{
		mIndices.insert(std::pair<std::string, unsigned int>(*itProf, uiProf));
	}

	// Create a table to receive
	mTable.resize(mIndices.size());
	std::vector<std::vector<interacting_set> >::iterator it;
	for(it = mTable.begin(); it != mTable.end(); ++ it)
	{
		it->resize(mIndices.size());
	}
}

InteractionTable::~InteractionTable()
{
	mIndices.clear();
	mTable.clear();
}

void InteractionTable::addInteraction(
	const annotate::CycleInfo& aRow,
	const annotate::CycleInfo& aColumn)
{
	std::string strRowProfile = aRow.getProfile().toString();
	std::string strColumnProfile = aColumn.getProfile().toString();

	std::map<std::string, unsigned int>::const_iterator itRow;
	itRow = mIndices.find(strRowProfile);
	std::map<std::string, unsigned int>::const_iterator itColumn;
	itColumn = mIndices.find(strColumnProfile);

	if(itRow != mIndices.end() && itColumn != mIndices.end())
	{
		unsigned int uiRow = itRow->second;
		unsigned int uiColumn = itColumn->second;
		interacting_pair interaction;
		if(aRow < aColumn)
		{
			interaction = interacting_pair(aRow, aColumn);
		}else
		{
			interaction = interacting_pair(aColumn, aRow);
		}
		mTable[uiRow][uiColumn].insert(interaction);
	}
}

const InteractionTable::interacting_set& InteractionTable::operator()(
		const std::string& astrRow,
		const std::string& astrColumn) const
{
	std::map<std::string, unsigned int>::const_iterator itRow;
	itRow = mIndices.find(astrRow);
	std::map<std::string, unsigned int>::const_iterator itColumn;
	itColumn = mIndices.find(astrColumn);

	if(itRow == mIndices.end())
	{
		std::string strMsg = "Cycle profile not found (";
		strMsg += astrRow;
		strMsg += ")";
		throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
	}
	else if(itColumn == mIndices.end())
	{
		std::string strMsg = "Cycle profile not found (";
		strMsg += astrRow;
		strMsg += ")";
		throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
	}
	unsigned int uiRow = itRow->second;
	unsigned int uiColumn = itColumn->second;
	return mTable[uiRow][uiColumn];
}

std::string InteractionTable::toString() const
{
	std::list<std::string> profileList = getProfileList();
	annotate::StringTable stringTable(mIndices.size() + 1);
	unsigned int uiTableIndex = 0;
	std::vector<string>& tableRow = stringTable.addRow();
	std::list<std::string>::const_iterator itProfile = profileList.begin();
	tableRow[0] = "";
	for(itProfile = profileList.begin(), uiTableIndex = 1; itProfile != profileList.end(); ++ itProfile, ++ uiTableIndex)
	{
		tableRow[uiTableIndex] = *itProfile;
	}
	std::vector<std::vector<interacting_set> >::const_iterator itRow;
	for(itRow = mTable.begin(), itProfile = profileList.begin(); itRow != mTable.end(); ++ itProfile, ++itRow)
	{
		std::vector<string>& tableRow2 = stringTable.addRow();
		tableRow2[0] = *itProfile;
		std::vector<interacting_set>::const_iterator itCol;
		for(itCol = itRow->begin(), uiTableIndex = 1; itCol != itRow->end(); ++ itCol, ++ uiTableIndex)
		{
			std::ostringstream oss;
			oss << itCol->size();
			tableRow2[uiTableIndex] = oss.str();
		}
	}
	return stringTable.toString(" : ");
}

std::list<std::string> InteractionTable::getProfileList() const
{
	std::list<std::string> profileList;
	profileList.push_back(std::string("3"));
	profileList.push_back(std::string("4"));
	profileList.push_back(std::string("5"));
	profileList.push_back(std::string("6"));
	profileList.push_back(std::string("2_2"));
	profileList.push_back(std::string("2_3"));
	profileList.push_back(std::string("2_4"));
	profileList.push_back(std::string("2_5"));
	profileList.push_back(std::string("2_6"));
	profileList.push_back(std::string("3_3"));
	profileList.push_back(std::string("3_4"));
	profileList.push_back(std::string("3_5"));
	profileList.push_back(std::string("4_4"));
	return profileList;
}
