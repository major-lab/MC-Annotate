/*
 * StringTable.cc
 *
 *  Created on: Aug 27, 2009
 *      Author: blanchmf
 */

#include "StringTable.h"

#include <algorithm>
#include <cassert>
#include <sstream>

namespace annotate
{

	StringTable::StringTable(unsigned int auiNbColumns)
	{
		if(0 < auiNbColumns)
		{
			mColumnWidth.resize(auiNbColumns, 0);
		}
	}

	std::vector<std::string>& StringTable::row(unsigned int auiRow)
	{
		assert(auiRow < mTable.size());
		return mTable[auiRow];
	}

	std::vector<std::string>& StringTable::addRow()
	{
		std::vector<std::string> row;
		row.resize(mColumnWidth.size());
		mTable.push_back(row);
		return mTable.back();
	}

	void StringTable::addRow(std::vector<std::string>& aRow)
	{
		assert(aRow.size() == mColumnWidth.size());
		mTable.push_back(aRow);
	}

	std::string StringTable::toString(const std::string& aSeparator)
	{
		std::ostringstream oss;

		updateColumnWidth();

		for(table::const_iterator itRow = mTable.begin(); itRow != mTable.end(); ++itRow)
		{
			for(std::size_t i = 0; i < mColumnWidth.size(); ++ i)
			{
				std::string strCell = (*itRow)[i];
				strCell.resize(mColumnWidth[i], ' ');

				if(0 < i)
				{
					oss << aSeparator;
				}
				oss << strCell;
			}
			oss << std::endl;
		}
		return oss.str();
	}

	void StringTable::updateColumnWidth()
	{
		mColumnWidth.clear();
		if(0 < mTable.size())
		{
			mColumnWidth.resize(mTable[0].size(), 0);
		}
		for(table::const_iterator itRow = mTable.begin(); itRow != mTable.end(); ++itRow)
		{
			for(unsigned int i = 0; i < mColumnWidth.size(); ++ i)
			{
				mColumnWidth[i] = std::max(mColumnWidth[i], (*itRow)[i].size());
			}
		}
	}
}
