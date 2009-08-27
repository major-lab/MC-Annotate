/*
 * StringTable.h
 *
 *  Created on: Aug 27, 2009
 *      Author: blanchmf
 */

#ifndef _annotate_StringTable_H_
#define _annotate_StringTable_H_

#include <string>
#include <vector>

namespace annotate
{
	class StringTable
	{
	public:


		StringTable(unsigned int uiNbColumns = 0);

		// ACCESS --------------------------------------------------------------
		std::vector<std::string>& row(unsigned int auiRow);

		// METHODS -------------------------------------------------------------
		std::vector<std::string>& addRow();
		void addRow(std::vector<std::string>& aRow);

		std::string toString(const std::string& aSeparator);

	private:
		typedef std::vector<std::string> table_row;
		typedef std::vector<table_row> table;
		table mTable;
		std::vector< std::size_t > mColumnWidth;

		void updateColumnWidth();
	};
}

#endif /* _annotate_StringTable_H_ */
