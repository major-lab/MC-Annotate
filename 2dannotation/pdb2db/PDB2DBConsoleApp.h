/*
 * PDB2DBConsoleApp.h
 *
 *  Created on: Jul 21, 2010
 *      Author: blanchmf
 */

#ifndef _PDB2DBCONSOLEAPP_H_
#define _PDB2DBCONSOLEAPP_H_

#include "pdb2db.h"

class PDB2DotBracketParamsConsole : public PDB2DotBracketParams
{
public:
	PDB2DotBracketParamsConsole(int argc, char * argv []);
private:
	void read_options (int argc, char* argv[]);

	void version () const;
	void usage () const;
	void help () const;
};

#endif /* _PDB2DBCONSOLEAPP_H_ */
