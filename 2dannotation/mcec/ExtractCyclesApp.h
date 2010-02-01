/*
 * ExtractCyclesApp.h
 *
 *  Created on: Jan 28, 2010
 *      Author: blanchmf
 */

#ifndef EXTRACTCYCLESAPP_H_
#define EXTRACTCYCLESAPP_H_

#include "CycleInfo.h"
#include <list>
#include <set>
#include <string>
#include <vector>

class ExtractCyclesApp
{
public:
	/**
	 * Constructor.
	 */
	ExtractCyclesApp(int argc, char * argv []);

	// Console specific methods
	void version() const;
	void usage() const;
	void help() const;

	std::string toString() const;

	// METHODS -----------------------------------------------------------------

private:
	std::string mstrCyclesFile;
	std::string mstrInputDirectory;
	std::string mstrOutputDirectory;
	std::set<annotate::CycleInfo> mCycles;

	void readOptions(int argc, char* argv[]);
	void readCyclesFile();

	annotate::CycleInfo readCycleFileLine(const std::string& astrLine) const;

	std::vector<std::vector<std::string> > getStrandResidues(
		const std::string& aResidues,
		const annotate::CycleProfile& aProfile) const;
	std::list<std::string> getResidues(
		const std::string& aResidues) const;
};

#endif /* EXTRACTCYCLESAPP_H_ */
