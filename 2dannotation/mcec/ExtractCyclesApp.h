/*
 * ExtractCyclesApp.h
 *
 *  Created on: Jan 28, 2010
 *      Author: blanchmf
 */

#ifndef EXTRACTCYCLESAPP_H_
#define EXTRACTCYCLESAPP_H_

#include "CycleInfo.h"

#include "GetModelRangeFunctor.h"

#include "mccore/GraphModel.h"
#include "mccore/Molecule.h"
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
	~ExtractCyclesApp();

	// Console specific methods
	void version() const;
	void usage() const;
	void help() const;

	std::string toString() const;

	// METHODS -----------------------------------------------------------------
	std::map<annotate::CycleInfo, mccore::GraphModel> extract() const;
	void writePDB(const std::map<annotate::CycleInfo, mccore::GraphModel>& aCycles) const;

private:
	std::list<std::string> mCyclesFiles;
	std::string mstrInputDirectory;
	std::string mstrOutputDirectory;
	std::set<annotate::CycleInfo> mCycles;
	std::map<annotate::CycleInfo, mccore::GraphModel> mCyclesModels;

	void readOptions(int argc, char* argv[]);
	void readCyclesFile();

	std::list<annotate::ModelInfo> getModels() const;

	mccore::Molecule* loadMoleculeFile (const string &filename) const;
	mccore::Molecule::const_iterator getModel(
		mccore::Molecule& aMolecule,
		unsigned int auiModel) const;
	std::map<annotate::CycleInfo, mccore::GraphModel> getCycleModels(
		const mccore::GraphModel& aMolecule,
		GetModelRangeFunctor<annotate::CycleInfo>::const_range& aRange) const;
	mccore::GraphModel getCycleModel(
		const mccore::GraphModel& aModel,
		const annotate::CycleInfo& aCycle) const;
};

#endif /* EXTRACTCYCLESAPP_H_ */
