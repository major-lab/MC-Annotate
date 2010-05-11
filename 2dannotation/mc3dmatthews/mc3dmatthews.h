/*
 * MC3DMatthews.h
 *
 *  Created on: May 3, 2010
 *      Author: blanchmf
 */

#ifndef _mc3dmatthews_H_
#define _mc3dmatthews_H_

#include <string>
#include <map>
#include <vector>

#include "AlgorithmExtra.h"

class MC3DMatthews
{
public:
	MC3DMatthews(int argc, char * argv []);

	// METHODS ---------------------------------------------
	void version () const;
	void usage () const;
	void help () const;
private:
	typedef std::pair<std::string, std::string> cycle_pair;
	typedef std::pair<unsigned int, unsigned int> inter_coords;
	typedef std::string inter_faces;
	typedef annotate::tuple3<cycle_pair, inter_faces, inter_coords> cycle_interaction;
	float mfAppVersion;
	std::string mstrAppName;

	std::string mstrHypotheticalInteractionsFile;
	std::multimap<float, cycle_interaction> mHypothesis;

	void read_options (int argc, char* argv[]);

	void readHypotheticalInteractionsFile();

	void displayHypothesis() const;
};
#endif // _mc3dmatthews_H_
