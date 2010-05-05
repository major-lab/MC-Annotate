/*
 * mc3dihg.cc
 *
 *  Created on: May 3, 2010
 *      Author: blanchmf
 */

#include "mc3dihf.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>

#include "StringTable.h"
#include "StringUtil.h"

static const char* gszShortopts = "Vhi:";

MC3DInteractionHypothesisFilter::MC3DInteractionHypothesisFilter(int argc, char * argv [])
: mfAppVersion(1.0f),
mstrAppName("MC3DInteractionHypothesisFilter")
{
	read_options (argc, argv);

	readHypotheticalInteractionsFile();

	displayHypothesis();
}


void MC3DInteractionHypothesisFilter::version () const
{
	std::cout << mstrAppName << " "; // Nom du logiciel
	std::cout << mfAppVersion << " ";			// Version du logiciel
	std::cout << "(" << __DATE__ << ")";	// Date de la compilation
	std::cout << std::endl;
}

void MC3DInteractionHypothesisFilter::usage () const
{
	std::cout << "usage: "
		<< " [-hV] -i <hypothetical interaction file>"
		<< std::endl;
}


void MC3DInteractionHypothesisFilter::help () const
{
	std::cout
		<< "This program filters the hypothetical interactions." << std::endl
		<< "  -h	print this information" << std::endl
		<< "  -V    print the software version info" << std::endl
		<< "  -i	File containing the scored hypothetical interactions." << std::endl;
}

void MC3DInteractionHypothesisFilter::read_options (int argc, char* argv[])
{
	int c;

	while ((c = getopt (argc, argv, gszShortopts)) != -1)
	{
		switch (c)
		{
		case 'V':
			version ();
			exit (EXIT_SUCCESS);
			break;
		case 'h':
			help();
			exit(EXIT_SUCCESS);
			break;
		case 'i':
			mstrHypotheticalInteractionsFile = optarg;
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if (0 == mstrHypotheticalInteractionsFile.size())
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

void MC3DInteractionHypothesisFilter::readHypotheticalInteractionsFile()
{
	std::map<cycle_pair, std::pair<cycle_interaction, float> > hypothesis;
	std::map<cycle_pair, std::pair<cycle_interaction, float> >::iterator it;
	std::ifstream infile;
	infile.open(mstrHypotheticalInteractionsFile.c_str(), std::ios_base::in);
	if(infile.good())
	{
		std::string strLine;
		while(std::getline(infile, strLine).good())
		{
			// Remove comments
			strLine = annotate::cutStringAfter(strLine, "//");
			if(0 < strLine.size())
			{
				// Remove whitespaces
				annotate::cleanString(strLine, ' ');

				// Read the data
				std::vector<std::string> fields = annotate::splitStringFields(strLine, ":");
				float fScore = atof(fields[0].c_str());
				std::string strFaces = fields[1];
				unsigned int uiPos1 = atoi(fields[2].c_str());
				unsigned int uiPos2 = atoi(fields[3].c_str());
				cycle_pair interaction;
				interaction.first = fields[4];
				interaction.second = fields[5];

				it = hypothesis.find(interaction);
				if(it != hypothesis.end())
				{
					if(it->second.second < fScore)
					{
						cycle_interaction inter(interaction, strFaces, inter_coords(uiPos1, uiPos2));
						it->second = std::pair<cycle_interaction, float>(inter, fScore);
					}
				}
				else
				{
					cycle_interaction inter(interaction, strFaces, inter_coords(uiPos1, uiPos2));
					std::pair<cycle_interaction, float> entry(inter, fScore);
					hypothesis.insert(std::pair<cycle_pair, std::pair<cycle_interaction, float> >(interaction, entry));
				}
			}
		}
	}
	else
	{
		// TODO : This should be an exception
		std::cout << "Error opening file " << mstrHypotheticalInteractionsFile << std::endl;
	}
	infile.close();

	for(it = hypothesis.begin(); it != hypothesis.end(); ++ it)
	{
		mHypothesis.insert(std::pair<float, cycle_interaction>(it->second.second, it->second.first));
	}
}

void MC3DInteractionHypothesisFilter::displayHypothesis() const
{
	std::map<float, cycle_interaction>::const_iterator it;
	for(it = mHypothesis.begin(); it != mHypothesis.end(); ++ it)
	{
		std::cout << it->first;
		std::cout << " : " << it->second.second;
		std::cout << " : " << it->second.third.first;
		std::cout << " : " << it->second.third.second;
		std::cout << " : " << it->second.first.first;
		std::cout << " : " << it->second.first.second;
		std::cout << std::endl;
	}
}

int main(int argc, char* argv[])
{
	MC3DInteractionHypothesisFilter theApp(argc, argv);
	return EXIT_SUCCESS;
}
