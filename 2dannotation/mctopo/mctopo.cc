//                              -*- Mode: C++ -*-
// mctopo.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Jan 18 15:03:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <iterator>
#include <string>
#include <sstream>
#include <unistd.h>
#include <vector>

#include "mccore/Binstream.h"
#include "mccore/Exception.h"
#include "mccore/Messagestream.h"
#include "mccore/ModelFactoryMethod.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"
#include "mccore/ResidueFactoryMethod.h"
#include "mccore/ResIdSet.h"
#ifdef HAVE_LIBRNAMLC__
#include "mccore/RnamlReader.h"
#endif
#include "mccore/Version.h"

#include "AnnotateModel.h"
#include "AnnotationChains.h"
#include "AnnotationInteractions.h"
#include "BaseStack.h"

bool binary = false;
unsigned int environment = 0;
bool oneModel = false;
unsigned int modelNumber = 0;  // 1 based vector identifier, 0 means all
unsigned char gucRelationMask =
	mccore::Relation::adjacent_mask
	| mccore::Relation::pairing_mask
	| mccore::Relation::stacking_mask
	| mccore::Relation::backbone_mask;
ResIdSet residueSelection;
const char* shortopts = "Vbe:f:hlr:vm:";

void
version ()
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void
usage ()
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-bhlvV] [-e num] [-f <model number>] [-r <residue ids>] [-m <relation mask>] <structure file> ..."
		<< std::endl;
}

void
help ()
{
	mccore::gOut (0)
		<< "This program annotate structures (and more)." << std::endl
		<< "  -b                read binary files instead of pdb files" << std::endl
		<< "  -m <mask>         annotation mask string: any combination of 'A' (adjacent), 'S' (stacking), 'P' (pairing) and 'B' (backbone). (default: all)" << std::endl
		<< "  -f model number   model to print" << std::endl
		<< "  -h                print this help" << std::endl
		<< "  -l                be more verbose (log)" << std::endl
		<< "  -v                be verbose" << std::endl
		<< "  -V                print the software version info" << std::endl;
}


void
read_options (int argc, char* argv[])
{
	int c;

	while ((c = getopt (argc, argv, shortopts)) != EOF)
	{
		switch (c)
		{
		case 'V':
			version ();
			exit (EXIT_SUCCESS);
			break;
		case 'b':
			binary = true;
			break;
		case 'f':
		{
			long int tmp;

			tmp = strtol (optarg, 0, 10);
			if (ERANGE == errno	|| EINVAL == errno || 0 > tmp)
			{
				mccore::gErr(0) << PACKAGE << ": invalid model value.";
				mccore::gErr(0) << std::endl;
				exit (EXIT_FAILURE);
			}
			modelNumber = tmp;
			oneModel = true;
			break;
		}
		case 'h':
			usage ();
			help ();
			exit (EXIT_SUCCESS);
			break;
		case 'l':
			mccore::gErr.setVerboseLevel (mccore::gErr.getVerboseLevel () + 1);
			break;
		case 'm':
    	{
    		std::string strMask = optarg;
    		gucRelationMask = 0;
    		if(strMask.find("A") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::adjacent_mask;
    		}
    		if(strMask.find("S") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::stacking_mask;
    		}
    		if(strMask.find("P") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::pairing_mask;
    		}
    		if(strMask.find("B") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::backbone_mask;
    		}
			break;
    	}
		case 'v':
			mccore::gOut.setVerboseLevel (mccore::gOut.getVerboseLevel () + 1);
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if (argc - optind < 1)
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}


mccore::Molecule* loadFile (const string &filename)
{
  Molecule *molecule;
  ResidueFM rFM;
  annotate::AnnotateModelFM aFM (residueSelection, environment, &rFM);

  molecule = 0;
  if (binary)
    {
      izfBinstream in;

      in.open (filename.c_str ());
      if (in.fail ())
	{
	  gErr (0) << PACKAGE << ": cannot open binary file '" << filename << "'." << endl;
	  return 0;
	}
      molecule = new Molecule (&aFM);
      in >> *molecule;
      in.close ();
    }
  else
    {
#ifdef HAVE_LIBRNAMLC__
      RnamlReader reader (filename.c_str (), &aFM);

      if (0 == (molecule = reader.read ()))
	{
#endif
	  izfPdbstream in;

	  in.open (filename.c_str ());
	  if (in.fail ())
	    {
	      gErr (0) << PACKAGE << ": cannot open pdb file '" << filename << "'." << endl;
	      return 0;
	    }
	  molecule = new Molecule (&aFM);
	  in >> *molecule;
	  in.close ();
#ifdef HAVE_LIBRNAMLC__
	}
#endif
    }
  return molecule;
}

std::string stripFilePath(const std::string& aFileName)
{
	std::string::size_type index;
	std::string filename = aFileName;
	if (std::string::npos != (index = filename.rfind ("/")))
	{
		filename.erase (0, index + 1);
	}
	return filename;
}

std::string getFilePrefix(const std::string& aFileName)
{
	std::string::size_type index;
	std::string filename = aFileName;
	if (std::string::npos != (index = filename.rfind ("/")))
    {
		filename.erase (0, index + 1);
    }
	if (string::npos != (index = filename.find (".")))
    {
		filename.erase (index, filename.size ());
    }
	return filename;
}

std::string stackingToString(const annotate::AnnotateModel &aModel)
{
	std::ostringstream oss;
	const annotate::AnnotationInteractions* pInteractions = NULL;
	pInteractions = aModel.getAnnotation<annotate::AnnotationInteractions>();
	assert(NULL != pInteractions);

	std::vector< annotate::BaseStack > nonAdjacentStacks;
	std::vector< annotate::BaseStack >::const_iterator bsit;

	for (bsit = pInteractions->stacks().begin (); pInteractions->stacks().end () != bsit; ++bsit)
	{
		const mccore::Relation *rel = aModel.internalGetEdge (bsit->first, bsit->second);
		if (rel->is (PropertyType::pAdjacent))
		{
			const std::set< const mccore::PropertyType* > &labels = rel->getLabels ();
			oss << bsit->fResId << "-" << bsit->rResId << " : " << "(adjacent stacking) ";
			std::copy (labels.begin (), labels.end (), std::ostream_iterator< const mccore::PropertyType* > (oss, " "));
			oss << endl;
		}
		else
		{
			nonAdjacentStacks.push_back (*bsit);
		}
	}

	for (bsit = nonAdjacentStacks.begin (); nonAdjacentStacks.end () != bsit; ++bsit)
	{
		const std::set< const mccore::PropertyType* > &labels = aModel.internalGetEdge (bsit->first, bsit->second)->getLabels ();

		oss << bsit->fResId << "-" << bsit->rResId << " : " << "(non-adj. stacking) ";
		std::copy (labels.begin (), labels.end (), std::ostream_iterator< const mccore::PropertyType* > (oss, " "));
		oss << endl;
	}
	return oss.str();
}

std::string conformationToString(const annotate::AnnotateModel &aModel)
{
	std::ostringstream oss;

	mccore::GraphModel::const_iterator it;

	for (it = aModel.begin (); it != aModel.end (); ++it)
	{
		oss << it->getResId () << " : ";
		if (it->getType ()->isNucleicAcid ())
		{
			oss << " " << it->getPucker();
			oss << " " << it->getGlycosyl();
		}
		oss << std::endl;
	}
	return oss.str();
}

std::string pairingToString(const annotate::AnnotateModel &aModel)
{
	std::ostringstream oss;
	const annotate::AnnotationInteractions* pInteractions = NULL;
	pInteractions = aModel.getAnnotation<annotate::AnnotationInteractions>();
	assert(NULL != pInteractions);

	std::vector< annotate::BasePair >::const_iterator bpit;
	for (bpit = pInteractions->pairs().begin (); pInteractions->pairs().end () != bpit; ++bpit)
	{
		// outputPair(oss, *bpit);
		const mccore::Relation &rel = *aModel.internalGetEdge (bpit->first, bpit->second);
		const std::set< const mccore::PropertyType* > &labels = rel.getLabels ();
		const std::vector< std::pair< const mccore::PropertyType*, const mccore::PropertyType* > > &faces = rel.getPairedFaces ();
		std::vector< pair< const mccore::PropertyType*, const mccore::PropertyType* > >::const_iterator pfit;

		oss << bpit->fResId << '-' << bpit->rResId << " : ";
		for (pfit = faces.begin (); faces.end () != pfit; ++pfit)
		{
			oss << *pfit->first << "/" << *pfit->second << ' ';
		}
		copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (oss, " "));
		oss << endl;
	}


	return oss.str();
}

std::string topology(const annotate::AnnotateModel &aModel)
{
	std::ostringstream oss;
	oss << conformationToString(aModel);
	oss << stackingToString(aModel);
	oss << pairingToString(aModel);
	return oss.str();
}

std::string symmetricTopology(const annotate::AnnotateModel &aModel)
{
	annotate::AnnotateModel model = aModel;
	annotate::AnnotateModel::iterator it;
	for(it = model.begin(); it != model.end(); ++ it)
	{
		mccore::ResId resId = it->getResId();
		if(resId.getChainId() == 'A')
		{
			resId.setChainId('B');
		}else
		{
			resId.setChainId('A');
		}
		it->setResId(resId);
	}
	return topology(model);
}

int main (int argc, char *argv[])
{
	std::map<std::string, std::list<std::string> > topologies;
	read_options (argc, argv);

	while (optind < argc)
    {
		Molecule *molecule;
		Molecule::iterator molIt;
		std::string filename = (std::string) argv[optind];

		molecule = loadFile (filename);
		if (0 != molecule)
		{
			for (molIt = molecule->begin (); molecule->end () != molIt; ++molIt)
			{
				if (0 != modelNumber)
				{
					--modelNumber;
				}
				else
				{
					annotate::AnnotateModel &am = (annotate::AnnotateModel&) *molIt;
					am.name(getFilePrefix(filename));
					annotate::AnnotationInteractions annInteractions;
		  			am.addAnnotation(annInteractions);
					am.annotate (gucRelationMask);
					std::string strTopology = topology(am);
					std::string strSymmetricTopology = symmetricTopology(am);
					std::map<std::string, std::list<std::string> >::iterator itTopo;
					itTopo = topologies.find(strTopology);
					if(itTopo == topologies.end())
					{
						itTopo = topologies.find(strSymmetricTopology);
					}

					if(itTopo != topologies.end())
					{
						itTopo->second.push_back(stripFilePath(filename));
					}else
					{
						std::pair<std::string, std::list<std::string> > entry;
						entry.first = strTopology;
						entry.second.push_back(stripFilePath(filename));
						topologies.insert(entry);
					}

					if (oneModel)
					{
						break;
					}
				}
			}
			delete molecule;
		}
		++optind;
	}

	std::map<std::string, std::list<std::string> >::const_iterator it;
	for(it = topologies.begin(); it != topologies.end(); ++ it)
	{
		mccore::gOut(0) << "----------------------------------------" << std::endl;
		mccore::gOut(0) << it->first;
		mccore::gOut(0) << "----------" << std::endl;
		std::list<std::string>::const_iterator itFile = it->second.begin();
		for(;itFile != it->second.end(); ++ itFile)
		{
			mccore::gOut(0) << *itFile << std::endl;
		}
	}
	mccore::gOut(0) << "Number of topologies : " << topologies.size() << std::endl;
	return EXIT_SUCCESS;
}
