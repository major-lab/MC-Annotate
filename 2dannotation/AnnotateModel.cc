//                              -*- Mode: C++ -*- 
// AnnotateModel.cc
// Copyright © 2001-06 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// $Revision: 58 $
// $Id: AnnotateModel.cc 58 2006-11-15 21:09:19Z larosem $


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <iterator>
#include <list>

#include "mccore/Binstream.h"
#include "mccore/Messagestream.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/UndirectedGraph.h"
#include "mccore/stlio.h"

#include "AnnotateModel.h"
#include "Annotation.h"
#include "AlgorithmExtra.h"

namespace annotate
{  
  AbstractModel* 
  AnnotateModelFM::createModel () const
  {
    return new AnnotateModel (residueSelection, environment, rFM);
  }


  AbstractModel*
  AnnotateModelFM::createModel (const AbstractModel &model) const
  {
    return new AnnotateModel (model, residueSelection, environment, rFM);
  }
  

  oBinstream& 
  AnnotateModelFM::write (oBinstream& obs) const
  {
    return obs;
  }
  
	AnnotateModel::AnnotateModel (
  		const ResIdSet &rs, 
  		unsigned int env, 
  		const ResidueFactoryMethod *fm)
  	: GraphModel (fm),
		residueSelection (rs),
		environment (env)
    {}
    
    AnnotateModel::AnnotateModel (
    	const AbstractModel &right, 
    	const ResIdSet &rs, 
    	unsigned int env, 
    	const ResidueFactoryMethod *fm)
    : GraphModel (right, fm),
		residueSelection (rs),
		environment (env)
    { }
    
    AnnotateModel::~AnnotateModel ()
    {
    	annotations.clear();
    }

	void AnnotateModel::addAnnotation(Annotation& aAnnotation) 
	{
		// Check dependencies
		std::set<std::string> requirements = aAnnotation.requires();
		
		std::set<std::string> match = SetIntersection<std::string>(
			requirements, 
			mProvidedAnnotations);
			
		if(match.size() == requirements.size())
		{
			annotations.push_back(&aAnnotation);
			mProvidedAnnotations.insert(aAnnotation.provides());
		}
		else
		{
			// TODO : throw an exception
		}
	}
	
	const Annotation* AnnotateModel::getAnnotation(
		const std::string& astrAnnotName) const
	{
		std::vector<Annotation *>::const_iterator it = annotations.begin();
		const Annotation* pAnnotation = NULL;
		for(; it != annotations.end() && NULL == pAnnotation; ++it)
		{
			if((*it)->provides() == astrAnnotName)
			{
				pAnnotation = *it;
			}
		}
		return pAnnotation;
	}


	void
	AnnotateModel::annotate (unsigned char aspb)
	{
		gOut (3) << "Computing basic annotation ..." << std::endl;
		GraphModel::annotate (aspb);
    
		// TODO : This should be moved into AnnotationCycle, but some const 
		// correctness work needs to be done in mccore.
		gOut (3) << "Computing minimum cycle bases union ..." << std::endl;
		unionMinimumCycleBases(mCyclesMolecule);

		// Compute all the requested annotations
		std::vector<Annotation*>::const_iterator it = annotations.begin();
		for(;it != annotations.end(); ++it)
		{
			gOut (3) << "Computing " << (*it)->provides() << " ..." << std::endl;
			(*it)->update(*this);
		}
		
		// Find the chains in the pdb file
		gOut (3) << "Computing chains ..." << std::endl;
		findChains();
	}
	 	
	void 
	AnnotateModel::findChains()
	{
		const_iterator i;
    	std::vector< const Residue *> current;
    
		for (i = begin (); i != end (); ++i)
		{			
			const Residue *res = &(*i);
			if(0 < current.size())
			{
				if(current.back()->getResId().getChainId() == res->getResId().getChainId())
				{
					current.push_back(res);
				}
				else
				{
					chains.push_back(current);
					current.clear();
				}
			}
			else
			{
				current.push_back(res);
			}
		}
		if(0 < current.size())
		{
			chains.push_back(current);
		}
	}

	void 
	AnnotateModel::dumpChains () const
	{
		std::vector< std::vector< const Residue * > >::const_iterator it;
		for(it = chains.begin(); it != chains.end(); ++ it)
		{
			gOut (0) << "Chain " << it->front()->getResId().getChainId() << " : ";
			std::vector< const Residue * >::const_iterator itRes;
			for(itRes = it->begin(); itRes != it->end(); )
			{
				gOut (0) << Pdbstream::stringifyResidueType ((*itRes)->getType ());
				++ itRes;
				if(itRes != it->end())
				{
					gOut (0) << "-";
				}
			}
			gOut (0) << endl;
		}		
	}

  
  void
  AnnotateModel::dumpConformations () const
  {
    const_iterator i;
    
    for (i = begin (); i != end (); ++i)
      {
	gOut(0) << i->getResId ()
		<< " : " << Pdbstream::stringifyResidueType (i->getType ());
	if (i->getType ()->isNucleicAcid ())
	  {
	    gOut (0) << " " << i->getPucker ()
		     << " " << i->getGlycosyl ();
	  }
	gOut (0) << endl;
      }
  }

  ostream&
  AnnotateModel::output (ostream &os) const
  {
    gOut (0) << "Residue conformations -------------------------------------------" << endl;
    dumpConformations ();
	gOut (0) << "Chains ----------------------------------------------------------" << endl;
	dumpChains();
	
	// Compute all the requested annotations
	std::vector<Annotation*>::const_iterator it = annotations.begin();
	for(;it != annotations.end(); ++it)
	{
		std::string strAnnotationName = (*it)->provides();
		strAnnotationName.append(" ");
		strAnnotationName.append(65 - strAnnotationName.size(), '-');
		gOut (0) << strAnnotationName << endl;
		std::string strOutput = (*it)->output();
		gOut (0) << strOutput;
	}
		
    return os;
  }


  iPdbstream&
  AnnotateModel::input (iPdbstream &is)
  {
    return GraphModel::input (is);
  }
  
  iBinstream&
  AnnotateModel::input (iBinstream &is)
  {
    return GraphModel::input (is);
  }
}

namespace std
{
  /**
   * Ouputs the residue to the stream.
   * @param os the output stream.
   * @param r the residue.
   * @return the used output stream.
   */
  ostream& 
  operator<< (ostream &os, const annotate::AnnotateModel &am)
  {
    return am.output (os);
  }  
}