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
	mccore::AbstractModel* AnnotateModelFM::createModel () const
	{
		return new AnnotateModel (residueSelection, environment, rFM);
	}


	mccore::AbstractModel* AnnotateModelFM::createModel(
		const mccore::AbstractModel &model) const
	{
		return new AnnotateModel (model, residueSelection, environment, rFM);
	}
  

	mccore::oBinstream& AnnotateModelFM::write (mccore::oBinstream& obs) const
	{
		return obs;
	}
  
	AnnotateModel::AnnotateModel (
  		const mccore::ResIdSet &rs, 
  		unsigned int env, 
  		const mccore::ResidueFactoryMethod *fm)
  	: mccore::GraphModel (fm),
		residueSelection (rs),
		environment (env)
    {
    	mucRelationMask = 0;
    }
    
    AnnotateModel::AnnotateModel (
    	const mccore::AbstractModel &right, 
    	const mccore::ResIdSet &rs, 
    	unsigned int env, 
    	const mccore::ResidueFactoryMethod *fm)
    : mccore::GraphModel (right, fm),
		residueSelection (rs),
		environment (env)
    {
    	mucRelationMask = 0;
    }
    
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
			mProvidedAnnotations.insert(aAnnotation.annotationName());
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
			if((*it)->annotationName() == astrAnnotName)
			{
				pAnnotation = *it;
			}
		}
		return pAnnotation;
	}


	void AnnotateModel::annotate (unsigned char aspb)
	{
		mucRelationMask = aspb;
		
		mccore::gOut(3) << "Computing basic annotation ..." << std::endl;
		GraphModel::annotate (aspb);
    
		// TODO : This should be moved into AnnotationCycle, but some const 
		// correctness work needs to be done in mccore.
		mccore::gOut(3) << "Computing minimum cycle bases union ..." << std::endl;
		unionMinimumCycleBases(mCyclesMolecule);

		// Compute all the requested annotations
		std::vector<Annotation*>::const_iterator it = annotations.begin();
		for(;it != annotations.end(); ++it)
		{
			mccore::gOut(3) << "Computing " << (*it)->annotationName();
			mccore::gOut(3) << " ..." << std::endl;
			(*it)->update(*this);
		}
	}
  
	void AnnotateModel::dumpConformations () const
	{
		const_iterator i;
    
		for (i = begin (); i != end (); ++i)
		{
			mccore::gOut(0) << i->getResId () << " : ";
			mccore::gOut(0) << mccore::Pdbstream::stringifyResidueType (i->getType ());
			if (i->getType ()->isNucleicAcid ())
			{
				mccore::gOut(0) << " " << i->getPucker();
				mccore::gOut(0) << " " << i->getGlycosyl();
			}
			mccore::gOut(0) << std::endl;
		}
	}

	ostream&
	AnnotateModel::output (ostream &os) const
	{
		mccore::gOut(0) << "Residue conformations -------------------------------------------" << endl;
		dumpConformations ();
	
		// Compute all the requested annotations
		std::vector<Annotation*>::const_iterator it = annotations.begin();
		for(;it != annotations.end(); ++it)
		{
			std::string strAnnotationName = (*it)->annotationName();
			strAnnotationName.append(" ");
			strAnnotationName.append(65 - strAnnotationName.size(), '-');
			mccore::gOut(0) << strAnnotationName << std::endl;
			mccore::gOut(0) << (*it)->output();
		}	
		return os;
	}

	mccore::iPdbstream& AnnotateModel::input (mccore::iPdbstream &is)
	{
		return GraphModel::input (is);
	}
  
	mccore::iBinstream& AnnotateModel::input (mccore::iBinstream &is)
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
	ostream& operator<< (ostream &os, const annotate::AnnotateModel &am)
	{
		return am.output (os);
	}
}