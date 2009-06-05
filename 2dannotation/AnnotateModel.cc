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
#include "mccore/Cycle.h"
#include "mccore/Messagestream.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/UndirectedGraph.h"
#include "mccore/stlio.h"

#include "AnnotateModel.h"
#include "Annotation.h"

//----------------------------------------------------------------------
template < class T >
set< T > SetIntersection( set< T > & set1, set< T > & set2 )
{
 set< T > setInter;
 insert_iterator< set< T > > iter( setInter, setInter.begin() );

 set_intersection( set1.begin(), set1.end(),
          set2.begin(), set2.end(),
          iter );

 return setInter;
};

//----------------------------------------------------------------------
template < class T >
set< T > SetDifference( set< T > & set1, set< T > & set2 )
{
 set< T > setDiff;
 insert_iterator< set< T > > iter( setDiff, setDiff.begin() );

 set_difference( set1.begin(), set1.end(),
        set2.begin(), set2.end(),
        iter );

 return setDiff;
};

namespace annotate
{
  static const unsigned int PAIRING_MARK = 1;
  
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
  AnnotateModel::annotate ()
  {
    basepairs.clear ();
    stacks.clear ();
    links.clear ();
    marks.clear ();

    GraphModel::annotate ();
    marks.resize (size (), 0);
    fillSeqBPStacks ();
    std::sort (basepairs.begin (), basepairs.end ());
    std::sort (stacks.begin (), stacks.end ());
    std::sort (links.begin (), links.end ());

	// Find the chains in the pdb file
	findChains();

	// Find the stems
	findStems();
	
	// Compute all the requested annotations
	std::vector<Annotation*>::const_iterator it = annotations.begin();
	for(;it != annotations.end(); ++it)
	{
		(*it)->update(*this);
	}
  }


  void
  AnnotateModel::fillSeqBPStacks ()
  {
    edge_iterator eit;
    
    for (eit = edge_begin (); edge_end () != eit; ++eit)
      {
	const Residue *ref;
	const Residue *res;

	if ((ref = (*eit)->getRef ())->getResId () < (res = (*eit)->getRes ())->getResId ())
	  {
	    GraphModel::label refLabel = getVertexLabel (const_cast< Residue* > (ref));
	    GraphModel::label resLabel = getVertexLabel (const_cast< Residue* > (res));
	
	    if ((*eit)->isPairing ())
	      {
		marks[refLabel] |= PAIRING_MARK;
		marks[resLabel] |= PAIRING_MARK;
		basepairs.push_back (BasePair (refLabel, ref->getResId (),
					       resLabel, res->getResId ()));
	      }
	    if ((*eit)->isStacking ())
	      {
		stacks.push_back (BaseStack (refLabel, ref->getResId (),
					     resLabel, res->getResId ()));
	      }
	    if ((*eit)->is (PropertyType::pAdjacent5p))
	      {
		links.push_back (BaseLink (refLabel, ref->getResId (),
					   resLabel, res->getResId ()));
	      }
	  }
      }
  }
  
  std::set< BasePair > AnnotateModel::getWWBasePairs()
  {
  	set< BasePair> oWWBasePairs;
  	vector< BasePair >::const_iterator bpit;

	for (bpit = basepairs.begin (); basepairs.end () != bpit; ++bpit)
	{
		const Relation &rel = *internalGetEdge (bpit->first, bpit->second);
		const vector< pair< const PropertyType*, const PropertyType* > > &faces = rel.getPairedFaces ();
		vector< pair< const PropertyType*, const PropertyType* > >::const_iterator pfit;

		for (pfit = faces.begin (); faces.end () != pfit; ++pfit)
	  	{
	  		const PropertyType* pProp1 = pfit->first;
	  		const PropertyType* pProp2 = pfit->second;
	  		std::string face1 = pProp1->toString();
	  		std::string face2 = pProp2->toString();
	  		if(	(0 < face1.size() && face1[0] == 'W') 
	  			&& (0 < face2.size() && face2[0] == 'W'))
	  		{
	  			BasePair oWWPair = *bpit;
	  			if(oWWPair.rResId < oWWPair.fResId)
	  			{
	  				oWWPair.reverse();	  				
	  			}
	  			oWWBasePairs.insert(oWWPair);
	  		}
		}
  	}
  	return oWWBasePairs;
  }
  
	void 
	AnnotateModel::findStems ()
	{
		std::set< BasePair > potentials;
		potentials = getWWBasePairs();
	
		Stem currentStem;
		std::set< BasePair >::const_iterator it = potentials.begin();
		for( ; it != potentials.end(); ++ it)
		{
			if(currentStem.continues(*it))
			{
				currentStem.push_back(*it);
			}
			else
			{
				if(1 < currentStem.size())
				{
					stems.push_back(currentStem);
				}
				currentStem.clear();
			}
		}
	
		if(1 < currentStem.size())
		{
			stems.push_back(currentStem);
		}
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

  
  void
  AnnotateModel::dumpStacks () const
  {
    vector< BaseStack > nonAdjacentStacks;
    vector< BaseStack >::const_iterator bsit;

    gOut(0) << "Adjacent stackings ----------------------------------------------" << endl;

    for (bsit = stacks.begin (); stacks.end () != bsit; ++bsit)
      {
	const Relation *rel = internalGetEdge (bsit->first, bsit->second);
	if (rel->is (PropertyType::pAdjacent))
	  {
	    const set< const PropertyType* > &labels = rel->getLabels ();

	    gOut (0) << bsit->fResId << "-" << bsit->rResId << " : ";
	    copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (gOut (0), " "));
	    gOut (0) << endl;
	  }
	else
	  {
	    nonAdjacentStacks.push_back (*bsit);
	  }
      }
    
    gOut(0) << "Non-Adjacent stackings ------------------------------------------" << endl;
    
    for (bsit = nonAdjacentStacks.begin (); nonAdjacentStacks.end () != bsit; ++bsit)
      {
	const set< const PropertyType* > &labels = internalGetEdge (bsit->first, bsit->second)->getLabels ();
	
	gOut (0) << bsit->fResId << "-" << bsit->rResId << " : ";
	copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (gOut (0), " "));
	gOut (0) << endl;
      }

    gOut(0) << "Number of stackings = " << stacks.size () << endl
	    << "Number of adjacent stackings = " << stacks.size () - nonAdjacentStacks.size () << endl
	    << "Number of non adjacent stackings = " << nonAdjacentStacks.size () << endl;
  }
  
	void
	AnnotateModel::dumpPair (const BasePair& aBasePair) const
	{
		const Relation &rel = *internalGetEdge (aBasePair.first, aBasePair.second);
		const set< const PropertyType* > &labels = rel.getLabels ();
		const vector< pair< const PropertyType*, const PropertyType* > > &faces = rel.getPairedFaces ();
		vector< pair< const PropertyType*, const PropertyType* > >::const_iterator pfit;

		gOut(0) << aBasePair.fResId << '-' << aBasePair.rResId << " : ";
		gOut(0) << Pdbstream::stringifyResidueType (rel.getRef ()->getType())
			<< "-"
			<< Pdbstream::stringifyResidueType (rel.getRes ()->getType ())
			<< " ";
		for (pfit = faces.begin (); faces.end () != pfit; ++pfit)
		{
			gOut (0) << *pfit->first << "/" << *pfit->second << ' ';
		}
		copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (gOut (0), " "));
		gOut (0) << endl;
	}
  

	void
	AnnotateModel::dumpPairs () const
	{
    	vector< BasePair >::const_iterator bpit;

		for (bpit = basepairs.begin (); basepairs.end () != bpit; ++bpit)
		{
			dumpPair(*bpit);
		}
	}
	
	void 
	AnnotateModel::dumpStems () const 
	{
  		int i=0;
    	vector< Stem >::const_iterator it;	
    
    	for (it = stems.begin (); stems.end () != it; ++it)
	    {
	    	ResId r1, r2, r3, r4;
	    	r1 = (*it).basePairs().front().fResId;
	    	r2 = (*it).basePairs().back().fResId;
	    	r3= (*it).basePairs().front().rResId;
	    	r4 = (*it).basePairs().back().rResId;
	    	gOut(0) << "Stem " << i << " : " << r1 << "-" << r2;
	    	gOut(0) << ", " << r3 << "-" << r4 << std::endl; 
	    	i ++;
    	}
	}
	
  ostream&
  AnnotateModel::output (ostream &os) const
  {
    gOut (0) << "Residue conformations -------------------------------------------" << endl;
    dumpConformations ();
    dumpStacks ();
    gOut (0) << "Base-pairs ------------------------------------------------------" << endl;
    dumpPairs ();
	gOut (0) << "Chains ----------------------------------------------------------" << endl;
	dumpChains();
	gOut (0) << "Stems -----------------------------------------------------------" << endl;
	dumpStems ();
	
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