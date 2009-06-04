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
	
    // Compute information on residues to be used for the annotations
    computeResidueInfos();
	
	// Find the linkers
	findLinkers();
	
	// Find the loops
	findLoops();
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
  
  	int AnnotateModel::getDirection(const StemConnection& aConnection) const
	{
		int iDirection = 0;
		if(aConnection.isValid())
		{
			switch(aConnection.getConnection())
			{
			case Stem::eFIRST_STRAND_FRONT_PAIR:
				iDirection = -1;
				break;
			case Stem::eFIRST_STRAND_BACK_PAIR:
				iDirection = 1;
				break;
			case Stem::eSECOND_STRAND_FRONT_PAIR:
				if(Stem::eANTIPARALLEL == aConnection.getStem().getOrientation())
				{
					iDirection = 1;
				}
				else
				{
					iDirection = -1;
				}
				break;
			case Stem::eSECOND_STRAND_BACK_PAIR:
				if(Stem::eANTIPARALLEL == aConnection.getStem().getOrientation())
				{
					iDirection = -1;
				}
				else
				{
					iDirection = 1;
				}
				break;
			case Stem::eUNDEFINED_CONNECTION:
				// TODO : Exception
				iDirection = 0;
			}
		}
		return iDirection;
	}
  
	Linker AnnotateModel::findLinker(const StemConnection& aConnection) const
	{
	  	Linker linker;
	  	
	  	// Find the possible strand from this end
	  	int iDirection = getDirection(aConnection);
	  	ResId resId = aConnection.getResidue();
	
	  	// Find the residue information
	  	std::vector<stResidueInfo>::const_iterator it = findResidueInfo(resId);
	  	
	  	std::vector<mccore::ResId> linkerResidues;
	  	
	  	// Set the start of the linker
	  	StemConnection linkerStart = aConnection;
	 	StemConnection linkerEnd;
	  	
	  	// Walk the chain until we find a stem not pseudo-knotted
	  	if(0 < iDirection)
	  	{
		  	std::advance(it, iDirection);
		  	while(it != mResidueInfos.end())
		  	{
		  		const Stem* pStem = (*it).pStem;
		  		if(NULL != pStem && !linkerStart.getStem().pseudoKnots(*pStem))
		  		{
		  			// Stem found
		  			Stem::enConnection eConnect = pStem->getConnection((*it).resId);
		  			linkerEnd = StemConnection(*pStem, eConnect);
		  			it = mResidueInfos.end();
		  		}
		  		else
		  		{
		  			linkerResidues.push_back((*it).resId);
		  			std::advance(it, iDirection);
		  		}
		  	}
	  	}
	  	linker = Linker(linkerResidues, linkerStart, linkerEnd);
	  	linker.order();
	  	
	  	return linker;
	}
  
  	void AnnotateModel::findLinker(
  		const Stem* apStem, 
  		const Stem::enConnection& aeConnect,
  		std::set<Linker>& outLinkerSet)
	{
	  	StemConnection connect(*apStem, aeConnect);
	  		
		Linker linker = findLinker(connect);
		
		if(!linker.isEmpty()) 
		{
			outLinkerSet.insert(linker); 
		}
  	}
  
  	void AnnotateModel::findLinkers()
  	{
	  	std::set<Linker> linkerSet;
	  	std::vector< Stem >::const_iterator it;
		for(it = stems.begin(); it != stems.end(); ++it)
		{
			findLinker((&*it), Stem::eFIRST_STRAND_FRONT_PAIR, linkerSet);
			findLinker((&*it), Stem::eFIRST_STRAND_BACK_PAIR, linkerSet);
			findLinker((&*it), Stem::eSECOND_STRAND_FRONT_PAIR, linkerSet);
			findLinker((&*it), Stem::eSECOND_STRAND_BACK_PAIR, linkerSet);		
		}
		
		std::set<Linker>::const_iterator itLinker = linkerSet.begin();
		for(;itLinker != linkerSet.end(); ++itLinker)
		{
			linkers.push_back(*itLinker);
		}
	}
	
	bool
	AnnotateModel::enclose(const BasePair& aBasePair, const Stem& aStem)
	{
		bool bEnclose = false;
		BasePair encEnd = aStem.basePairs().front();
		if(aBasePair.fResId < encEnd.fResId && encEnd.rResId < aBasePair.rResId)
		{
			bEnclose = true;
		}
		return bEnclose;		
	}
	
	std::vector<const Stem*> 
	AnnotateModel::getEnclosedStems(const BasePair& aBasePair)
	{
		std::vector< const Stem* > enclosed;
		std::vector< Stem >::const_iterator itEnc;
		for(itEnc = stems.begin(); itEnc != stems.end(); ++itEnc)
		{
			if(enclose(aBasePair, *itEnc))
			{
				enclosed.push_back(&(*itEnc));
			}
		}
		return enclosed;
	}
	
	void AnnotateModel::computeResidueInfos()
	{
		int i = 0;
		mResidueInfos.resize(GraphModel::size());
		const_iterator it = GraphModel::begin();
		for(;it != GraphModel::end(); ++ it)
		{
			mResidueInfos[i].resId = (*it).getResId();
			mResidueInfos[i].pResidue = &(*it);
			mResidueInfos[i].pStem = NULL;
			std::vector<Stem>::const_iterator stemIt;
			for(stemIt = stems.begin(); stemIt != stems.end(); ++stemIt)
			{
				if((*stemIt).contains((*it).getResId()))
				{
					if(NULL != mResidueInfos[i].pStem)
					{
						gOut(0) << "Residue associated with more than one stem" << endl;
					}
					mResidueInfos[i].pStem = &(*stemIt);
				}
			}
			++ i;
		}	
	}
	
	std::vector<AnnotateModel::stResidueInfo>::const_iterator 
	AnnotateModel::findResidueInfo(const ResId& aResId) const
	{
		std::vector<AnnotateModel::stResidueInfo>::const_iterator it;
		for(it = mResidueInfos.begin(); it != mResidueInfos.end(); ++it)
		{
			if(aResId == (*it).resId)
			{
				break;
			}
		}
		return it;
	}

	std::map< const Residue*, const Stem* >
	AnnotateModel::getResidueStemAssociation(unsigned int iChain) const
	{
		std::map< const Residue*, const Stem*> association;
		
		if(iChain < chains.size() && 0 < chains[iChain].size()) // This should be an assert
		{
			std::vector<const Residue*>::const_iterator itRes;
			itRes = chains[iChain].begin();
			for(; itRes != chains[iChain].end(); ++itRes)
			{
				std::pair< const Residue*, const Stem*> pairStem;
				const Stem* pStem = NULL;
				pairStem.first = *itRes;
				std::vector<Stem>::const_iterator itStem = stems.begin();
				for(; itStem != stems.end() && NULL == pStem; ++itStem)
				{
					if(itStem->contains(*(*itRes)))
					{
						pStem = &(*itStem);
					}
				}
				pairStem.second = pStem;
			}
		}
		return association;
	}
	
	ResId 
	AnnotateModel::nextId(
		const Stem& aStem, 
		const StemConnection& aConnection) const
	{
		ResId id;
		switch(aConnection.getConnection())
		{
		case Stem::eFIRST_STRAND_FRONT_PAIR:
			id = aStem.basePairs().front().rResId;
			break;
		case Stem::eSECOND_STRAND_FRONT_PAIR:
			id = aStem.basePairs().front().fResId;
			break;
		case Stem::eFIRST_STRAND_BACK_PAIR:
			id = aStem.basePairs().back().rResId;
			break;
		case Stem::eSECOND_STRAND_BACK_PAIR:
			id = aStem.basePairs().back().fResId;
			break;
		default:
			// TODO : This should be an exception
			break;
		}
		
		return id;
	}
	
	Linker AnnotateModel::nextLinker(
		const Linker& aLinker,
		const std::map<mccore::ResId, const Linker*>& aResidueLinkerMap)
	{
		Linker nextLinker;
		const Stem* pNextStem = &aLinker.getEnd().getStem();
		if(NULL != pNextStem)
		{
			ResId resId = nextId(
				aLinker.getEnd().getStem(), 
				aLinker.getEnd());
			std::map<mccore::ResId, const Linker*>::const_iterator it;
			it = aResidueLinkerMap.find(resId);
			if(it != aResidueLinkerMap.end())
			{
				const Linker* pLinker = (*it).second;
				nextLinker = *pLinker;
				
				if(nextLinker.getEnd().isValid())
				{
					ResId endId = nextLinker.getEnd().getResidue();
					if(it->first == endId)
					{
						nextLinker.reverse();	
					}					
				}				
			}
		}
		return nextLinker;
	}
	
	std::map<mccore::ResId, const Linker*>
	AnnotateModel::getResidueLinkerMap() const
	{
		std::map<mccore::ResId, const Linker*> residueLinkerMap;
		std::vector<Linker>::const_iterator it;
		for(it = linkers.begin();it != linkers.end(); ++it)
		{
			if((*it).getStart().isValid())
			{
				const StemConnection* pConnect = &(*it).getStart();
				ResId resId = pConnect->getResidue();
				std::pair<ResId, const Linker*> mapping(resId, &(*it));
				residueLinkerMap.insert(mapping);	
			}
			
			if((*it).getEnd().isValid())
			{
				const StemConnection* pConnection = &(*it).getEnd();
				ResId resId = pConnection->getResidue();
				std::pair<ResId, const Linker*> mapping(resId, &(*it));
				residueLinkerMap.insert(mapping);
			}
		}
		return residueLinkerMap;
	}
	
	void
	AnnotateModel::removeLinker(
		std::map<mccore::ResId, const Linker*>& aResidueLinkerMap,
		const Linker& aLinker)
	{
		if(aLinker.getEnd().isValid())
		{
			ResId resId = aLinker.getEnd().getResidue();
			aResidueLinkerMap.erase(resId);
		}
		
		if(aLinker.getStart().isValid())
		{
			ResId resId = aLinker.getStart().getResidue();
			aResidueLinkerMap.erase(resId);
		}
	}
	
	void
	AnnotateModel::findLoops()
	{
		std::map<ResId, const Linker*> residueLinkerMap;
		residueLinkerMap = getResidueLinkerMap();
		while(!residueLinkerMap.empty())
		{
			const Linker* pFirstLinker = NULL;
			std::vector<Linker> potentialLoop;
			
			// Pick the first linker
			std::pair<ResId, const Linker*> mapping = *(residueLinkerMap.begin());
			
			// Navigate until next is stem is NULL or we loop
			pFirstLinker = mapping.second;
			potentialLoop.push_back(*pFirstLinker);
			
			// Remove the linker from the map
			removeLinker(residueLinkerMap, *pFirstLinker);
			
			Linker currentLinker = nextLinker(
				*pFirstLinker, 
				residueLinkerMap);
			
			while(!currentLinker.isEmpty() && currentLinker != *pFirstLinker)
			{
				potentialLoop.push_back(currentLinker);
				
				// Remove the linker from the map
				removeLinker(residueLinkerMap, currentLinker);
				
				// Get the next strands
				currentLinker = nextLinker(
					currentLinker, 
					residueLinkerMap);				
			}
			
			if(0 < potentialLoop.size())
			{
				if(potentialLoop.front().getStart().isValid() && potentialLoop.back().getEnd().isValid())
				{
					const Stem* startStem = &potentialLoop.front().getStart().getStem();
					const Stem* endStem = &potentialLoop.back().getEnd().getStem();
					if(startStem == endStem)
					{
						loops.push_back(Loop(potentialLoop));
					}
				}
			}
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
	
	void
	AnnotateModel::dumpLoop(const Loop& aLoop) const
	{
		std::vector< Linker >::const_iterator it;
		for(it = aLoop.getLinkers().begin(); 
			it != aLoop.getLinkers().end(); 
			++ it)
		{
			gOut (0) << "{";
			dumpLinker(*it);
			gOut (0) << "}, ";
		}
	}
	
	void
	AnnotateModel::dumpLoops() const
	{
		int i = 0;
		std::vector< Loop >::const_iterator it;
		for(it = loops.begin(); it != loops.end(); ++it)
		{
			gOut (0) << "Loop " << i << " : ";
			dumpLoop(*it);
			gOut (0) << std::endl;
			++ i;
		}
	}
	
	void
	AnnotateModel::dumpLinker(const Linker& aLinker) const
	{
		gOut (0) << aLinker.getResidues().front();
		gOut (0) << "-"; 
		gOut (0) << aLinker.getResidues().back();
	}
	
	void 
	AnnotateModel::dumpLinkers() const
	{
		int i = 0;
		std::vector<Linker>::const_iterator it;
		for(it = linkers.begin(); it != linkers.end(); ++it)
		{
			gOut (0) << "Linker " << i << " : ";
			dumpLinker(*it);
			gOut (0) << std::endl;
			++ i;
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
	gOut (0) << "Strands ---------------------------------------------------------" << endl;
	dumpLinkers();
	gOut (0) << "Loops -----------------------------------------------------------" << endl;
	dumpLoops ();
		
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