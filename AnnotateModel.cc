//                              -*- Mode: C++ -*- 
// AnnotateModel.cc
// Copyright © 2001-06 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// $Revision$
// $Id$


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <cmath>
#include <iterator>
#include <list>
#include <time.h>

#include "mccore/Binstream.h"
#include "mccore/Messagestream.h"
#include "mccore/Pdbstream.h"
#include "mccore/UndirectedGraph.h"
#include "mccore/stlio.h"

extern "C" {
#include "cliquer.h"
#include "graph.h"
#include "set.h"
}

#include "AnnotateModel.h"



namespace annotate
{
  
  static const unsigned int MIN_HELIX_SIZE = 3;
  static const unsigned int PAIRING_MARK = 1;
  static const unsigned int LHELIX = 2;
  static const unsigned int RHELIX = 4;
  static unsigned int weight_GC = 100;
  static unsigned int weight_AU = 99;
  static unsigned int weight_GU = 95;
  static unsigned int weight_GA = 90;
  static unsigned int weight_AC = 90;
  static unsigned int weight_CU = 90;
  static unsigned int weight_GG = 90;
  static unsigned int weight_AA = 90;
  static unsigned int weight_CC = 90;
  static unsigned int weight_UU = 90;
  
  
  static void
  print_clique (set_t s, graph_t *g)
  {
    unsigned int i;

    printf ("size=%d, weight=%d:  ",set_size(s),graph_subgraph_weight(g,s));
    for (i=0; i<SET_MAX_SIZE(s); i++) {
      if (SET_CONTAINS(s,i)) {
	printf(" %d",i);
      }
    }
    printf("\n");
    return;
  }

  
  static boolean
  print_clique_func (set_t s,graph_t *g,clique_options *opts)
  {
    print_clique(s,g);
    return true;
  }

  
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
  
  
  bool 
  areLinked (const pair< AnnotateModel::label, AnnotateModel::label> &p , const pair< AnnotateModel::label, AnnotateModel::label> &q)
  {
    return (((p.first <= q.first) && (q.first <= p.second) && (p.second <= q.second))
	    || ((q.first <= p.first) && (p.first <= q.second) && (q.second <= p.second)));
  }


  const ResidueType*
  generalizeType (const ResidueType *type)
  {
    return (type->isA ()
	    ? ResidueType::rA
	    : (type->isC ()
	       ? ResidueType::rC
	       : (type->isG ()
		  ? ResidueType::rG
		  : (type->isU () || type->isT ()
		     ? ResidueType::rU
		     : ResidueType::rNull))));		     
  }


  unsigned int
  getWeight (const ResidueType *tf, const ResidueType *ts)
  {
    return ((tf == ResidueType::rG && ts == ResidueType::rC)
	    || (tf == ResidueType::rC && ts == ResidueType::rG)
	    ? weight_GC
	    : ((tf == ResidueType::rA && ts == ResidueType::rU)
	       || (tf == ResidueType::rU && ts == ResidueType::rA)
	       ? weight_AU
	       : ((tf == ResidueType::rG && ts == ResidueType::rU)
		  || (tf == ResidueType::rU && ts == ResidueType::rG)
		  ? weight_GU
		  : ((tf == ResidueType::rG && ts == ResidueType::rA)
		     || (tf == ResidueType::rA && ts == ResidueType::rG)
		     ? weight_GA
		     : ((tf == ResidueType::rC && ts == ResidueType::rA)
			|| (tf == ResidueType::rA && ts == ResidueType::rC)
			? weight_AC
			: ((tf == ResidueType::rC && ts == ResidueType::rU)
			   || (tf == ResidueType::rU && ts == ResidueType::rC)
			   ? weight_CU
			   : (tf == ts
			      ? (tf == ResidueType::rG
				 ? weight_GG
				 : (tf == ResidueType::rA
				    ? weight_AA
				    : (tf == ResidueType::rC
				       ? weight_CC
				       : (tf == ResidueType::rU
					  ? weight_UU
					  : 0))))
			      : 0)))))));
  }
  
  
  void
  AnnotateModel::s_secondaire (set< pair< AnnotateModel::label, AnnotateModel::label > > &stable)
  {
    vector< pair< AnnotateModel::label, AnnotateModel::label > > paire;
    vector< pair< AnnotateModel::label, AnnotateModel::label > >::iterator pit;
    AnnotateModel::label ref;
    graph_t *g;
    int i;
  
    stable.clear ();
    for (ref = 0; size () != ref; ++ref)
      {
	list< AnnotateModel::label > neighboors = internalNeighborhood (ref);
	list< AnnotateModel::label >::iterator rit;
      
	for (rit = neighboors.begin (); neighboors.end () != rit; ++rit)
	  {
	    AnnotateModel::label res = *rit;

	    if (ref < res && internalGetEdge (ref, res)->isPairing ())
	      {
		paire.push_back (make_pair (ref, res));
	      }
	  }
      }
    
    // creation du graphe de compatibilite
    g = graph_new (paire.size ());
  
    for (i = 0, pit = paire.begin (); paire.end () != pit; ++i, ++pit)
      {
// 	const ResidueType *ref = generalizeType (internalGetVertex (pit->first)->getType ());
// 	const ResidueType *res = generalizeType (internalGetVertex (pit->second)->getType ());
	
// 	g->weights[i] = getWeight (ref, res);
	g->weights[i] = 1;
      }
  
    unsigned int p = 0;
    
    for (pit = paire.begin (); paire.end () != pit; ++pit, ++p)
      {
	vector< pair< AnnotateModel::label, AnnotateModel::label > >::iterator qit = pit;
	unsigned int q;

	for (++qit, q = qit - paire.begin (); paire.end () != qit; ++qit, ++q)
	  {
	    if (! areLinked (*pit, *qit))
	      {
		GRAPH_ADD_EDGE (g, p, q);
	      }
	  }
      }

//     graph_print (g);
//     fflush (stdout);
//     cout << endl;

    // recherche de stable
    static int min_weight = 0;
    static int max_weight = 0;
    static boolean maximal = FALSE;
    clique_options *opts;
    
    opts = (clique_options*) malloc (sizeof (clique_options));
    opts->time_function = NULL;
    opts->output = stderr;
    opts->reorder_function = reorder_by_default;
    opts->reorder_map = NULL;
    opts->user_function = print_clique_func;
    opts->user_data = NULL;
    opts->clique_list = NULL;
    opts->clique_list_length = 0;
    
    set_t s = clique_find_single (g, min_weight, max_weight, maximal, opts);
  
    for (pit = paire.begin (), p = 0; paire.end () != pit; ++pit, ++p)
      {
	if (SET_CONTAINS (s, p))
	  {
	    stable.insert (*pit);
	  }
      }  
  }
  
  
  void
  AnnotateModel::annotate ()
  {
    time_t t;
    
    basepairs.clear ();
    stacks.clear ();
    links.clear ();
// //     helices.clear ();
// //     bulges.clear ();
// //     loops.clear ();
// //     internalloops.clear ();
// //     multiloops.clear ();
// //     singlestrands.clear ();
    marks.clear ();

// //     GraphModel::annotate (residueSelection);
    time (&t);
    GraphModel::annotate ();
    gOut (3) << "annotate " << time (0) - t << "s" << endl;

    time (&t);
    marks.resize (size (), 0);
    fillBPStacks ();
    gOut (3) << "fillBPStacks " << time (0) - t << "s" << endl;
    
    std::sort (basepairs.begin (), basepairs.end ());
    std::sort (stacks.begin (), stacks.end ());
//     buildSequences ();
    set< pair< AnnotateModel::label, AnnotateModel::label > > stable;
    set< pair< AnnotateModel::label, AnnotateModel::label > >::const_iterator sit;
    
    time (&t);
    s_secondaire (stable);
    gOut (3) << "s_secondaire " << time (0) - t << "s" << endl;
    
    gOut (0) << "Structure Secondaire:";
    for (sit = stable.begin (); stable.end () != sit; ++sit)
      {
	const pair< AnnotateModel::label, AnnotateModel::label > &paire = *sit;

	gOut (0) << " (" << internalGetVertex (paire.first)->getResId () << ' '
		 << internalGetVertex (paire.second)->getResId () << ")";
      }
    gOut (0) << endl;

//     findHelices ();
// //     findLoops ();
// //     findInternalLoops ();
// //     findMultiLoops ();
// //     findSingleStrands ();
  }


  void
  AnnotateModel::fillBPStacks ()
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
	
	    if ((*eit)->isPairing () || (*eit)->isBHbond ())
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
	  }
      }
  }


  bool
  AnnotateModel::isHelixPairing (const Relation &r)
  {
    return r.is (PropertyType::pPairing);
  }
    

  void 
  AnnotateModel::findHelices ()
  {
    list< BasePair > hPC;
    vector< BasePair >::const_iterator bpit;

    for (bpit = basepairs.begin (); basepairs.end () != bpit; ++bpit)
      {
	if (isHelixPairing (*internalGetEdge (bpit->first, bpit->second)))
	  {
	    hPC.push_back (*bpit);
	  }
      }
//     while (! hPC.empty ())
//       {
// 	Helix helix;
// 	list< BasePair >::iterator hpcIt;
// 	label hfst;
// 	label hsnd;
// 	label tfst;
// 	label tsnd;

// 	hpcIt = hPC.begin ();
// 	helix.push_back (*hpcIt);
// 	hpcIt = hPC.erase (hpcIt);
// 	hfst = helix.front ().first;
// 	hsnd = helix.front ().second;
// 	tfst = helix.back ().first;
// 	tsnd = helix.back ().second;
// 	while (hPC.end () != hpcIt)
// 	  {
// 	    label fst = hpcIt->first;
// 	    label snd = hpcIt->second;
// 	    bool inserted = false;
	    
// 	    if (((internalAreConnected (tfst, fst)
// 		  && internalGetEdge (tfst, fst)->isAdjacent ())
// 		 && (internalAreConnected (tsnd, snd)
// 		     && (internalGetEdge (tsnd, snd)->isAdjacent ())))
// 		|| ((internalAreConnected (tfst, snd)
// 		     && internalGetEdge (tfst, snd)->isAdjacent ())
// 		    && (internalAreConnected (tsnd, fst)
// 			&& internalGetEdge (tsnd, fst)->isAdjacent ())
// 		    && (hpcIt->reverse (), std::swap (fst, snd), true)))
// 	      {
// 		helix.push_back (*hpcIt);
// 		tfst = fst;
// 		tsnd = snd;
// 		inserted = true;
// 	      }
// 	    else if (((internalAreConnected (fst, hfst)
// 		       && internalGetEdge (fst, hfst)->isAdjacent ())
// 		      && (internalAreConnected (snd, hsnd)
// 			  && internalGetEdge (snd, hsnd)->isAdjacent ()))
// 		     || ((internalAreConnected (fst, hsnd)
// 			  && internalGetEdge (fst, hsnd)->isAdjacent ())
// 			 && (internalAreConnected (snd, hfst)
// 			     && internalGetEdge (snd, hfst)->isAdjacent ())
// 			 && (hpcIt->reverse (), std::swap (fst,snd), true)))
// 	      {
// 		helix.insert (helix.begin (), *hpcIt);
// 		hfst = fst;
// 		hsnd = snd;
// 		inserted = true;
// 	      }
// 	    if (inserted)
// 	      {
// 		hpcIt = hPC.erase (hpcIt);
// 		if (hPC.end () == hpcIt)
// 		  {
// 		    hpcIt = hPC.begin ();
// 		  }
// 	      }
// 	    else
// 	      {
// 		++hpcIt;
// 	      }
// 	  }
// 	if (MIN_HELIX_SIZE <= helix.size ())
// 	  {
// 	    Helix::iterator hit;
	    
// 	    helix.setId (helices.size ());
// 	    for (hit = helix.begin (); helix.end () != hit; ++hit)
// 	      {
// 		marks[hit->first] = LHELIX;
// 		marks[hit->second] = RHELIX;
// 	      }
// 	    helices.push_back (helix);
// 	  }
//       }
  }

  //     li = 0;
  // //     while (gi != end ())
  //     while (li < size ())
  //       {
  // 	// Find a pairing involving gi
  // // 	neighbor = neighborhood (&*gi);
  // 	neighbor = internalNeighborhood (li);
  // 	for (nborIt = neighbor.begin (); neighbor.end () != nborIt; ++nborIt)
  // 	  {
  // // 	    if (isHelixPairing (getEdge ((Residue*) &*gi, *nborIt)))
  // 	    if (isHelixPairing (*internalGetEdge (li, *nborIt)))
  // 	      {
  // 		break;
  // 	      }
  // 	  }

  // 	if (nborIt != neighbor.end ())
  // 	  {
  // 	    Residue *ri = internalGetVertex (li);
  // 	    Residue *rj = internalGetVertex (*nborIt);
	    
  // 	    if (ri->getResId () < rj->getResId () 
  // 		&& marks.end () != marks.find (rj)
  // 		&& marks[rj] != '('
  // 		&& marks[rj] != ')')
  // 	      {
  // 		helix.push_back (make_pair (ri, rj));
          
  // 		// Extending with bulge detection
  // 		while (true)
  // 		  {
  // 		    lip++ = li;
  // 		    if (gip == size ()
  // 		      {
  // 			break;
  // 		      }
  // 		    if (gk == begin ())
  // 		      {
  // 			break;
  // 		      }
  // 		    gkp-- = gk;
        
  // 		    if (sequence_mask[&*gip] == sequence_mask[&*gi]
  // 			&& sequence_mask[&*gkp] == sequence_mask[&*gk])
  // 		      {
  // 			neighbor = neighborhood ((Residue*) &*gip);
  // 			if ((nborIt = std::find (neighbor.begin (), neighbor.end (), (Residue*) &*gkp)) != neighbor.end ()
  // 			    && isHelixPairing (getEdge ((Residue*) &*gip, (Residue*) &*gkp)))
  // 			  {
  // 			    helix.push_back (make_pair ((Residue*) &*gip, (Residue*) &*gkp));
  // 			    gi = gip;
  // 			    gk = gkp;
  // 			  }
  // 			else
  // 			  {
  // 			    break;
  // 			  }
  // 		      }
  // 		    else
  // 		      {
  // 			break;
  // 		      }
  // 		  }
		
  // 		if ((int) helix.size () >= MIN_HELIX_SIZE)
  // 		  {
  // 		    Helix::iterator i;
		    
  // 		    for (i = helix.begin (); i != helix.end (); ++i)
  // 		      {
  // 			helix_mask[i->first] = helices.size ();
  // 			marks[i->first] = '(';
  // 			helix_mask[i->second] = helices.size ();
  // 			marks[i->second] = ')';
  // 		      }
  // 		    helices.push_back (helix);
  // 		  } 
  // 		helix.clear ();
  // 	      }
  // 	  }
  // 	++gi;
  //       }
  //   }

  /**
   * Output helices in text representation
   * Ex:
   * H0, length = 14
   *      A1-UUAUAUAUAUAUAA-A14
   *      B14-AAUAUAUAUAUAUU-B1
   */
  void
  AnnotateModel::dumpHelices () const 
  {
    vector< Helix >::const_iterator i;
    
    for (i = helices.begin (); helices.end () != i; ++i)
      {
	Helix::const_iterator hIt;
	
	// Helix index and length
	gOut (0) << "H" << i - helices.begin () << ", length = " << i->size ()
		 << endl;
	
	// First strand
	gOut (0).setf (ios::right, ios::adjustfield);
	gOut (0) << setw (6) << internalGetVertex (i->front ().first)->getResId () << "-";
	for (hIt = i->begin (); i->end () != hIt; ++hIt)
	  {
	    const ResidueType *type = internalGetVertex (hIt->first)->getType ();
	    
	    gOut (0) << (type->isNucleicAcid ()
			 ? Pdbstream::stringifyResidueType (type)
			 : "X");
	  }
	--hIt;
	gOut (0) << "-" << internalGetVertex (hIt->first)->getResId () << endl;

	// Second strand 
	gOut (0).setf (ios::right, ios::adjustfield);
	gOut (0) << setw (6) << internalGetVertex (i->front ().second)->getResId () << "-";
	for (hIt = i->begin (); i->end () != hIt; ++hIt)
	  {
	    const ResidueType *type = internalGetVertex (hIt->second)->getType ();
	    
	    gOut (0) << (type->isNucleicAcid ()
			 ? Pdbstream::stringifyResidueType (type)
			 : "X");
	  }
	--hIt;
	gOut (0) << "-" << internalGetVertex (hIt->second)->getResId() << endl;
      }
    gOut (0).setf (ios::left, ios::adjustfield);
  }

  
  void
  AnnotateModel::findStrands ()
  {

//     const_iterator gi, gk, gl;
//     map< const Residue*, int >::iterator gf;
//     OStrand strand;

//     for (gi=begin (); gi!=end (); )
//     {
//       if (helix_mask.find(&(*gi)) == helix_mask.end()) {
      
//         gk = gi;
//         while (gk != end () &&
//             helix_mask.find(&(*gk)) == helix_mask.end() &&
//             sequence_mask[&(*gi)] == sequence_mask[&(*gk)]
//             )
//           ++gk;
//         gk--;
//         strand.first = &(*gi);
//         strand.second = &(*gk);
//         strand.type = OTHER;
//         for (gl = gi; gl<=gk; ++gl) { 
//           strand_mask[&(*gl)] = strands.size ();
//         }
//         strands.push_back (strand);
//         gi = gk;
//       }
//       gi++;
//     }
  }

    
  void 
  AnnotateModel::buildStrands ()
  {
//     list< edge_const_iterator > unsorted;
//     edge_const_iterator edgeIt;

//     sequences.clear ();
//     for (edgeIt = edge_begin (); edge_end () != edgeIt; ++edgeIt)
//       {
// 	if ((*edgeIt)->is (PropertyType::pAdjacent5p))
// 	  {
// 	    unsorted.push_back (edgeIt);
// 	  }
//       }
//     while (! unsorted.empty ())
//       {
// 	edge_const_iterator sj;
// 	list< edge_const_iterator > sorted;
// 	list< edge_const_iterator >::iterator it;
    
// 	sj = unsorted.front ();
// 	unsorted.pop_front ();
// 	sorted.push_back (sj);
// 	for (it = unsorted.begin (); unsorted.end () != it;)
// 	  {
// 	    if ((**it)->getRes ()->getResId () == (*sj)->getRef ()->getResId ())
// 	      {
// 		sj = *it;
// 		sorted.push_front (sj);
// 		unsorted.erase (it);
// 		it = unsorted.begin ();
// 	      }
// 	    else
// 	      {
// 		++it;
// 	      }
// 	  }
// 	sj = sorted.back ();
// 	for (it = unsorted.begin (); unsorted.end () != it;)
// 	  {
// 	    if ((**it)->getRef ()->getResId () == (*sj)->getRes ()->getResId ())
// 	      {
// 		sj = *it;;
// 		sorted.push_back (sj);
// 		unsorted.erase (it);
// 		it = unsorted.begin ();
// 	      }
// 	    else
// 	      {
// 		++it;
// 	      }
// 	  }
// 	sequences.push_back (Strand ());
// 	Strand &seq = sequences.back ();
    
// 	it = sorted.begin ();
// 	seq.push_back ((**it)->getRef ());
// 	while (sorted.end () != it)
// 	  {
// 	    seq.push_back ((**it)->getRes ());
// 	    ++it;
// 	  }
//       }

//     // Verify that the sequences have the same chain id and have consecutive
//     // residue number.
//     StrandSet::iterator seqIt;
//     int k = 0;

//     for (seqIt = sequences.begin (); sequences.end () != seqIt; ++seqIt, ++k)
//       {
// 	Strand &resVec = *seqIt;
// 	Strand::iterator resVecIt;
// 	char chain;
// 	int number;
// 	int len = 1;
    
// 	resVecIt = resVec.begin ();
// 	chain = (*resVecIt)->getResId ().getChainId ();
// 	number = (*resVecIt)->getResId ().getResNo ();
// 	sequence_mask[*resVecIt] = k;
// 	++resVecIt;
// 	while (resVec.end () != resVecIt)
// 	  {
// 	    if ((*resVecIt)->getResId ().getChainId () != chain
// 		|| (*resVecIt)->getResId ().getResNo () != number + 1)
// 	      {
// 		Strand rest;
      
// 		if (++resVec.begin () == resVecIt)
// 		  {
// 		    rest.insert (rest.end (), resVecIt, resVec.end ());
// 		    resVec.erase (++resVecIt, resVec.end ());
// 		  }
// 		else
// 		  {
// 		    rest.insert (rest.end (), resVecIt - 1, resVec.end ());
// 		    resVec.erase (resVecIt, resVec.end ());
// 		  }
// 		++seqIt;
// 		if (1 < rest.size ())
// 		  {
// 		    seqIt = --sequences.insert (seqIt, rest);
// 		  }
// 		break;
// 	      }
// 	    else
// 	      {
// 		sequence_mask[*resVecIt] = k;
// 		number += 1;
// 		++resVecIt;
// 		len++;
// 	      }
// 	  }
// 	sequence_length.push_back(len);
//       }
  }

//   void AnnotateModel::classifyStrands ()  {
//     vector< OStrand >::iterator j;
//     const_iterator first, last, prev, next;
  
//     for (j=strands.begin (); j!=strands.end (); ++j) {
  
//       prev = find((*(j->first)).getResId());
//       next = find((*(j->second)).getResId());
//       --prev; ++next;
//       first = begin();
//       last = end();
//       --last;
  
//       if (&(*(j->first)) == &(*first) || &(*(j->second)) == &(*last) ||
//         sequence_mask[&(*prev)] != sequence_mask[&(*next)])
//       {
//         j->type = OTHER;
//       }
//       else if (helix_mask[&(*prev)] == helix_mask[&(*next)])
//       {
//         const Residue *i = 0;
//         // Find pairs of the extrimities.
//         if (helices[helix_mask[&(*prev)]].front ().first == &(*prev))
//           i = helices[helix_mask[&(*prev)]].front ().second;
//         else if (helices[helix_mask[&(*prev)]].front ().second == &(*prev))
//           i = helices[helix_mask[&(*prev)]].front ().first;
//         else if (helices[helix_mask[&(*prev)]].back ().first == &(*prev))
//           i = helices[helix_mask[&(*prev)]].back ().second;
//         else if (helices[helix_mask[&(*prev)]].back ().second == &(*prev))
//           i = helices[helix_mask[&(*prev)]].back ().first;
//         else 
//           cerr << "Pairing not found for (a) " << *prev << endl;
  
//         if (i == &(*next))
//           j->type = LOOP;
//         else 
//           j->type = BULGE_OUT;
//       }
//       else
//       { 
//         const Residue *i = 0, *k = 0;
    
//         // Find pairs.
//         if (helices[helix_mask[j->first-1]].front ().first == j->first-1)
//           i = helices[helix_mask[j->first-1]].front ().second;
//         else if (helices[helix_mask[j->first-1]].front ().second == j->first-1)
//           i = helices[helix_mask[j->first-1]].front ().first;
//         else if (helices[helix_mask[j->first-1]].back ().first == j->first-1)
//           i = helices[helix_mask[j->first-1]].back ().second;
//         else if (helices[helix_mask[j->first-1]].back ().second == j->first-1)
//           i = helices[helix_mask[j->first-1]].back ().first;
//         else 
//           cerr << "Pairing not found for (b) " << j->first-1 << endl;
      
//         if (helices[helix_mask[j->second+1]].front ().first == j->second+1)
//           k = helices[helix_mask[j->second+1]].front ().second;
//         else if (helices[helix_mask[j->second+1]].front ().second == j->second+1)
//           k = helices[helix_mask[j->second+1]].front ().first;
//         else if (helices[helix_mask[j->second+1]].back ().first == j->second+1)
//           k = helices[helix_mask[j->second+1]].back ().second;
//         else if (helices[helix_mask[j->second+1]].back ().second == j->second+1)
//           k = helices[helix_mask[j->second+1]].back ().first;
//         else 
//           cerr << "Pairing not found for (c) " << j->first-1 << endl;
      
//         if (sequence_mask[i] != sequence_mask[k]) {
//           j->type = OTHER;
//         } else if (i==k+1 || i==k-1) {
//           j->type = BULGE;
//         } else {
//           // Find the other part of the internal loop
//           vector< OStrand >::iterator s;
//           for (s=strands.begin (); s!=strands.end (); ++s) {
//             if (s->first == k+1 && s->second == i-1) {
//               j->type = INTERNAL_LOOP;
//               j->ref = s-strands.begin ();  // FIXIT!
//               break;
//             } else {
//               j->type = OTHER;
//             }
//           }
//         }
//       }    
//     }
//   }
  
//   void AnnotateModel::dumpStrands ()
//   {
//     int j;
//     const_iterator k;
  
//     for (j=0; j<(int)strands.size (); ++j) {
      
//       gOut(0).setf (ios::left, ios::adjustfield);
//       gOut(0) << "S" << setw (5) << j << " ";
      
//       gOut(0).setf (ios::left, ios::adjustfield);
      
//       if (strands[j].type == OTHER) gOut(0) << setw (16) << "single strand:";
//       else if (strands[j].type == BULGE) gOut(0) << setw (16) << "bulge:";
//       else if (strands[j].type == BULGE_OUT) gOut(0) << setw (16) << "bulge out:";
//       else if (strands[j].type == LOOP) gOut(0) << setw (16) << "loop:";
//       else if (strands[j].type == INTERNAL_LOOP) gOut(0) << setw (16) << "internal loop:";
      
//       gOut(0) << (strands[j].first)->getResId() << "-";
      
//       for (k = find((strands[j].first)->getResId()); k <= find((strands[j].second)->getResId()); ++k) {
//         gOut(0) << ((k)->getType()->isNucleicAcid () ? Pdbstream::stringifyResidueType ((k)->getType()) : "X");
//       }
//       gOut(0) << "-" << (strands[j].second)->getResId();
  
//       if (strands[j].type == INTERNAL_LOOP) {
//         gOut(0) << " -- ";
//       // FIXIT!    
//   /*      gOut(0) << (strands[strands[j].ref].first)->getResId() << "-";
//         for (k=strands[strands[j].ref].first; k<=strands[strands[j].ref].second; ++k) {
//           gOut(0) << ((k)->getType()->isNucleicAcid () ? Pdbstream::stringifyResidueType ((k)->getType()) : "X");
//         }
//         gOut(0) << "-" << (strands[strands[j].ref].second)->getResId();
//   */      
//       }
//       gOut(0) << endl;
//     }
//     gOut(0).setf (ios::left, ios::adjustfield);  
//   }

  
  void
  AnnotateModel::findKissingHairpins ()
  {
//     const_iterator gi;
//     list< Residue* > neighbor;
//     list< Residue* >::iterator nborIt;
//     int check = 0;
//     int nb_helical_bp = 0;
//     set< const PropertyType* >::iterator k;
  
//     gOut(0).setf (ios::left, ios::adjustfield);
  
//     for (gi=begin (); gi!=end (); ++gi) {
//       neighbor = neighborhood (const_cast<Residue*> (&(*gi)));  
//       for (nborIt = neighbor.begin (); nborIt != neighbor.end (); ++nborIt) {      
//         if (gi->getResId() < (*nborIt)->getResId() && isPairing (getEdge ((Residue*)(&(*gi)), *nborIt))) {
          
//           if (helix_mask.find(&(*gi)) != helix_mask.end() &&
//               helix_mask [&(*gi)] == helix_mask [*nborIt]) {
//             // Simple helical base-pair.
//             gOut(0) << setw (20) << "helix bp: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
//             check++;
//             nb_helical_bp++;
//           } 
      
//           if (helix_mask.find(&(*gi)) != helix_mask.end() &&
//               helix_mask [&(*gi)] != helix_mask [*nborIt]) {
//             gOut(0) << setw (20) << "inter helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
//             tertiary_mask[&(*gi)] = 1;
//             tertiary_mask[*nborIt] = 1;
//             check++;
//           }
  
//           if (strand_mask.find(&(*gi)) != strand_mask.end() &&
//               helix_mask.find(&(*gi)) != helix_mask.end()) {
//             if (strands[strand_mask [&(*gi)]].type == OTHER) {
//               gOut(0) << setw (20) << "strand/helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
// 	      tertiary_mask[&(*gi)] = 1;
// 	      tertiary_mask[*nborIt] = 1;
//               check++;
//             } else if (strands[strand_mask [&(*gi)]].type == BULGE ||
// 		       strands[strand_mask [&(*gi)]].type == BULGE_OUT) {
//               gOut(0) << setw (20) << "bulge/helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
// 	      tertiary_mask[&(*gi)] = 1;
// 	      tertiary_mask[*nborIt] = 1;
//               check++;
//             } else if (strands[strand_mask [&(*gi)]].type == LOOP) {
//               gOut(0) << setw (20) << "loop/helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
// 	      tertiary_mask[&(*gi)] = 1;
// 	      tertiary_mask[*nborIt] = 1;
//               check++;
//             } else if (strands[strand_mask [&(*gi)]].type == INTERNAL_LOOP) {
//               gOut(0) << setw (20) << "internal loop/helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
// 	      tertiary_mask[&(*gi)] = 1;
// 	      tertiary_mask[*nborIt] = 1;
//               check++;
//             } else {
//               gOut(0) << "Other A: " << " : ";
//             }
//           }
  
//           if (strand_mask.find(&(*gi)) != strand_mask.end() &&
//               strand_mask.find(*nborIt) != strand_mask.end()) {
            
//             if (strand_mask [&(*gi)] == strand_mask [*nborIt]) {
//               gOut(0) << setw (20) << "intraloop: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
// 	      tertiary_mask[&(*gi)] = 1;
// 	      tertiary_mask[*nborIt] = 1;
//               check++;
//             } else {             
//               int size = 0;
//               if (strands[strand_mask [&(*gi)]].type == OTHER) {
//                 gOut(0) << "strand/";
//                 size = 20-7;	    
//               } else if (strands[strand_mask [&(*gi)]].type == BULGE || 
// 			 strands[strand_mask [&(*gi)]].type == BULGE_OUT) {
//                 gOut(0) << "bulge/";
//                 size = 20-6;
//               } else if (strands[strand_mask [&(*gi)]].type == LOOP) {
//                 gOut(0) << "loop/";
//                 size = 20-5;
//               } else if (strands[strand_mask [&(*gi)]].type == INTERNAL_LOOP) {
//                 gOut(0) << "internal loop/";
//                 size = 20-14;
//               }              
//               if (strands[strand_mask [*nborIt]].type == OTHER) {
//                 gOut(0) << setw (size) << "strand: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
//                 tertiary_mask[&(*gi)] = 1;
//                 tertiary_mask[*nborIt] = 1;
//                 check++;
// 	      } else if (strands[strand_mask [*nborIt]].type == BULGE ||
// 			 strands[strand_mask [*nborIt]].type == BULGE_OUT) {
//                 gOut(0) << setw (size) << "bulge: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
//                 tertiary_mask[&(*gi)] = 1;
//                 tertiary_mask[*nborIt] = 1;
//                 check++;
//               } else if (strands[strand_mask [*nborIt]].type == LOOP) {
//                 gOut(0) << setw (size) << "loop: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
//                 tertiary_mask[&(*gi)] = 1;
//                 tertiary_mask[*nborIt] = 1;
//                 check++;
//               } else if (strands[strand_mask [*nborIt]].type == INTERNAL_LOOP) {
//                 gOut(0) << setw (size) << "internal loop: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
//                 tertiary_mask[&(*gi)] = 1;
//                 tertiary_mask[*nborIt] = 1;
//                 check++;
//               } else {
//                 gOut(0) << "Other B:" << " : ";
//               }
//             }
//           }
         
//           gOut(0) << Pdbstream::stringifyResidueType (gi->getType()) << "-"
//                   << Pdbstream::stringifyResidueType ((*nborIt)->getType()) << " ";
          
//           Relation *e = getEdge ((Residue*) &(*gi), *nborIt);
                   
//           if (isPairing(e))
//             gOut(0) << e->getRefFace () << "/" << e->getResFace () << " ";
 
//           const set< const PropertyType* > &labels = e->getLabels();
//           copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (gOut(0), " "));          
//           gOut(0) << endl;       
//         }
//       }
//     }
//     gOut(0) << "Number of base pairs = " << nb_pairings << endl;
//     gOut(0) << "Number of helical base pairs = " << nb_helical_bp << endl;
  
//     if (nb_pairings != check) 
//       gOut(0) << "Missing interactions: " << check << "/" << nb_pairings << endl;
  }


  void
  AnnotateModel::findPseudoknots ()
  {
//   vector< Helix >::iterator i, j;
  
//   for (i=helices.begin (); i!=helices.end (); ++i) {
//     if (strand_mask[i->front ().first] == 
// 	strand_mask[i->front ().second]) { 
//       const Residue *a, *b, *c, *d;
//       a = (*i).front ().first;
//       b = (*i).back ().first;
//       c = (*i).back ().second;
//       d = (*i).front ().second;
//       for (j=i; j!=helices.end (); ++j) {
// 	if (i!=j && strand_mask[j->front ().first] == 
// 	    strand_mask[j->front ().second]) {
// 	  const Residue *ap, *bp, *cp, *dp;
// 	  ap = (*j).front ().first;
// 	  bp = (*j).back ().first;
// 	  cp = (*j).back ().second;
// 	  dp = (*j).front ().second;
	  
// 	  if (b->getResId() < ap->getResId() && bp->getResId() < c->getResId() && d->getResId() < cp->getResId()) {
// 	    gOut(0) << "Pseudoknot : " << "H" << helix_mask[a]+1 << " " 
// 		 << "H" << helix_mask[ap]+1 << endl;
// 	  }
// 	}
//       }
//     }
//   }
  }
 
  void
  AnnotateModel::dumpSequences (bool detailed) const
  {
    vector< Sequence >::const_iterator it;
    unsigned int currseq;

    gOut (0) << "Sequences -------------------------------------------------------"
	     << endl;
    for (it = sequences.begin (), currseq = 0; sequences.end () != it; ++it, ++currseq)
      {
	const Sequence &sequence = *it;
	Sequence::const_iterator sit;

	gOut (0) << "Sequence " << currseq << " (length = " 
		 << sequence.size () << "): " << endl;
	gOut(0).setf (ios::right, ios::adjustfield);
	gOut(0) << setw (2) << " ";
	gOut(0).setf (ios::left, ios::adjustfield);
	gOut(0) << setw (7) << internalGetVertex (*sequence.begin ())->getResId();
	for (sit = sequence.begin (); sequence.end () != sit; ++sit)
	  {
	    gOut (0) << Pdbstream::stringifyResidueType (internalGetVertex (*sit)->getType ());
	  }
	gOut (0) << endl;
      }
//   iterator i;
//   int j = 0;
//   int currseq = -1;
//   int size = 0;
  
//   if (!detailed) {
// //    for (i=0; i<(int)conformations.size (); ++i) {
// //      if (i==0 || sequence_mask[i] != sequence_mask[i-1])
// //	gOut(0) << " " << setw (8) << getResId (i) << " ";
// //      gOut(0) << getType (i)->toString ();
// //    }	
//     return;
//   } 
  
//   for (i=begin (); i!=end (); ) {
    
//     currseq = sequence_mask[(Residue*) &(*i)];
//     size = sequence_length[currseq];

//     gOut(0) << "Sequence " << currseq << " (length = " 
// 	 << size << "): " << endl;
//     int pos = 0;
//     while (pos < size) 
//       {
// 	for (j=pos; (j<size && (j+1)%50!=0); ++j) {
// 	  if ((j)%50==0) {
// 	    gOut(0).setf (ios::right, ios::adjustfield);
// 	    gOut(0) << setw (3) << (*(i+j))->getResId().getChainId ();
// 	    gOut(0).setf (ios::left, ios::adjustfield);
// 	    gOut(0) << setw (5) << (*(i+j))->getResId().getResNo () << " ";
// 	  }
	  
// 	  gOut(0) << ((*(i+j))->getType()->isNucleicAcid () ? Pdbstream::stringifyResidueType ((*(i+j))->getType()) : "X");
	  
// 	  if ((j+1)%10==0) gOut(0) << " ";
// 	}
// 	gOut(0) << endl;
	
// 	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
// 	  if ((j)%50==0) gOut(0) << setw (8) << " " << " ";
//  	  if (marks.find(*(i+j)) == marks.end()) gOut(0) << "-";
//           else gOut(0) << marks[*(i+j)];
// 	  if ((j+1)%10==0) gOut(0) << " ";
// 	}
// 	gOut(0) << endl;

// //  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
// //  	  if ((j)%50==0)
// //          gOut(0) << setw (6) << "seq" << " ";
// //  	  if (sequence_mask[i+j] == -1)
// //          gOut(0) << "-";
// //  	  else gOut(0) << sequence_mask[i+j];
// //  	  if ((j+1)%10==0)
// //          gOut(0) << " ";
// //  	}
// //  	gOut(0) << endl;

//  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
//  	  if ((j)%50==0) gOut(0) << setw (8) << "helix" << " ";
//  	  if (helix_mask.find(*(i+j)) == helix_mask.end()) gOut(0) << "-";
//  	  else gOut(0) << helix_mask[*(i+j)];
//  	  if ((j+1)%10==0) gOut(0) << " ";
//  	}
//  	gOut(0) << endl;

// //  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
// //  	  if ((j)%50==0) gOut(0) << setw (6) << "str" << " ";
// //  	  if (strand_mask[i+j] == -1) gOut(0) << "-";
// //  	  else gOut(0) << strand_mask[i+j];
// //  	  if ((j+1)%10==0) gOut(0) << " ";
// //  	}
// //  	gOut(0) << endl;

// // 	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
// // 	  if ((j)%50==0) gOut(0) << setw (6) << "ter" << " ";
// // 	  if (tertiary_mask[i+j] == -1) gOut(0) << "-";
// // 	  else gOut(0) << 'X';
// // 	  if ((j+1)%10==0) gOut(0) << " ";
// // 	}

// // 	gOut(0) << endl << endl;

// 	pos = j+1;
//       }
//     i += pos-1;
//   }
//   gOut(0).setf (ios::left, ios::adjustfield);
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
// 	    << "Number of helical stackings = " << nb_helical_stacks << endl
	    << "Number of adjacent stackings = " << stacks.size () - nonAdjacentStacks.size () << endl
	    << "Number of non adjacent stackings = " << nonAdjacentStacks.size () << endl;
  }
  

  void
  AnnotateModel::dumpPairs () const
  {
    vector< BasePair >::const_iterator bpit;

    for (bpit = basepairs.begin (); basepairs.end () != bpit; ++bpit)
      {
	const Relation &rel = *internalGetEdge (bpit->first, bpit->second);
	const set< const PropertyType* > &labels = rel.getLabels ();
	const vector< pair< const PropertyType*, const PropertyType* > > &faces = rel.getPairedFaces ();
	vector< pair< const PropertyType*, const PropertyType* > >::const_iterator pfit;

	gOut(0) << bpit->fResId << '-' << bpit->rResId << " : ";
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
  }

  
  void
  AnnotateModel::dumpTriples ()
  {
//     const_iterator gi;
//     list< Residue* > neighbor;
//     list< Residue* >::iterator nborIt;

//     int count;
//     int nb = 1;
//     set< node > nodes;
//     set< node >::iterator k, l;
//     vector< bool > treated;

//     treated.resize (size (), false);

//     for (gi=begin (); gi!=end (); ++gi) {
//       count = 0;
//       nodes.clear ();
//       if (!treated[*gi]) {
// 	nodes.insert (*gi);
      
// 	neighbor = getNeighbors (*gi);

// 	for (nborIt=neighbor.begin (); nborIt!=neighbor.end (); ++nborIt) {
// 	  if (isPairing (getEdge (*gi,*nborIt))) {
// 	    count++;
// 	    nodes.insert (*nborIt);
// 	  }
// 	}
// 	if (count > 1) {
// 	  gOut(0).setf (ios::left, ios::adjustfield);
// 	  gOut(0) << "T" << setw (5) << nb++ << " "; 
// 	  for (k=nodes.begin (); k!=nodes.end (); ++k) {	
// 	    for (l=k; l!=nodes.end (); ++l) {
// 	      if (k!=l && isPairing (*k, *l)) {
// 		gOut(0) << getResId (*k) << "-"
// 			<< getResId (*l) << " ";
// 		treated[*k] = true;
// 		treated[*l] = true;
// 	      }
// 	    }
// 	  }
// 	  gOut(0) << endl;
// 	}
//       }
//     }
//     gOut(0).setf (ios::left, ios::adjustfield);
  }

  ostream&
  AnnotateModel::output (ostream &os) const
  {
    gOut (0) << "Residue conformations -------------------------------------------" << endl;
    dumpConformations ();
    dumpStacks ();
    gOut (0) << "Base-pairs ------------------------------------------------------" << endl;
// //     findKissingHairpins ();
    dumpPairs ();
//     gOut (0) << "Triples ---------------------------------------------------------" << endl;
// //     dumpTriples ();
//     gOut (0) << "Helices ---------------------------------------------------------" << endl;
//     dumpHelices ();
//     gOut (0) << "Strands ---------------------------------------------------------" << endl;
//     dumpStrands ();
//     gOut (0) << "Various features ------------------------------------------------" << endl;
// //     findPseudoknots ();
    dumpSequences ();
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
  
//   oBinstream&
//   AnnotateModel::output (oBinstream &os) const
//   {
//     return os;
//   }
 
//   iBinstream&
//   operator>> (iBinstream &is, AnnotateModel &model)
//   {
//     return model.input (is);
//   }

//   oBinstream&
//   operator<< (oBinstream &os, const AnnotateModel &model)
//   {
//     return model.output (os);
//   }

}

namespace std
{
  ostream &
  operator<< (ostream &out, const annotate::Strand &t)
  {
    return t.output (out);
  }

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
  
// namespace annotate {

// void AnnotateModel::PDF_drawLoop2 (PDF *p, int li, int source_jct)
// {
//   int font = PDF_findfont (p, "Helvetica", "host", 0);

//   Vector rect (helix_width,unit_length*4);

//   int rad = unit_length*3;

//   // Drawing current loop ------------------------------------------------------
//   Block *loop = &l_blocks[li];
//   loop->output (cout); cout << endl;

//   PDF_circle (p, 0, 0, rad);
//   PDF_stroke (p);
  
//   for (int i=0; i!=(int)loop->front().size (); ++i)
//     {
//       char b[2];
//       sprintf (b, "%c", getType (loop->front ()[i])->getOneLetterRep());
      
//       PDF_set_text_pos (p, rad*sin (i*2*M_PI/loop->front ().size ()), rad*cos (i*2*M_PI/loop->front ().size ()));

//       PDF_show (p, b);
//     }
  
//   if (source_jct>=0) { PDF_rotate (p, 180);
//     cout << "Rotating by " << 180 << endl;
//   }
  
//   // Draw the connected helices ------------------------------------------------
//   pair< node, node > *jip;
//   int ji;
//   float theta_step = 360/loop->junctions.size ();
//   float l = rect.y/2 + sqrt (rad*rad - (rect.x/2)*(rect.x/2));
  
//   for (ji=0; ji<(int)loop->junctions.size (); ++ji) {      
//     int actual_ji = (source_jct==-1)?ji:(ji+source_jct)%loop->junctions.size ();
//     jip = &(loop->junctions[actual_ji]);
    
//     if (actual_ji!=source_jct)  {
//       Block* helix;
//       cout << "Looking for helix (" << getResId (jip->first) << ", " << getResId (jip->second) << ")" << endl;
//       cout << "Found helix " << helix_mask[jip->first] << " :"; 
      
//       helix = &h_blocks[helix_mask[jip->first]];
//       helix->output (cout); cout << endl;
      
//       // Move to center of helix.
//       cout << "Rotating by " << 360-(theta_step * ji) << endl;
//       PDF_rotate (p, 360-(theta_step * ji));      
//       PDF_translate (p, 0, l);
//       PDF_rect (p, -rect.x/2, -rect.y/2, rect.x, rect.y);
//       PDF_stroke (p);	

//       cout << "RECT = " << rect.x << " " <<  rect.y << endl;
      
//       for (int i=0; i!=(int)helix->front().size (); ++i)
// 	{
// 	  char b[2];
// 	  sprintf (b, "%c", getType (helix->front ()[i])->getOneLetterRep());
// 	  cout << b << " " << -rect.y/2 + i*rect.y/helix->front ().size () << endl;
// 	  cout << PDF_stringwidth (p, b, font, 8) << endl;;

// 	  PDF_set_text_pos (p, -rect.x/2 - PDF_stringwidth (p, b, font, 8)/2 , -rect.y/2 + i*rect.y/helix->front ().size ());
// 	  PDF_show (p, b);
// 	}
      
//       // Move to center of new loop.
//       PDF_translate (p, 0, l);	
      
//       pair< node, node > junc;
//       if ((helix->junctions[0].first == jip->second &&
// 	   helix->junctions[0].second == jip->first)) {
// 	junc = helix->junctions[1];
// 	cout << "first case" << endl;
//       } else if (helix->junctions[0].first == jip->first &&
// 		 helix->junctions[0].second == jip->second) {
// 	cerr << "Looking for helix (" << getResId (jip->first) << ", " << getResId (jip->second) << ")" << endl;
// 	cerr << "but found " << getResId (helix->junctions[0].first) << ", " 
// 	     << getResId (helix->junctions[0].second) << ")" << endl;	  
//       } else {
// 	cout << "second case" << endl;
// 	junc = helix->junctions[0];
//       }
      
//       // Finding loop at the other extrimity and do a recursive call.
      
//       cout << "Looking for loop (" << getResId (junc.first) << ", " << getResId (junc.second) << ")" << endl;
      
//       vector< Block >::iterator bi;
//       vector< pair< node, node > >::iterator bj;
//       for (bi=l_blocks.begin (); bi!=l_blocks.end (); ++bi) {
// 	for (bj=bi->junctions.begin (); bj!=bi->junctions.end (); ++bj) {
// 	  if((bj->first==junc.first && bj->second==junc.second) ||
// 	     (bj->first==junc.second && bj->second==junc.first)) {	   
// 	      PDF_drawLoop2 (p, bi-l_blocks.begin (), bj-bi->junctions.begin ());
// 	      //PDF_circle (p, 0, 0, rad/2);
// 	      //PDF_stroke (p);	      
// 	  }
// 	}
//       }
      
//       cout << "DeRotating by " << (theta_step * ji) << endl;
//       PDF_translate (p, 0, -l);
//       PDF_translate (p, 0, -l);
//       PDF_rotate (p, theta_step * ji);  
      
//     }
//   }
//   if (source_jct>=0) {
//     PDF_rotate (p, 180);
//     cout << "DeRotating by " << 180 << endl;
//   }
// }



// void AnnotateModel::PDF_drawHelix (PDF *p, int hi, int li_ref, int x, int y)
// {
//   int font = PDF_findfont (p, "Helvetica", "host", 0);
  
//   Block *helix = &h_blocks[hi];
//   cout << "Drawing helix " << h_blocks[hi] << endl;

//   // Calculate size of helix
//   int length = helix->front ().size ();
//   Vector rect (helix_width,unit_length*(length-1));
  
//   // Find radius of reference loop
//   int rad = 30;

//   // Calculate distance from center of reference loop
//   float l = 0;
//   if (li_ref!=-1) {
//     l = rect.y/2 + sqrt (rad*rad - (rect.x/2)*(rect.x/2));
//   }
  
//   // Draw helix
//   PDF_translate (p, 0, l);
//   PDF_rect (p, -rect.x/2, -rect.y/2, rect.x, rect.y);
//   PDF_stroke (p);

//   //  float spacing = rect.y/(helix->front ().size ()-1);

//   {
//     char b[8];
//     sprintf (b, "H%d", hi);
//     PDF_set_text_pos (p, 0, 0);
//     PDF_show (p, b);
//   }
// //   for (int i=0; i!=(int)helix->front().size (); ++i)
// //     {
// //       char b[2];
// //       sprintf (b, "%c", getType (helix->front ()[i])->getOneLetterRep());
// //       PDF_set_text_pos (p, -rect.x/2 - PDF_stringwidth (p, b, font, 8)/2 , -rect.y/2 + i*unit_length);
// //       PDF_show (p, b);

// //       sprintf (b, "%c", getType (helix->back ()[i])->getOneLetterRep());
// //       PDF_set_text_pos (p, rect.x/2 - PDF_stringwidth (p, b, font, 8)/2 , rect.y/2 - i*unit_length);
// //       PDF_show (p, b);
// //     }
  
//   // Draw connected loops

//   vector< Block >::iterator b;
//   vector< pair< node, node > >::iterator i, j;
//   for (i=helix->junctions.begin (); i!=helix->junctions.end (); ++i) {
//     for (b=l_blocks.begin (); b!=l_blocks.end (); ++b) {
//       if (b-l_blocks.begin () != li_ref) {
// 	for (j=b->junctions.begin (); j!=b->junctions.end (); ++j) {
// 	  if ((j->first==i->first && j->second==i->second) ||
// 	      (j->first==i->second && j->second==i->first)) {
// 	    PDF_drawLoop (p, b-l_blocks.begin (), hi);
// 	  }
// 	}
//       }
//     }
//   }

//   PDF_translate (p, 0, -l);
  
//   cout << "END" << endl;
// }



// void AnnotateModel::PDF_drawLoop (PDF *p, int li, int hi_ref, int x, int y)
// {
//   Block *loop = &l_blocks[li];
//   cout << "Drawing loop " << l_blocks[li] << endl;

//   // Calculate size of loop
//   int rad = 30;

//   // Find length of reference helix
//   int length;
//   Vector rect;
//   if (hi_ref!=-1) {    
//     length = h_blocks[hi_ref].front ().size ()-1;
//     rect = Vector (helix_width,unit_length*length);
//   }

//   // Calculate distance from middle of reference helix
//   float l = 0;
//   if (hi_ref!=-1) {
//     l = rect.y/2 + sqrt (rad*rad - (rect.x/2)*(rect.x/2));
//   }

//   PDF_translate (p, 0, l);
//   PDF_circle (p, 0, 0, rad);
//   PDF_stroke (p);

//   char b[8];
//   sprintf (b, "L%d", li);
//   PDF_set_text_pos (p, 0, 0);
//   PDF_show (p, b);

// //   for (int i=0; i!=(int)loop->front().size (); ++i)
// //   {
// //     char b[2];
// //     sprintf (b, "%c", getType (loop->front ()[i])->getOneLetterRep());
    
// //     float a = (float)rect.x / rad;
// //     cout << rect.x << " " << rad << " " << (float)rect.x/rad << endl;

// //     PDF_set_text_pos (p, rad*sin (i*2*M_PI/loop->front ().size ()), rad*cos (i*2*M_PI/loop->front ().size ()));    
// //     PDF_show (p, b);
// //   }

//   // Find which helix to draw
//   int ji;
//   float theta_step = 360/loop->junctions.size ();

//   for (ji=0; ji<(int)loop->junctions.size (); ++ji) {
//     if (helix_mask[loop->junctions[ji].first] != hi_ref) {
//       cout << helix_mask[loop->junctions[ji].first] << " " << flush 
// 	   << h_blocks[helix_mask[loop->junctions[ji].first]] << endl;
      
//       PDF_rotate (p, 360-(theta_step * ji));      
//       PDF_drawHelix (p, helix_mask[loop->junctions[ji].first], li);
//       PDF_rotate (p, theta_step * ji);  
//     }
//   }
//   PDF_translate (p, 0, -l);
// }



// void AnnotateModel::dumpPDF (const char* pdfname)
// {
//   float width = 500;//letter_width;
//   float height = 500;//letter_height;
//   Vector v (width/2, height/2);

//   cout << "Writing PDF document" << endl;

//   PDF *p = PDF_new (); 
 
//   PDF_open_file (p, pdfname);

//   int font = PDF_findfont (p, "Helvetica", "host", 0);

//   PDF_begin_page (p, width, height);
//   PDF_setfont (p, font, 8);

//   vector< Block >::iterator bi, bj;

//   PDF_save (p);
//   PDF_translate (p, v.x, v.y);
//   //  PDF_drawLoop (p, 2);
//   PDF_drawHelix (p, 0);
//   PDF_restore (p);

// //   if (l_blocks.size () > 0) {
// //     // Draw the largest loop in the center first!
// //     bj = l_blocks.begin ();
// //     int lsize = bj->junctions.size ();
// //     for (bi=l_blocks.begin (); bi!=l_blocks.end (); ++bi) {
// //       if ((int)bi->junctions.size () > lsize) {
// // 	bj = bi;
// // 	lsize = bi->junctions.size ();
// //       }
// //     }
    
// //     Vector v (width/2, height/2);
    
// //     PDF_save (p);
// //     PDF_translate (p, v.x, v.y);
    
// //     PDF_drawLoop2 (p, bj-l_blocks.begin (), -1);
    
// //     PDF_restore (p);
// //   }
// //   else if (h_blocks.size ()>0) {
// //     cerr <<  "There is only an helix!" << endl;

// //     Vector rect (helix_width, unit_length*4);
// //     Vector v (width/2, height/2);
// //     PDF_save (p);
// //     PDF_translate (p, v.x, v.y);
// //     PDF_rect (p, -rect.x/2, -rect.y/2, rect.x, rect.y);
// //     PDF_stroke (p);	
// //     PDF_restore (p);
// //   }  

//   PDF_end_page (p);
//   PDF_close (p);
//   PDF_delete (p);
// }



// // void AnnotateModel::PDF_drawGraph (PDF *p, node curr, node prev, int x, int y, int a)
// // {
// //   int font = PDF_findfont (p, "Helvetica", "host", 0);

// //   char b[2];
// //   sprintf (b, "%c", getType (curr)->getOneLetterRep ());
// //   PDF_set_text_pos (p, x, y);
// //   PDF_show (p, b);

// // //   PDF_setfont (p, font, 6);
// // //   char q[8];
// // //   cout << getResId (curr) << endl;
// // //   sprintf (q, "%s", (const char*)getResId (curr));
// // //   PDF_show (p, q);
// // //   PDF_setfont (p, font, 12);

// //   graph.unmarkNode (curr);
  
// //   Graph< node, edge >::adjlist n;
// //   Graph< node, edge >::adjlist::iterator i;
// //   n = graph.getNeighborIds (curr);
  
// //   for (i=n.begin (); i!=n.end (); ++i) {
// // //     cout << getType (i->first)->getOneLetterRep () << " " 
// // // 	 << graph.isMarked (i->first) << endl;

// //     if (graph.isMarked (i->first)) {
// //       if (isPairing (i->first, curr)) {
// // 	//PDF_drawGraph (p, i->first, curr, x+10, y, a);
// //       } else if (isAdjacent (graph.getEdge(i->first, curr))){
// // 	PDF_drawGraph (p, i->first, curr, x, y+10, a);
// //       }
// //     }
// //   }

// //   graph.markNode (curr);
// // }


// // void AnnotateModel::dumpSimplePDF (const char* pdbname, const char* pdfname)
// // {
// //   float width = 500;//letter_width;
// //   float height = 500;//letter_height;

// //   cout << "Writing PDF document" << endl;

// //   PDF *p = PDF_new (); 
 
// //   PDF_open_file (p, pdfname);

// //   int font = PDF_findfont (p, "Helvetica", "host", 0);

// //   PDF_begin_page (p, width, height);
// //   PDF_setfont (p, font, 12);

// //   cout << pdbname << endl;
// //   PDF_set_text_pos (p, 2, height-10);
// //   PDF_show (p, pdbname);
// //   PDF_show (p, pdbname);

// //   PDF_continue_text (p, pdbname);
  

// //   int i;
// //   for (i=0; i<(int)conformations.size (); ++i) {
// //     if (i==0 || sequence_mask[i] != sequence_mask[i-1]) {
// //       PDF_continue_text (p, getResId (i));
// //     }
// //     PDF_continue_text (p, *getType (i));
// //   }	
  

// //   PDF_drawGraph (p, 0, -1, width/2, height/2, 0);

// //   PDF_end_page (p);
// //   PDF_close (p);
// //   PDF_delete (p);
// // }

// }
