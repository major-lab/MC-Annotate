//                              -*- Mode: C++ -*- 
// AnnotateModel.cc
// Copyright � 2001-04 Laboratoire de Biologie Informatique et Th�orique.
//                     Universit� de Montr�al
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// $Revision$
// $Id$
//


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// #include <iomanip>
#include <algorithm>
#include <list>
#include <utility>

#include "mccore/Algo.h"
// #include "mccore/Exception.h"
#include "mccore/Pdbstream.h"
#include "mccore/stlio.h"

#include "AnnotateModel.h"



namespace annotate
{
  
  bool
  Contact (Residue &a, Residue &b, float min_cut, float max_cut)
  {
    min_cut *= min_cut;
    max_cut *= max_cut;
    
    for (Residue::iterator i = a.begin (); i != a.end (); ++i)
      {
	for (Residue::iterator j = b.begin (); j != b.end (); ++j)
	  {
	    float d = i->squareDistance(*j);
	    if (d < min_cut) return true;
	    else if (d > max_cut) return false;
	  }
      }
    return false;
  }


  AnnotateModel::AnnotateModel (string &name, PdbFileHeader &header, Model &m)
    : name (name), fileHeader (header), nb_pairings (0), nb_connect (0)
  {
    Model::iterator i;
    unsigned int a;
    map< unsigned int, ResId > corrMap;
    vector< pair< Model::iterator, Model::iterator > > contacts;
    vector< pair< Model::iterator, Model::iterator > >::iterator l;
    list< unsigned int > seq;
    vector< vector< Residue* > > sequences;
  
    // Standardization of ResIDs
    // This is needed for the extractContacts algo to work correctly when
    // there are multiple residues of the same resid.
    for (i = m.begin (), a = 1; m.end () != i; ++i, ++a)
      {
	ResId resid (a);

	insert (i->clone ());
	corrMap[a] = i->getResId ();  
	internalGetNode (a)->setResId (resid);
      }
    // Extraction of possible relations based on the Axis Aligned Bounding
    // Box method
    contacts = Algo::extractContacts (begin (), end (), 5.0);
    for (l = contacts.begin (); contacts.end () != l; ++l)
      {
	Relation rel (&*l->first, &*l->second);
	
	if (rel.annotate ())
	  {
	    if (rel.is (PropertyType::pDIR_3p))
	      {
		connect (&*l->first, &*l->second, rel.invert ()->clone ());
		seq.push_back (edgeSize () - 1);
	      }
	    else
	      {
		connect (&*l->first, &*l->second, rel->clone ());
		if (rel.is (PropertyType::pDIR_5p))
		  {
		    seq.push_back (edgeSize () - 1);	  	    
		  }
	      }
	  }
      }
    // Selection sort of the sequence.
    sequenceSort (seq, sequences);
    findLoneResidues (sequences);
    rebuildGraph (sequences);
  }


  void
  AnnotateModel::sequenceSort (list< unsigned int > &seq, vector< vector< Residue* > > &sequences)
  {
    while (! seq.empty ())
      {
	unsigned int sj;
	list< unsigned int > sorted_seq;
	list< unsigned int >::iterator it;
	
	sj = seq.front ();
	seq.pop_front ();
	sorted_seq.push_back (sj);
	for (it = seq.begin (); seq.end () != it; ++it)
	  {
	    if (internalGetEdge (*it)->getRes ()->getResId () == internalGetEdge (sj)->getRef ()->getResId ())
	      {
		sj = *it;
		sorted_seq.push_front (sj);
		seq.erase (it);
		it = seq.begin ();
	      }
	  }

	sj = sorted_seq.back ();
	for (it = seq.begin (); seq.end () != it; ++it)
	  {
	    if (internalGetEdge (*it)->getRef ()->getResId () == internalGetEdge (sj)->getRes ()->getResId ())
	      {
		sj = *it;
		sorted_seq.push_back (sj);
		seq.erase (it);
		it = seq.begin ();
	      }
	  }

	sequences.push_back (vector< Residue* > ());
	vector< Residue* > &actual_seq = sequences.back ();
	
	it = sorted_seq.begin ();
	actual_seq.push_back (internalGetEdge (*it)->getRef ());
	while (sorted_seq.end () != it)
	  {
	    actual_seq.push_back (internalGetEdge (*it)->getRes ());
	    it++;
	  }
      }
  }


  void
  AnnotateModel::findLoneResidues (vector< vector< Residue* > > &sequences)
  {
    Model::iterator i;
    ResIdSet modelSet;
    vector< vector< Residue* > >::iterator vvit;
    ResIdSet sequenceSet;
    ResIdSet tmp;
    ResIdSet::iterator risIt;

    for (i = begin (); end () != i; ++i)
      {
	modelSet.insert ((*i)->getResId ());
      }
    for (vvit = sequences.begin (); sequences.end () != vvit; ++vvit)
      {
	vector< Residue* >::iterator vit;

	for (vit = vvit->begin (); vvit->end () != vit; ++vit)
	  {
	    sequenceSet.insert ((*vit)->getResId ());
	  }
      }
    set_difference (modelSet.begin (), modelSet.end (),
		    sequenceSet.begin (), sequenceSet.end (),
		    inserter (tmp, tmp.begin ()));

    for (risIt = tmp.begin (); tmp.end () != risIt; ++risIt)
      {
	sequences.push_back (vector< Residue* > ());
	vector< Residue* > &seq = sequences.back ();
	seq.push_back (&internalGetNode (risIt->getResNo ()));
      }
  }


  void
  AnnotateModel::rebuildGraphModel (vector< vector< Residue* > > &sequences)
  {
    unsigned int a;
    vector< vector< Residue* > >::iterator vvit;
    GraphModel newGraph;
    map< unsigned int, ResId > &newMap;
    vector< Relation* >::iterator rit;

    for (vvit = sequences.begin (), a = 1; sequences.end () != vvit; ++vvit)
      {
	vector< Residue* >::iterator vit;

	for (vit = vvit->begin (); vvit->end () != vit; ++vit, ++a)
	  {
	    Residue *res = *vit;
	    ResId id (a);
	    
	    newMap[a] = int2ResIdMap[res->getResId ().getResNo ()];
	    newGraph.insert (res);
	    res->setResId (id);
	    sequence_mask.push_back (vvit - sequences.begin ());
	  }
      }
    
    marks.resize (size (), '-');
    helix_mask.resize (size (), -1);
    strand_mask.resize (size (), -1);
    tertiary_mask.resize (size (), -1);

    for (rit = edges.begin (); edges.end () != rit; ++rit)
      {
	newGraph.connect ((*rit)->getRef (), (*rit)->getRes (), *rit);

	if ((*rit)->is (PropertyType::pDIR_3p)
	    || (*rit)->is (PropertyType::pDIR_5p))
	  {
	    nb_connect++;
	  }
	if ((*rit)->is (PropertyType::pPairing))
	  {
	    marks[(*rit)->getRef ()->getResId ().getResNo ()] = 'b';
	    marks[(*rit)->getRes ()->getResId ().getResNo ()] = 'b';
	    nb_pairings++;
	  }
      }
    *this = newGraph;
    int2ResIdMap = newMap;
  }
  
  
  void 
  AnnotateModel::findHelices ()
  {
    UndirectedGraph< node, edge >::iterator gi, gk, gip, gkp;
    list< node >::iterator gj, gjp;
    Helix helix;
    list< node > neigh;
    
    int min_size = 3;
    
    gi = graph.begin ();
    
    while (gi!=graph.end ())
      {      
	//cout << "( " << getResId (*gi) << " " << flush;
      
	// Find a pairing involving gi
	neigh = graph.getNeighbors (*gi);
	for (gj=neigh.begin (); gj!=neigh.end (); ++gj)
	  {
	    if (isHelixPairing (graph.getEdge (*gi, *gj))) break;
	  }
      
	if (gj != neigh.end ()) 
	  {
	    if (*gj > *gi && marks[*gj] != '(' && marks[*gj] != ')')
	      {
		// Find counterpart of gi.
		gk = graph.find (*gj);

		//cout << getResId (*gk) << ")" << endl;

		helix.push_back (make_pair (*gi, *gk));
	  
		// Extending with bulge detection
		while (true)
		  {
		    gip = gi;  ++gip;
		    if (gip == graph.end ()) break;
		    if (gk == graph.begin ()) break;
		    gkp = gk;  --gkp;
	    
		    //cout << "?  " << getResId(*gip) << " " << getResId(*gkp) << endl;

		    if (sequence_mask[*gip] == sequence_mask[*gi] &&
			sequence_mask[*gkp] == sequence_mask[*gk])
		      {
	      
			neigh = graph.getNeighbors (*gip);
			if ((gjp = std::find (neigh.begin (), neigh.end (), *gkp)) != neigh.end ()
			    && isHelixPairing (graph.getEdge (*gip, *gkp)))
			  {
			    helix.push_back (make_pair (*gip, *gkp));
			    gi = gip;
			    gk = gkp;
		
			    //cout << "+ " << getResId(*gi) << " " << getResId(*gk) << endl;

			  }
			else
			  {
			    break;
			  }
// This is for allowing bulges in the helix...  Doesn't work for the following case:
// '0'110-'0'51 : C-G Ww/Ww pairing cis XIX 
// '0'111-'0'50 : C-G Ww/Ww pairing cis XIX 
// '0'112-'0'49 : G-A Ss/Hh pairing trans XI 
// '0'113-'0'47 : A-G Hh/Ss pairing trans XI 
// '0'114-'0'48 : A-A Hh/Hh pairing trans II
// 	      } else {

// 		cout << "Trying bulge in second part" << endl;
// 		// There may be a bulge in one of the two chains
//  		Graph< node, edge >::adjlist::iterator ga, gb;
// 		bool found = false;
// 		for (ga=gip->second.begin (); ga!=gip->second.end (); ++ga) {

//  		  cout << getResId (ga->first) << " " << getResId (gk->first) << " " 
//  		       << "bulge length " << gk->first - ga->first << endl;

// 		  if (isHelixPairing (ga->second) &&
// 		      gk->first - ga->first <= min_size && 
// 		      (gb = gk->second.find (ga->first)) != gk->second.end () &&
// 		      relations[gb->second].isStacking ()) {
// 		    helix.push_back (make_pair (-1, ga->first-gk->first));
// 		    helix.push_back (make_pair (gip->first, gb->first));
// 		    gi = gip;
// 		    gk = graph.find (gb->first);

//   		    cout << "+ " << getResId(gi->first) << " " << getResId(gk->first) << endl;

// 		    found = true;
// 		    break;
// 		  }
// 		}
//   		cout << "Trying bulge in first part" << endl;
// 		if (!found) {
// 		  for (ga=gkp->second.begin (); ga!=gkp->second.end (); ++ga) {

//  		    cout << getResId (ga->first) << " " << getResId (gi->first) << " " 
//  			 << ga->first - gi->first << endl;

// 		    if (isHelixPairing (ga->second) &&
// 			ga->first - gi->first <= min_size && 
// 			(gb = gi->second.find (ga->first)) != gi->second.end () &&
// 			relations[gb->second].isStacking ()) {
// 		      helix.push_back (make_pair (gi->first-ga->first, -1));
// 		      helix.push_back (make_pair (gb->first, gkp->first));
// 		      gi = graph.find (gb->first);
// 		      gk = gkp;
//   		      cout << "+ " << getResId(gi->first) << " " << getResId(gk->first) << endl;
// 		      found = true;
// 		      break;
// 		    }
// 		  }
// 		}
// 		if (!found) {
// 		  break;
// 		}
// 	      }
		      }
		    else
		      {
			break;
		      }
		  }
	  
		if ((int)helix.size () >= min_size)
		  {
		    Helix::iterator i;
		    for (i=helix.begin (); i!=helix.end (); ++i)
		      {
			if (i->first >= 0)
			  {
			    helix_mask[i->first] = helices.size ();
			    marks[i->first] = '(';
			    helix_mask[i->second] = helices.size ();
			    marks[i->second] = ')';
			  }
		      }
		    helices.push_back (helix);
		  } 
		helix.clear ();
	      }
	  }
	++gi;
      }
  }


  void
  AnnotateModel::findStrands ()
  {
    int i;
    UndirectedGraph< node, edge >::iterator gi, gk;
    Strand strand;

    for (gi=graph.begin (); gi!=graph.end (); )
      {
	if (helix_mask[*gi] == -1)
	  {
	    gk = gi;
	    while (gk != graph.end ()
		   && sequence_mask[*gi] == sequence_mask[*gk]
		   && helix_mask[*gk] == -1)
	      {
		++gk;
	      }
	    gk--;
	    strand.first = *gi;
	    strand.second = *gk;
	    strand.type = OTHER;
	    for (i=*gi; i<=*gk; ++i)
	      {
		strand_mask[i] = strands.size ();
	      }
	    strands.push_back (strand);
	    gi = gk;
	  }
	gi++;
      }
  }


  void
  AnnotateModel::classifyStrands ()
  {
    vector< Strand >::iterator j;

    for (j=strands.begin (); j!=strands.end (); ++j)
      {
	if (j->first == 0
	    || j->second == (int)graph.size () - 1
	    || sequence_mask[j->first-1] != sequence_mask[j->second+1])
	  {
	    j->type = OTHER;
	  }
	else if (helix_mask[j->first-1] == helix_mask[j->second+1])
	  {
	    node i = 0;
	    // Find pairs of the extrimities.
	    if (helices[helix_mask[j->first-1]].front ().first == j->first-1)
	      i = helices[helix_mask[j->first-1]].front ().second;
	    else if (helices[helix_mask[j->first-1]].front ().second == j->first-1)
	      i = helices[helix_mask[j->first-1]].front ().first;
	    else if (helices[helix_mask[j->first-1]].back ().first == j->first-1)
	      i = helices[helix_mask[j->first-1]].back ().second;
	    else if (helices[helix_mask[j->first-1]].back ().second == j->first-1)
	      i = helices[helix_mask[j->first-1]].back ().first;
	    else 
	      cerr << "Pairing not found for (a) " << j->first-1 << endl;

	    if (i==j->second+1)
	      j->type = LOOP;
	    else 
	      j->type = BULGE_OUT;
	  }
	else
	  { 
	    node i = 0, k = 0;
	
	    // Find pairs.
	    if (helices[helix_mask[j->first-1]].front ().first == j->first-1)
	      i = helices[helix_mask[j->first-1]].front ().second;
	    else if (helices[helix_mask[j->first-1]].front ().second == j->first-1)
	      i = helices[helix_mask[j->first-1]].front ().first;
	    else if (helices[helix_mask[j->first-1]].back ().first == j->first-1)
	      i = helices[helix_mask[j->first-1]].back ().second;
	    else if (helices[helix_mask[j->first-1]].back ().second == j->first-1)
	      i = helices[helix_mask[j->first-1]].back ().first;
	    else 
	      cerr << "Pairing not found for (b) " << j->first-1 << endl;

	    if (helices[helix_mask[j->second+1]].front ().first == j->second+1)
	      k = helices[helix_mask[j->second+1]].front ().second;
	    else if (helices[helix_mask[j->second+1]].front ().second == j->second+1)
	      k = helices[helix_mask[j->second+1]].front ().first;
	    else if (helices[helix_mask[j->second+1]].back ().first == j->second+1)
	      k = helices[helix_mask[j->second+1]].back ().second;
	    else if (helices[helix_mask[j->second+1]].back ().second == j->second+1)
	      k = helices[helix_mask[j->second+1]].back ().first;
	    else 
	      cerr << "Pairing not found for (c) " << j->first-1 << endl;

	    if (sequence_mask[i] != sequence_mask[k])
	      {
		j->type = OTHER;
	      }
	    else if (i==k+1 || i==k-1)
	      {
		j->type = BULGE;
	      }
	    else
	      {
		// Find the other part of the internal loop
		vector< Strand >::iterator s;
		for (s=strands.begin (); s!=strands.end (); ++s)
		  {
		    if (s->first == k+1 && s->second == i-1)
		      {
			j->type = INTERNAL_LOOP;
			j->ref = s-strands.begin ();
			break;
		      }
		    else
		      {
			j->type = OTHER;
		      }
		  }
	      }
	  }    
      }
  }

  
  void
  AnnotateModel::findKissingHairpins ()
  {
    UndirectedGraph< node, edge >::iterator gi;
    list< node >::iterator gj;
    list< node > neigh;
    int check = 0;
    int nb_helical_bp = 0;
    set< const PropertyType* >::iterator k;

    cout.setf (ios::left, ios::adjustfield);

    for (gi=graph.begin (); gi!=graph.end (); ++gi)
      {
	neigh = graph.getNeighbors (*gi);
	for (gj=neigh.begin (); gj!=neigh.end (); ++gj)
	  {
	    if (*gi < *gj && isPairing (graph.getEdge (*gi, *gj)))
	      {
		if (helix_mask [*gi] != -1
		    && helix_mask [*gi] == helix_mask [*gj])
		  {
		    // Simple helical base-pair.
		    cout << setw (20) << "helix bp: " 
			 << getResId (*gi) << "-" 
			 << getResId (*gj) << " : ";
		    check++;
		    nb_helical_bp++;
		  } 
	
		if (helix_mask [*gi] != -1
		    && helix_mask [*gi] != helix_mask [*gj])
		  {
		    cout << setw (20) << "inter helix: " 
			 << getResId (*gi) << "-" 
			 << getResId (*gj) << " : ";
		    tertiary_mask[*gi] = 1;
		    tertiary_mask[*gj] = 1;
		    check++;
		  }

		if (strand_mask [*gi] != -1
		    && helix_mask [*gj] != -1)
		  {
		    if (strands[strand_mask [*gi]].type == OTHER)
		      {
			cout << setw (20) << "strand/helix: " 
			     <<  getResId (*gi) << "-" 
			     << getResId (*gj) << " : ";
			tertiary_mask[*gi] = 1;
			tertiary_mask[*gj] = 1;
			check++;
		      }
		    else if (strands[strand_mask [*gi]].type == BULGE
			     || strands[strand_mask [*gi]].type == BULGE_OUT)
		      {
			cout << setw (20) << "bulge/helix: " 
			     <<  getResId (*gi) << "-" 
			     << getResId (*gj) << " : ";
			tertiary_mask[*gi] = 1;
			tertiary_mask[*gj] = 1;
			check++;
		      }
		    else if (strands[strand_mask [*gi]].type == LOOP)
		      {
			cout << setw (20) << "loop/helix: " 
			     <<  getResId (*gi) << "-" 
			     << getResId (*gj) << " : ";
			tertiary_mask[*gi] = 1;
			tertiary_mask[*gj] = 1;
			check++;
		      }
		    else if (strands[strand_mask [*gi]].type == INTERNAL_LOOP)
		      {
			cout << setw (20) << "internal loop/helix: " 
			     <<  getResId (*gi) << "-" 
			     << getResId (*gj) << " : ";
			tertiary_mask[*gi] = 1;
			tertiary_mask[*gj] = 1;
			check++;
		      }
		    else
		      {
			cout << "Other A: " << " : ";
		      }
		  }
	
		if (strand_mask [*gi] != -1
		    && strand_mask [*gj] != -1)
		  {
	  
		    if (strand_mask [*gi] == strand_mask [*gj])
		      {
			cout << setw (20) << "intraloop: " 
			     << getResId (*gi) << "-"
			     << getResId (*gj) << " : ";
			tertiary_mask[*gi] = 1;
			tertiary_mask[*gj] = 1;
			check++;
		      }
		    else
		      {
			int size = 0;
			if (strands[strand_mask [*gi]].type == OTHER)
			  {
			    cout << "strand/";
			    size = 20-7;	    
			  }
			else if (strands[strand_mask [*gi]].type == BULGE
				 || strands[strand_mask [*gi]].type == BULGE_OUT)
			  {
			    cout << "bulge/";
			    size = 20-6;
			  }
			else if (strands[strand_mask [*gi]].type == LOOP)
			  {
			    cout << "loop/";
			    size = 20-5;
			  }
			else if (strands[strand_mask [*gi]].type == INTERNAL_LOOP)
			  {
			    cout << "internal loop/";
			    size = 20-14;
			  }
			if (strands[strand_mask [*gj]].type == OTHER)
			  {
			    cout << setw (size) << "strand: " 
				 <<  getResId (*gi) << "-" 
				 << getResId (*gj) << " : ";
			    tertiary_mask[*gi] = 1;
			    tertiary_mask[*gj] = 1;
			    check++;
			  }
			else if (strands[strand_mask [*gj]].type == BULGE
				 || strands[strand_mask [*gj]].type == BULGE_OUT)
			  {
			    cout << setw (size) << "bulge: " 
				 <<  getResId (*gi) << "-" 
				 << getResId (*gj) << " : ";
			    tertiary_mask[*gi] = 1;
			    tertiary_mask[*gj] = 1;
			    check++;
			  }
			else if (strands[strand_mask [*gj]].type == LOOP)
			  {
			    cout << setw (size) << "loop: " 
				 <<  getResId (*gi) << "-" 
				 << getResId (*gj) << " : ";
			    tertiary_mask[*gi] = 1;
			    tertiary_mask[*gj] = 1;
			    check++;
			  }
			else if (strands[strand_mask [*gj]].type == INTERNAL_LOOP)
			  {
			    cout << setw (size) << "internal loop: " 
				 <<  getResId (*gi) << "-" 
				 << getResId (*gj) << " : ";
			    tertiary_mask[*gi] = 1;
			    tertiary_mask[*gj] = 1;
			    check++;
			  }
			else
			  {
			    cout << "Other B:" << " : ";
			  }
		      }
		  }
	
		edge e = graph.getEdge (*gi, *gj);
	
		cout << Pdbstream::stringifyResidueType (getType (*gi)) << "-" << flush
		     << Pdbstream::stringifyResidueType (getType (*gj)) << " " << flush;
	
		if (relations[e].getRefFace ())
		  cout << relations[e].getRefFace () << "/" << flush
		       << relations[e].getResFace () << " " << flush;
	
		for (k=relations[e].getLabels ().begin (); 
		     k!=relations[e].getLabels ().end (); ++k)
		  {
		    cout << **k << " " ;
		  }
		cout << endl;       
	      }
	  }
      }

    cout << "Number of base pairs = " << nb_pairings << endl;
    cout << "Number of helical base pairs = " << nb_helical_bp << endl;

    if (nb_pairings != check) 
      cout << "Missing interactions: " << check << "/" << nb_pairings << endl;
  }


  void
  AnnotateModel::findPseudoknots ()
  {
    vector< Helix >::iterator i, j;
  
    for (i=helices.begin (); i!=helices.end (); ++i)
      {
	if (strand_mask[i->front ().first] == strand_mask[i->front ().second])
	  {
	    node a, b, c, d;
      
	    a = i->front ().first;
	    b = i->back ().first;
	    c = i->back ().second;
	    d = i->front ().second;
	    for (j=i; j!=helices.end (); ++j)
	      {
		if (i!=j
		    && strand_mask[j->front ().first] ==  strand_mask[j->front ().second])
		  {
		    node ap, bp, cp, dp;
	  
		    ap = j->front ().first;
		    bp = j->back ().first;
		    cp = j->back ().second;
		    dp = j->front ().second;
	  
		    if (ap > b && bp < c && cp > d)
		      {
			cout << "Pseudoknot : " << "H" << helix_mask[a]+1 << " " 
			     << "H" << helix_mask[ap]+1 << endl;
		      }
		  }
	      }
	  }
      }
  }


  ResIdSet
  AnnotateModel::extract (ResIdSet &seed, int size)
  {
    set< int > query;
    set< int > result;
    set< int > visited;
    set< int >::iterator s;
    ResIdSet::iterator r;
    map< ResId, ResId >::iterator i;
    set< int > tmp;
    ResIdSet outset;
  
    for (r = seed.begin (); r != seed.end (); ++r)
      {
	for (i = translation.begin (); i != translation.end (); ++i)
	  {
	    if (i->second == *r) break;
	  }
    
	if (i != translation.end ())
	  {
	    query.insert (i->first.getResNo ());
	    visited.insert (i->first.getResNo ());
	  }
      }

    while (size-- > 0)
      {
	result.clear ();
	for (s = query.begin (); s != query.end (); ++s)
	  {
	    int id = *s;
	    list< node > neigh;
	    list< node >::iterator i;
      
	    neigh = graph.getNeighbors (id);
	    for (i = neigh.begin (); i != neigh.end (); ++i)
	      {
		result.insert (*i);
	      }
	  }

	set_difference (result.begin (), result.end (),
			visited.begin (), visited.end (),
			inserter (tmp, tmp.begin ()));
	query = tmp;
	tmp.clear ();
    
	set_union (visited.begin (), visited.end (),
		   result.begin (), result.end (),
		   inserter (tmp, tmp.begin ()));
	visited = tmp;
      }
  
    for (s=visited.begin (); s!=visited.end (); ++s)
      outset.insert (getResId (*s));
  
    cout << outset << endl;

    return outset;
  }



// void
// AnnotateModel::decompose ()
// {
//   list< Block > s_blocks_tmp;
//   vector< Block >::iterator bi;
//   list< Block >::iterator bk, bl;
  
//   {
//     // Find helices...
//     vector< Helix >::iterator i;
//     Helix::iterator j;
//     Helix::reverse_iterator rj;
    
//     for (i=helices.begin (); i!=helices.end (); ++i) {
//       Block b (this);
//       vector< node > v;
//       vector< node > vi;
//       for (j=i->begin (); j!=i->end (); ++j) {
// 	v.push_back (j->first);
//       }
//       for (rj=i->rbegin (); rj!=i->rend (); ++rj) {
// 	vi.push_back (rj->second);
//       }    
//       b.push_back (v);    
//       b.push_back (vi);
//       b.junctions.push_back (make_pair (v.front (), vi.back ()));
//       b.junctions.push_back (make_pair (vi.front (), v.back ()));
      
//       h_blocks.push_back (b);
//     }
//   }
//   {
//     // Find strands...
//     vector< Strand >::iterator i;
    
//     for (i=strands.begin (); i!=strands.end (); ++i) {
//       Block b (this);
//       vector< node > v;

//       if (i->first>0 && sequence_mask[i->first] == sequence_mask[i->first-1]) {
// 	v.push_back (i->first-1);
// 	strand_mask[i->first-1] = strand_mask[i->first];
//       }
//       for (int k=i->first; k<=i->second; ++k) {
// 	v.push_back (k);
//       }
//       if (i->second<(int)sequence_mask.size ()-1 && sequence_mask[i->second] == sequence_mask[i->second+1]) {
// 	v.push_back (i->second+1);
// 	strand_mask[i->second+1] = strand_mask[i->second];
//       }
//       b.push_back (v);
      
//       s_blocks_tmp.push_back (b);
//     }
    
//     // Hack: Find regions where helices are connected with no strand...
    
//     for (int j=0; j<(int)conformations.size ()-1; ++j) {
//       if (helix_mask[j] != -1 && helix_mask[j+1] != -1 &&
// 	  helix_mask[j] != helix_mask[j+1] && 
// 	  sequence_mask[j] == sequence_mask[j+1]) {
// 	Block b (this);
// 	vector< node > v;
// 	v.push_back (j);
// 	v.push_back (j+1);
// 	b.push_back (v);
// 	b.output (cout);
// 	cout << endl;
// 	s_blocks_tmp.push_back (b);
// 	strand_mask[j] = s_blocks_tmp.size ()-1;
// 	strand_mask[j+1] = s_blocks_tmp.size ()-1;
//       }
//     }
//   }
  
//   {
//     // Find loops and internal loops...
//     bk = s_blocks_tmp.begin ();
//     while (bk!=s_blocks_tmp.end ()) {
//       node s, e;
//       int hs, he;
      
      
//       bk->output (cout); cout << endl;

//       s = bk->begin ()->front ();
//       e = bk->begin ()->back ();
      
//       hs = helix_mask[s];
//       he = helix_mask[e];

//       if (hs==he) {
// 	bk->junctions.push_back (make_pair (s, e));
// 	l_blocks.push_back (*bk);
// 	bk = s_blocks_tmp.erase (bk);
//       } else if (hs ==-1 || he == -1) {
// 	s_blocks.push_back (*bk);
// 	bk = s_blocks_tmp.erase (bk);
//       } else {
// 	node ep = s;
// 	node sp = s;

// 	while (ep!=e) {
// 	  ep = 
// 	    (h_blocks[hs].junctions.front ().first == s) ? h_blocks[hs].junctions.front ().second :
// 	    ((h_blocks[hs].junctions.front ().second == s) ? h_blocks[hs].junctions.front ().first :
// 	     ((h_blocks[hs].junctions.back ().first == s) ? h_blocks[hs].junctions.back ().second : 
// 	      h_blocks[hs].junctions.back ().first));
// 	  // Find new s.
// 	  if (ep!=e) {
// 	    bl = s_blocks_tmp.begin ();
// 	    while (bl!=s_blocks_tmp.end ()) {
	      
// 	      if (bl->front ().front () == ep) { 
// 		sp = bl->front ().back ();
// 		break;
// 	      }
// 	      else if (bl->front ().back () == ep) {
// 		sp = bl->front ().front ();
// 		break;
// 	      }
// 	      bl++;
// 	    }
	    
// 	    if (bl==s_blocks_tmp.end ()) cerr << "Woa!  This should not happen" << endl;

// 	    // Add bl to bk and erase bl
// 	    bk->push_back (bl->front ());
// 	    s_blocks_tmp.erase (bl);
// 	  }
	  
// 	  bk->junctions.push_back (make_pair (s, ep));
// 	  s = sp;
// 	  hs = helix_mask[s];
// 	}
	
// 	l_blocks.push_back (*bk);
// 	bk = s_blocks_tmp.erase (bk);
//       }
//     }
//   }
  
//   cout << "Helices" << endl;
//   for (bi=h_blocks.begin (); bi!=h_blocks.end (); ++bi) {
//     bi->output (cout);
//     cout << endl;
//   }

//   cout << "Strand" << endl;
//   for (bk=s_blocks_tmp.begin (); bk!=s_blocks_tmp.end (); ++bk) {
//     bk->output (cout);
//     cout << endl;
//   }
//   cout << "Loops" << endl;
//   for (bi=l_blocks.begin (); bi!=l_blocks.end (); ++bi) {
//     bi->output (cout);
//     cout << endl;
//   }
// }


  void
  AnnotateModel::dumpSequences (bool detailed) 
  {
    int i, j;
    int currseq = -1;
  
    if (!detailed)
      {
	for (i = 0; i < (int) conformations.size (); ++i)
	  {
	    if (i == 0 || sequence_mask[i] != sequence_mask[i-1])
	      cout << " " << setw (8) << getResId (i) << " " << flush;
	    cout << getType (i)->toString () << flush;
	  }
	return;
      } 

    for (i = 0; i < (int) conformations.size (); )
      {
	int pos = 0;

	currseq = sequence_mask[i];
	cout << "Sequence " << currseq << " (length = " 
	     << sequence_length[currseq] << "): " << endl;
    
	while (pos < sequence_length[currseq]) 
	  {
	    for (j = pos; (j < sequence_length[currseq] && (j+1)%50!=0); ++j)
	      {
		if ((j) % 50 == 0)
		  {
		    cout.setf (ios::right, ios::adjustfield);
		    cout << setw (3) << getResId (i+j).getChainId () << flush;
		    cout.setf (ios::left, ios::adjustfield);
		    cout << setw (5) << getResId (i+j).getResNo () << " " << flush;
		  }
	  
		cout << (getType (i + j)->isNucleicAcid ()
			 ? Pdbstream::stringifyResidueType (getType (i+j))
			 : "X")
		     << flush;
	  
		if ((j + 1) % 10 == 0)
		  cout << " " << flush;
	      }
	    cout << endl;
	
	    for (j = pos; (j < sequence_length[currseq] && (j + 1) % 50 != 0); ++j)
	      {
		if ((j) % 50 == 0)
		  cout << setw (8) << " " << " " << flush;
		cout << marks[i + j] << flush;
		if ((j + 1) % 10 == 0)
		  cout << " " << flush;
	      }
	    cout << endl;

//  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
//  	  if ((j)%50==0) cout << setw (6) << "seq" << " " << flush;
//  	  if (sequence_mask[i+j] == -1) cout << "-" << flush;
//  	  else cout << sequence_mask[i+j] << flush;
//  	  if ((j+1)%10==0) cout << " " << flush;
//  	}
//  	cout << endl;

	    for (j = pos; (j < sequence_length[currseq] && (j + 1) % 50 != 0); ++j)
	      {
		if ((j) % 50 == 0)
		  cout << setw (8) << "helix" << " " << flush;
		if (helix_mask[i + j] == -1)
		  cout << "-" << flush;
		else
		  cout << helix_mask[i + j] << flush;
		if ((j + 1) % 10 == 0)
		  cout << " " << flush;
	      }
	    cout << endl;

//  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
//  	  if ((j)%50==0) cout << setw (6) << "str" << " " << flush;
//  	  if (strand_mask[i+j] == -1) cout << "-" << flush;
//  	  else cout << strand_mask[i+j] << flush;
//  	  if ((j+1)%10==0) cout << " " << flush;
//  	}
//  	cout << endl;

// 	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
// 	  if ((j)%50==0) cout << setw (6) << "ter" << " " << flush;
// 	  if (tertiary_mask[i+j] == -1) cout << "-" << flush;
// 	  else cout << 'X' << flush;
// 	  if ((j+1)%10==0) cout << " " << flush;
// 	}

// 	cout << endl << endl;

	    pos = j + 1;
	  }
	i += pos - 1;
      }
    cout.setf (ios::left, ios::adjustfield);
  }


// void AnnotateModel::dumpGraph () 
// {
//   Graph< node, edge >::adjgraph::iterator i;
//   Graph< node, edge >::adjlist::iterator j;
//   PropertySet::iterator k;

//   cout << "Graph of relations ----------------------------------------------" << endl;

//   for (i=graph.begin (); i!=graph.end (); ++i) {
//     for (j=i->second.begin (); j!=i->second.end (); ++j) {
//       cout << getResId (i->first) << "-" 
// 	   << getResId (j->first) << " : ";
//       for (k=relations[j->second].getGlobalProp ().begin (); 
// 	   k!=relations[j->second].getGlobalProp ().end (); ++k) {
// 	cout << (const char*)**k << " " ;
//       }
//       cout << endl;
//     }
//   } 
// }


  void
  AnnotateModel::dumpConformations () 
  {
    UndirectedGraph< node, edge >::iterator i;
  
    for (i = graph.begin (); i != graph.end (); ++i)
      {
	cout << getResId (*i) << " : " 
	     << (const char*) *getType (*i) << " " 
	     << conformations[*i]->getPucker () << " " 
	     << conformations[*i]->getGlycosyl () << endl; 
      }
  }


  void
  AnnotateModel::dumpPairs () 
  {
    UndirectedGraph< node, edge >::iterator i;
    list< node >::iterator j;
    list< node > neigh;
    set< const PropertyType* >::iterator k;

    for (i = graph.begin (); i != graph.end (); ++i)
      {

	neigh = graph.getNeighbors (*i);

	for (j = neigh.begin (); j != neigh.end (); ++j)
	  {
	    if (*i < *j && isPairing (*i, *j))
	      {
		cout << getResId (*i) << "-" 
		     << getResId (*j) << " : ";

		edge e = graph.getEdge (*i, *j);

		cout << Pdbstream::stringifyResidueType (getType (*i)) << "-"
		     << Pdbstream::stringifyResidueType (getType (*j)) << " ";
		if (relations[e].getRefFace ())
		  cout << relations[e].getRefFace () << "/"
		       << relations[e].getResFace () << " ";
	
		for (k = relations[e].getLabels ().begin (); 
		     k!=relations[e].getLabels ().end (); ++k)
		  {
		    cout << **k << " " ;
		  }
		cout << endl;
	      }
	  } 
      }
  }


  void
  AnnotateModel::dumpStacks () 
  {
    UndirectedGraph< node, edge >::iterator i, k;
    list< node >::iterator j;
    list< node > neigh;
    set< const PropertyType* >::iterator l;

    int nb_stacks = 0;
    int nb_helical_stacks = 0;
    int nb_adjacent_stacks = 0;

    cout << "Adjacent stackings ----------------------------------------------" << endl;

    for (i = graph.begin (), k = i, ++k; i != graph.end () && k != graph.end (); ++i, ++k)
      {
	neigh = graph.getNeighbors (*i);

	if ((j = std::find (neigh.begin (), neigh.end (), *k)) != neigh.end ()
	    && relations[graph.getEdge (*i, *j)].is (PropertyType::pAdjacent))
	  {
	    if (relations[graph.getEdge (*i, *j)].is (PropertyType::pStack))
	      {
		nb_adjacent_stacks++;
		nb_stacks++;
		if (helix_mask[*i] == helix_mask[*j] && helix_mask[*i] != -1)
		  nb_helical_stacks++;
	      }
      
	    cout << getResId (*i) << "-" << getResId (*j) << " : ";
	    for (l = relations[graph.getEdge (*i, *j)].getLabels ().begin (); 
		 l != relations[graph.getEdge (*i, *j)].getLabels ().end (); ++l)
	      {
		cout << **l << " " ;
	      }
	    cout << endl;
	  }
      }
    cout << endl;

    cout << "Non-Adjacent stackings ------------------------------------------" << endl;
    for (i = graph.begin (); i != graph.end (); ++i)
      {
	neigh = graph.getNeighbors (*i);
    
	for (j = neigh.begin (); j != neigh.end (); ++j)
	  {
	    if (*i < *j
		&& relations[graph.getEdge (*i, *j)].is (PropertyType::pStack)
		&& ! relations[graph.getEdge (*i, *j)].is (PropertyType::pAdjacent))
	      {
		nb_stacks++;
		cout << getResId (*i) << "-" << getResId (*j) << " : ";
		for (l = relations[graph.getEdge (*i, *j)].getLabels ().begin (); 
		     l != relations[graph.getEdge (*i, *j)].getLabels ().end (); ++l)
		  {
		    //if (*l != PropertyType::pReverse)
		    cout << (const char*)**l << " ";
		  }
		cout << endl;
	      }
	  }
      }

    cout << "Number of stackings = " << nb_stacks << endl
	 << "Number of helical stackings = " << nb_helical_stacks << endl
	 << "Number of adjacent stackings = " << nb_adjacent_stacks << endl
	 << "Number of non adjacent stackings = "
	 << nb_stacks - nb_adjacent_stacks << endl;
  }


  void
  AnnotateModel::dumpHelices () 
  {
    vector< Helix >::iterator i;
    //   Helix::iterator j;

     
    for (i = helices.begin (); i != helices.end (); ++i)
      {
//     char type;
//     int tmpsize = 0;
//     if (conformations[i->begin ()->first]->getPucker ()->is (PropertyType::pC3p_endo) &&
// 	conformations[i->begin ()->first]->getGlycosyl ()->is (PropertyType::pAnti))
//       type = 'A';
//     else if (conformations[i->begin ()->first]->getPucker ()->is (PropertyType::pC2p_endo) &&
// 	     conformations[i->begin ()->first]->getGlycosyl ()->is (PropertyType::pAnti))
//       type = 'B';
//     else type = 'X';
    
//     for (j=i->begin (); j!=i->end (); ++j) {
//       if (j->first >= 0) {
// 	char tmptype = 'X';
// 	if (conformations[j->first]->getPucker ()->is (PropertyType::pC3p_endo) &&
// 	    conformations[j->first]->getGlycosyl ()->is (PropertyType::pAnti))
// 	  tmptype = 'A';
// 	else if (conformations[j->first]->getPucker ()->is (PropertyType::pC2p_endo) &&
// 		 conformations[j->first]->getGlycosyl ()->is (PropertyType::pAnti))
// 	  tmptype = 'B';
// 	else tmptype = 'X';
	
// 	if (type != tmptype) type = 'X';
	
// 	if (conformations[j->second]->getPucker ()->is (PropertyType::pC3p_endo) &&
// 	    conformations[j->second]->getGlycosyl ()->is (PropertyType::pAnti))
// 	  tmptype = 'A';
// 	else if (conformations[j->second]->getPucker ()->is (PropertyType::pC2p_endo) &&
// 		 conformations[j->second]->getGlycosyl ()->is (PropertyType::pAnti))
// 	  tmptype = 'B';
// 	else tmptype = 'X';
	
// 	if (type != tmptype) type = 'X';
	
// 	tmpsize++;
//       }
//     }
    

	cout << "H" << i-helices.begin () << ", length = " 
	     << i->size () << "" << endl;
    
//     cout << "type = " << type << " : " << endl;
    
	int k;


	cout.setf (ios::right, ios::adjustfield);
	cout << setw (6) << getResId ((*i)[0].first) << "-" << flush;
	for (k = 0; k < (int) i->size (); ++k)
	  {
	    if ((*i)[k].first >= 0)
	      {
		cout << (getType ((*i)[k].first)->isNucleicAcid ()
			 ? Pdbstream::stringifyResidueType (getType ((*i)[k].first))
			 : "X")
		     << flush;
	      }
	    else
	      {
		if (- (*i)[k].first -1 == 0)
		  {
		    cout << "-" << setw (max ((int) floor (log10 ((double)-(*i)[k].first)),
					      (int) floor (log10 ((double)-(*i)[k].second))) + 1) 
			 << setfill ('-') << "-" << "-" << flush;
		  }
		else
		  {
		    cout << "-" << setw (max ((int) floor (log10 ((double)-(*i)[k].first)),
					      (int) floor (log10 ((double)-(*i)[k].second))) + 1) 
			 << - (*i)[k].first -1 << "-" << flush;
		  }
	      }
	  }
	cout << "-" << getResId ((*i)[i->size ()-1].first) << endl;
    
	cout.setf (ios::right, ios::adjustfield);
	cout << setw (6) << getResId ((*i)[0].second) << "-" << flush;

	for (k = 0; k < (int) i->size (); ++k)
	  {
	    if ((*i)[k].first >= 0)
	      {
		cout << (getType ((*i)[k].second)->isNucleicAcid ()
			 ? Pdbstream::stringifyResidueType (getType ((*i)[k].second))
			 : "X")
		     << flush;
	      }
	    else
	      {
		if (-(*i)[k].second -1 == 0)
		  {
		    cout << "-" << setw (max ((int) floor (log10 ((double)-(*i)[k].first)),
					      (int) floor (log10 ((double)-(*i)[k].second))) + 1) 
			 << setfill ('-') << "-" << "-" << flush;
		  }
		else
		  {  
		    cout << "-" << setw (max ((int) floor (log10 ((double)-(*i)[k].first)),
					      (int) floor (log10 ((double)-(*i)[k].second))) + 1) 
			 << - (*i)[k].second -1 << "-" << flush;
		  }
	      }
	  }
	cout << "-" << getResId ((*i)[i->size ()-1].second)
	     << setfill (' ') << endl;
      }
    cout.setf (ios::left, ios::adjustfield);
  }


  void
  AnnotateModel::dumpTriples () 
  {
    UndirectedGraph< node, edge >::iterator gi;
    list< node >::iterator gj;
    list< node > neigh;
    int count;
    int nb = 1;
    set< node > nodes;
    set< node >::iterator k, l;
    vector< bool > treated;

    treated.resize (graph.size (), false);

    for (gi = graph.begin (); gi != graph.end (); ++gi)
      {
	count = 0;
	nodes.clear ();
	if (! treated[*gi])
	  {
	    nodes.insert (*gi);
      
	    neigh = graph.getNeighbors (*gi);

	    for (gj = neigh.begin (); gj != neigh.end (); ++gj)
	      {
		if (isPairing (graph.getEdge (*gi,*gj)))
		  {
		    count++;
		    nodes.insert (*gj);
		  }
	      }
	    if (count > 1)
	      {
		cout.setf (ios::left, ios::adjustfield);
		cout << "T" << setw (5) << nb++ << " "; 
		for (k = nodes.begin (); k != nodes.end (); ++k)
		  {
		    for (l = k; l != nodes.end (); ++l)
		      {
			if (k != l && isPairing (*k, *l))
			  {
			    cout << getResId (*k) << "-"
				 << getResId (*l) << " ";
			    treated[*k] = true;
			    treated[*l] = true;
			  }
		      }
		  }
		cout << endl;
	      }
	  }
      }
    cout.setf (ios::left, ios::adjustfield);
  }

  
  void 
  AnnotateModel::dumpStrands () 
  {
    int j;
    int k;
  
    for (j = 0; j < (int) strands.size (); ++j)
      {
	cout.setf (ios::left, ios::adjustfield);
	cout << "S" << setw (5) << j << " " <<  flush;
    
	cout.setf (ios::left, ios::adjustfield);
    
	if (strands[j].type == OTHER)
	  cout << setw (16) << "single strand:";
	else if (strands[j].type == BULGE)
	  cout << setw (16) << "bulge:";
	else if (strands[j].type == BULGE_OUT)
	  cout << setw (16) << "bulge out:";
	else if (strands[j].type == LOOP)
	  cout << setw (16) << "loop:";
	else if (strands[j].type == INTERNAL_LOOP)
	  cout << setw (16) << "internal loop:";
    
	cout << getResId (strands[j].first) << "-" << flush;
	for (k = strands[j].first; k <= strands[j].second; ++k)
	  {
	    cout << (getType (k)->isNucleicAcid ()
		     ? Pdbstream::stringifyResidueType (getType (k))
		     : "X")
		 << flush;
	  }
	cout << "-" << getResId (strands[j].second);
    
	if (strands[j].type == INTERNAL_LOOP)
	  {
	    cout << " -- " << flush;
	    cout << getResId (strands[strands[j].ref].first) << "-" << flush;
	    for (k = strands[strands[j].ref].first; k <= strands[strands[j].ref].second; ++k)
	      {
		cout << (getType (k)->isNucleicAcid ()
			 ? Pdbstream::stringifyResidueType (getType (k))
			 : "X")
		     << flush;
	      }
	    cout << "-" << getResId (strands[strands[j].ref].second);
	  }
	cout << endl;
      }
  
    //   for (k=strands[j].first; k<=strands[j].second; ++k) {
    //     cout << *getType (k) << flush;
    //   }
    //   cout << " " << flush;
    // }
    cout.setf (ios::left, ios::adjustfield);
  }


  void 
  AnnotateModel::dumpMcc (const char* pdbname)
  {
    char str[256];
  
    cout << "//" << endl;
    cout << "// Annotation results ------------------------------------------------" << endl;
    cout << "//" << endl;
    gethostname (str, 255);
    cout << "// Author         : " << getenv ("LOGNAME") << '@' << str << endl;
    cout << "// Structure file : " << pdbname << endl;
    cout << "// Structural annotation generated by " << PACKAGE << " "  << VERSION << endl;
    cout << "// ";
    cout << endl;
  
    cout << "// Sequences ---------------------------------------------------------" << endl;
    cout << "// The distance between atoms O3' and P of two residues must be       " << endl;
    cout << "// inferior to 20nm for them to be considered adjacent.               " << endl;
    cout << "//" << endl;
    {
      int i, j;
      int currseq = -1;
    
      for (i = 0; i < (int) conformations.size (); )
	{
	  int pos = 0;

	  currseq = sequence_mask[i];
	  cout << "sequence( r " << endl;
	  while (pos < sequence_length[currseq]) 
	    {
	      for (j = pos; (j < sequence_length[currseq] && (j + 1) % 50 !=0 ); ++j)
		{
		  if (j == 0)
		    cout << setw (8) << getResId (i+j) << " ";
		  cout << getType (i+j)->toString ();
		  if ((j + 1) % 10 == 0) cout << " ";
		}
	      cout << endl;
	      pos = j + 1;
	    }
      
	  cout << ")" << endl;
      
	  i += pos - 1;
	}
    }

    cout << "//" << endl;
    cout << "// Nucleotide conformations ------------------------------------------" << endl;
    cout << "// The distance between atoms O3' and P of two residues must be       " << endl;
    cout << "// inferior to 20nm for them to be considered adjacent.               " << endl;
    cout << "//" << endl;
    {
      UndirectedGraph< node, edge >::iterator i;

      cout << "residue(" << endl;
      for (i = graph.begin (); i != graph.end (); ++i)
	{
	  cout << setw(8) << getResId (*i) << " { ";
	  cout.setf (ios::left, ios::adjustfield);
	  cout << setw(8) << *conformations[*i]->getPucker () << " && " 
	       << setw(4) << *conformations[*i]->getGlycosyl () << " }"
	       << "  100%" << endl; 
	  cout.setf (ios::right, ios::adjustfield);
	}
      cout << ")" << endl;
    }

    cout << "//" << endl;
    cout << "// Non-adjacent base-pairs and stackings ------------------------" << endl;
    cout << "//" << endl;
 
    {
      UndirectedGraph< node, edge >::iterator i;
      list< node >::iterator j;
      list< node > neigh;
      set< const PropertyType* >::iterator k;
    
      if (nb_pairings > 0)
	{
	  cout << "pair(" << endl;
	  for (i = graph.begin (); i != graph.end (); ++i)
	    {
	      neigh = graph.getNeighbors (*i);
	      for (j = neigh.begin (); j != neigh.end (); ++j)
		{
		  edge e = graph.getEdge (*i, *j);
	  
		  if (*i < *j && ! isAdjacent (e))
		    {
		      cout << setw(8) << getResId (*i) 
			   << setw(8) << getResId (*j) << " { ";
	    
		      bool useand = false;
		      if (isPairing (e)
			  && ! relations[e].is (PropertyType::parseType ("unclassified")))
			{
			  cout << (const char*)(*relations[e].getRefFace ()) << "/"
			       << (const char*)(*relations[e].getResFace ()) << " ";
			  useand = true;
			}
		      for (k = relations[e].getLabels ().begin (); 
			   k != relations[e].getLabels ().end (); ++k)
			{
			  if ((*k) == PropertyType::pReverse
			      || (*k) == PropertyType::pCis
			      || (*k) == PropertyType::pTrans
			      || (*k) == PropertyType::pStack)
			    {
			      if (useand)
				cout << "&& ";
			      cout << (const char*) **k << " " ;
			      useand = true;
			    }
			}
		      cout << "}" << "  100%" << endl; 
		    }
		}
	    }
	  cout << ")" << endl;
	}
    }
  
    cout << "//" << endl;
    cout << "// Adjacent relations -------------------------------------------" << endl;
    cout << "//" << endl;

    {
      UndirectedGraph< node, edge >::iterator i;
      list< node >::iterator j;
      list< node > neigh;
      set< const PropertyType* >::iterator k;
  
      if (nb_connect > 0)
	{
	  cout << "connect(" << endl;
	  for (i = graph.begin (); i != graph.end (); ++i)
	    {
	      neigh = graph.getNeighbors (*i);
	      for (j = neigh.begin (); j != neigh.end (); ++j)
		{
		  edge e = graph.getEdge (*i, *j);
	  
		  if (*i < *j && isAdjacent (e))
		    {
		      cout << setw(8) << getResId (*i) 
			   << setw(8) << getResId (*j) << " { ";
	    
		      if (isPairing (e))
			{
			  cout << (const char*)(*relations[e].getRefFace ()) << "/"
			       << (const char*)(*relations[e].getResFace ())
			       << " && ";
			}
	    
		      if (isStacking (e))
			cout << "stack" << " ";
		      else 
			cout << "!stack" << " ";
	    
		      for (k = relations[e].getLabels ().begin (); 
			   k != relations[e].getLabels ().end (); ++k)
			{
			  if ((*k) == PropertyType::pReverse
			      || (*k) == PropertyType::pCis
			      || (*k) == PropertyType::pTrans
			      || (*k) == PropertyType::pPairing)
			    {
			      cout << "&& " << (const char*)**k << " " ;
			    }
			}
		      cout << "}" << "  100%" << endl; 
		    }
		} 
	    }
	  cout << ")" << endl;    
	}
    }
  
    cout << "//" << endl;
    cout << "// Construction order -------------------------------------------" << endl;
    cout << "// This section defines a spanning tree of minimal height that   " << endl;
    cout << "// connects all residues of the molecule and states the order in " << endl;
    cout << "// which they are placed in the modeling process."                 << endl;
    cout << "//" << endl;

  // Here, we need to find connected components of the initial graph and build a backtrack
  // for each of them.  The generated script is not guaranteed to work when an annotated
  // PDB file contains many fragments of RNA...  What should we do???????

    {
      cout << "global = backtrack(" << endl;

      UndirectedGraph< node, edge >::iterator i;
      list< node >::iterator j;
      list< node > neigh;
    
      for (i = graph.begin (); i != graph.end (); ++i)
	{
	  neigh = graph.getNeighbors (*i);
	  for (j = neigh.begin (); j != neigh.end (); ++j)
	    {
	      edge e = graph.getEdge (*i, *j);
	
	      if (isPairing (e)
		  && ! relations[e].is (PropertyType::parseType ("unclassified")))
		graph.setEdgeWeight (*i, *j, 1);
	      else
		graph.setEdgeWeight (*i, *j, 2);
	    }
	}

      vector< pair< node, node > > edges;
      vector< pair< node, node > >::reverse_iterator k;
      list< list< node > > treated;
      list< list< node > >::iterator x;
      list< node >::iterator y;

      edges = graph.minimumSpanningTree ();
      for (k = edges.rbegin (); k != edges.rend (); ++k)
	{
	  node a, b;
      
	  a = *graph.find (k->first);
	  b = *graph.find (k->second);

	  bool done = false;
	  for (x = treated.begin (); x != treated.end (); ++x)
	    {
	      if (x->front () == b)
		{
		  x->push_front (a);
		  done = true;
		}
	      else if (x->back () == a)
		{
		  x->push_back (b);
		  done = true;
		}
	      if (done)
		break;
	    }
	  if (!done)
	    {
	      list< node > tmp;
	
	      tmp.push_back (a);
	      tmp.push_back (b);
	      treated.push_front (tmp);
	    }
	}

      // correction pass to keep a valid backtrack statement
      size_t placed_sz = 0;
      set< int > placed;

      for (y = treated.begin ()->begin (); y != treated.begin ()->end (); ++y)
	placed.insert (*y);

      if (treated.size () > 1)
	{
	  x = treated.begin ();
	  x++;
	  for (; x != treated.end (); ++x) 
	    {
	      placed_sz = placed.size ();
	      placed.insert (*x->begin ());

	      if (placed.size () > placed_sz)
		{
		  // oups! reference residue not placed yet!
		  // What if we just flip the sub-list over...
		  placed_sz = placed.size ();
		  placed.insert (x->back ());

		  if (placed.size () > placed_sz)
		    {
		      // hum...something is really wrong here!
		      cerr << "Fatal Error: unable to build a valid backtrack statement" << endl;
		      exit (EXIT_FAILURE);
		    }

		  // update placed residues set, then flip the sub-list over
		  list< int > tmp;
		  for (y = x->begin (); y != x->end (); ++y)
		    {
		      placed.insert (*y);
		      tmp.push_front (*y);
		    }
		  x->clear ();
		  for (y = tmp.begin (); y != tmp.end (); ++y)
		    x->push_back (*y);
		}
	      else
		{
		  // just update placed residues set
		  for (y = x->begin (); y != x->end (); ++y)
		    placed.insert (*y);
		}
	    }
	}
    
      for (x = treated.begin (); x != treated.end (); ++x)
	{
	  cout << "  ( ";
	  for (y = x->begin (); y != x->end (); ++y)
	    {
	      cout << getResId (*y) << " ";
	    }
	  cout << ")" << endl;
	}
      cout << ")" << endl;
    }
  
    cout << "//" << endl;
    cout << "// Molecule cache -----------------------------------------------" << endl;
    cout << "// A cache is normally placed on top of the backtrack so         " << endl;
    cout << "// generated models are kept only if they are different enough   " << endl;
    cout << "// from previously generated models."                              << endl;
    cout << "//" << endl;

    cout << "global_cache = cache(" << endl;
    cout << "  global" << endl;
    cout << "  rmsd (1.0 align base_only  no_hydrogen)" << endl;
    cout << ")" << endl;
  
    cout << "//" << endl;
    cout << "// Constraints --------------------------------------------------" << endl;
    cout << "//" << endl;
    cout << "adjacency( " << endl
	 << "  global 1.0 5.0" << endl
	 << ")" << endl << endl;
    cout << "res_clash(" << endl
	 << "  global" << endl
	 << "  fixed_distance 1.0" << endl
	 << "  all no_hydrogen" << endl
	 << ")" << endl;
    cout << "//" << endl;
    cout << "// Exploration type ---------------------------------------------" << endl;
    cout << "//" << endl;
    cout << "explore(" << endl
	 << "  global_cache" << endl
	 << "  file_pdb (\"global-%05d.pdb\" zipped)" << endl
	 << ")" << endl;

    cout << "//" << endl
	 << "// --------------------------------------------------------------" << endl
	 << "//" << endl;
  }


// void AnnotateModel::dumpCt (const char* pdbname)
// {
// //fprintf(ofp, "%5d %c   %5d %4d %4d %4d\n",

//   //  char str[256];
//   int i, j;

//   cout << setw (5) << (int)conformations.size () 
//        << " ENERGY = -0.0 " << pdbname << endl;

//   int modelId = sequence_mask[0];

//   for (i=0; i<(int)conformations.size (); ++i) {
//     if (sequence_mask[i] != modelId) return;
//     cout << setw (6) << i+1 << " "
// 	 << *getType (i) << "   " 
// 	 << setw (5) << ((i>0)?i-1:0) << " " 
// 	 << setw (4) << ((i<(int)conformations.size ()-1)?i+2:0) << " ";
//     if ((j=helix_mask[i]) != -1) {

//       for (int k=0; k<(int)helices[j].size (); ++k) {
// 	if (helices[j][k].first == i)
// 	  cout << setw (4) << helices[j][k].second+1 << " ";
// 	else if (helices[j][k].second == i)
// 	  cout << setw (4) <<helices[j][k].first+1 << " ";
//       }
//     }
//     else cout << setw (4) << 0 << " ";
//     cout << setw (4) << i+1 << endl;
//   }


// //    for (i=0; i<(int)conformations.size (); ++i) {
// //      cout << setw (6) << getResId (i).GetResNo () << " "
// //  	 << *getType (i) << "   " 
// //  	 << setw (5) << ((i>0)?getResId (i-1).GetResNo ():0) << " " 
// //  	 << setw (4) << ((i<(int)conformations.size ()-1)?getResId (i+1).GetResNo ():0) << " ";
// //      if ((j=helix_mask[i]) != -1) {

// //        for (int k=0; k<helices[j].size (); ++k) {
// //  	if (helices[j][k].first == i)
// //  	  cout << setw (4) << getResId (helices[j][k].second).GetResNo () << " ";
// //  	else if (helices[j][k].second == i)
// //  	  cout << setw (4) << getResId (helices[j][k].first).GetResNo () << " ";
// //        }
// //      }
// //      else cout << setw (4) << 0 << " ";
// //      cout << setw (4) << getResId (i).GetResNo () << endl;
// //    }

// }


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
