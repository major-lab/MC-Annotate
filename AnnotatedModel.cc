//                              -*- Mode: C++ -*- 
// AnnotatedModel.cc
// Copyright © 2001, 2002 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// Last Modified By : 
// Last Modified On : 
// Update Count     : 0
// Status           : Unknown.
// 



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iomanip.h>
#include <iostream.h>
#include <unistd.h>
#include <vector.h>
#include <list.h>

#include "AnnotatedModel.h"
#include "Graph.h"

#include "mccore/Algo.h"
#include "mcpl/mcpl.h"
#include "mcpl/Relation.h"
#include "mcpl/Conformation.h"
#include "mcpl/PairingPattern.h"
#include "mcpl/PropertyType.h"
#include "mcpl/MaximumFlow.h"
#include "mcpl/PairAnnote.h"


AnnotatedModel::AnnotatedModel (Model *m)
  : model (m)
{
  nb_pairings = 0;
  nb_connect = 0;
}

AnnotatedModel::~AnnotatedModel ()
{
  //delete model;
}


void 
AnnotatedModel::annotate ()
{
  Model::iterator i, j;
  list< edge > seq;
  list< edge >::iterator sj;
  vector< Conformation > tmpconformations;

  map< CResId, CResId > table;
  map< CResId, Model::iterator > corresp;
  
  vector< vector< Model::iterator > > sequences;

  // Standardization of ResIDs
  // This is needed for the ExtractContact_AABB algo to work correctly when
  // there are multiple residues of the same resid.
  { 
    Model::iterator i;
    int a;
    for (i=model->begin (), a=1; i!=model->end (); ++i, ++a) {
      CResId resid (a);
      table[resid] = (CResId&)*i;  
      *i = resid;
      corresp[resid] = i;
    }
  }

  // Extraction of possible relations based on the Axis Aligned Bounding Box method
  vector< pair< Model::iterator, Model::iterator > > possibleContacts;
  vector< pair< Model::iterator, Model::iterator > >::iterator l;
    
  possibleContacts = 
    Algo::ExtractContact_AABB (model->begin (), model->end (), 5.0);

  for (l=possibleContacts.begin (); l!=possibleContacts.end (); ++l)
    {
      i = l->first;
      j = l->second;
      
      Relation rel(&*i, &*j);
      rel.isPairing ();
      rel.isStacking ();
      rel.isAdjacent ();

      if (rel.isAdjacent () || rel.isPairing () || rel.isStacking ()) {
	if (rel.getDirection () == p_DIR_5p) {
	  relations.push_back (rel);
	  seq.push_back (relations.size () - 1);	  	    
	  Relation reli(&*j, &*i);
	  reli.isAdjacent ();
	  reli.isPairing ();
	  reli.isStacking ();
	  relations.push_back (reli);
	} else if (rel.getDirection () == p_DIR_3p) {
	  relations.push_back (rel);
	  Relation reli(&*j, &*i);
	  reli.isAdjacent ();
	  reli.isPairing ();
	  reli.isStacking ();
	  relations.push_back (reli);
	  seq.push_back (relations.size () - 1);
	} else {
	  relations.push_back (rel);
	} 
      }
    }

  // Selection sort of the sequence.
  while (seq.size () > 0) {
    list< edge > sorted_seq;
    sj = seq.begin ();
    sorted_seq.push_back (*sj);
    seq.erase (sj);
    
    while (true) {
      for (sj=seq.begin (); sj!=seq.end (); ++sj) {
	if ((const CResId&)(relations[*sj].getRes ()) == 
	    (const CResId&)(relations[sorted_seq.front ()].getRef ())) {
	  break;
	} 
      }
      if (sj != seq.end ()) {
	sorted_seq.push_front (*sj);
	seq.erase (sj);
      } else {
	break;
      }
    }
    
    while (true) {
      for (sj=seq.begin (); sj!=seq.end (); ++sj) {
	if ((const CResId&)(relations[*sj].getRef ()) == 
	    (const CResId&)(relations[sorted_seq.back ()].getRes ())) {
	  break;
	} 
      }
      if (sj != seq.end ()) {
	sorted_seq.push_back (*sj);
	seq.erase (sj);
      } else {
	break;
      }
    }

    vector< Model::iterator > actual_seq;

    sj=sorted_seq.begin ();
    actual_seq.push_back (corresp[(const CResId&)(relations[*sj].getRef ())]);
    while (sj!=sorted_seq.end ()) {
      actual_seq.push_back (corresp[(const CResId&)(relations[*sj].getRes ())]);
      sj++;
    }

    sequences.push_back (actual_seq);
  }

  // Dump lone residues.
  {
    Model::iterator i;
    for (i=model->begin (); i!=model->end (); ++i) {
      bool found = false;
      for (int j=0; j<(int)sequences.size (); ++j) {
	if (find (sequences[j].begin (), sequences[j].end (), i)
	    != sequences[j].end ()) {
	  found = true;
	  break;
	}
      }
      if (!found) { 
	vector< Model::iterator > actual_seq;
	actual_seq.push_back (i);
  	sequences.push_back (actual_seq);
      }
    }
  }

  map< CResId, node > newcorresp;

  // Create graph...
  {
    vector< vector< Model::iterator > >::iterator i;
    vector< Model::iterator >::iterator j;
    int k = 0;

    for (i=sequences.begin (); i!=sequences.end (); ++i) {
      
      sequence_length.push_back (i->size ());

      for (j=i->begin (); j!=i->end (); ++j) {

    	CResId resid (k++);
 	translation[resid] = table[(CResId&)**j];  
	**j = resid;
	
	Conformation confo (&**j);

  	conformations.push_back (confo);
	graph.addNode (conformations.size () - 1);
	
	newcorresp[(const CResId&)(conformations.back ().getRes ())] = 
	  conformations.size () - 1;
  	sequence_mask.push_back (i-sequences.begin ());
      }
    }

    marks.resize ((int)conformations.size (), '-');
    helix_mask.resize ((int)conformations.size (), -1);
    strand_mask.resize ((int)conformations.size (), -1);
    tertiary_mask.resize ((int)conformations.size (), -1);

    for (k=0; k<(int)relations.size (); ++k) {
      if (relations[k].getDirection () == p_DIR_ANY) {
	graph.addEdge (newcorresp[(const CResId&)relations[k].getRef ()],
		       newcorresp[(const CResId&)relations[k].getRes ()],		      
		       k, false);	
      } else {
  	graph.addEdge (newcorresp[(const CResId&)relations[k].getRef ()],
  		       newcorresp[(const CResId&)relations[k].getRes ()],
  		       k, true);
	nb_connect++;
      }
      if (relations[k].isPairing ()) {
	marks[newcorresp[(const CResId&)relations[k].getRef ()]] = 'b';
	marks[newcorresp[(const CResId&)relations[k].getRes ()]] = 'b';
	nb_pairings++;
      }
    }
  }
}





void 
AnnotatedModel::findHelices ()
{
  Graph< node, edge >::adjgraph::iterator gi, gk, gip, gkp;
  Graph< node, edge >::adjlist::iterator gj, gjp;
  Helix helix;

  gi = graph.begin ();

  while (gi!=graph.end ()) {
    
//      cout << "( " << getResId (gi->first) << " " << flush;

    // Find a pairing involving gi
    for (gj=gi->second.begin (); gj!=gi->second.end (); ++gj) {
      if (isHelixPairing (gj->second)) break;
    }
      
    if (gj != gi->second.end ()) 
      {
	if (gj->first > gi->first &&
	    marks[gj->first] != '(' && marks[gj->first] != ')') {
	  // Find counterpart of gi.
	  gk = graph.find (gj->first);

//  	  cout << getResId (gk->first) << ")" << endl;

	  helix.push_back (make_pair (gi->first, gk->first));
	  
	  // Extending with bulge detection
	  while (true) {
	    gip = gi;  ++gip;
	    if (gip == graph.end ()) break;
	    if (gk == graph.begin ()) break;
	    gkp = gk;  --gkp;
	    
//  	    cout << "?  " << getResId(gip->first) << " " << getResId(gkp->first) << endl;

	    if (sequence_mask[gip->first] == sequence_mask[gi->first] &&
		sequence_mask[gkp->first] == sequence_mask[gk->first]) {
	      
	      if ((gjp = gip->second.find (gkp->first)) != gip->second.end () &&
		  isHelixPairing (gjp->second)) {
		helix.push_back (make_pair (gip->first, gkp->first));
		gi = gip;
		gk = gkp;

//  		cout << "+ " << getResId(gi->first) << " " << getResId(gk->first) << endl;

	      } else {

		//		cout << "Trying bulge in second part" << endl;
		// There may be a bulge in one of the two chains
 		Graph< node, edge >::adjlist::iterator ga, gb;
		bool found = false;
		for (ga=gip->second.begin (); ga!=gip->second.end (); ++ga) {

//  		  cout << getResId (ga->first) << " " << getResId (gk->first) << " " 
//  		       << "bulge length " << gk->first - ga->first << endl;

		  if (isHelixPairing (ga->second) &&
		      gk->first - ga->first <= 3 && 
		      (gb = gk->second.find (ga->first)) != gk->second.end () &&
		      relations[gb->second].isStacking ()) {
		    helix.push_back (make_pair (-1, ga->first-gk->first));
		    helix.push_back (make_pair (gip->first, gb->first));
		    gi = gip;
		    gk = graph.find (gb->first);

//  		    cout << "+ " << getResId(gi->first) << " " << getResId(gk->first) << endl;

		    found = true;
		    break;
		  }
		}
//  		cout << "Trying bulge in first part" << endl;
		if (!found) {
		  for (ga=gkp->second.begin (); ga!=gkp->second.end (); ++ga) {

//  		    cout << getResId (ga->first) << " " << getResId (gi->first) << " " 
//  			 << ga->first - gi->first << endl;

		    if (isHelixPairing (ga->second) &&
			ga->first - gi->first <= 3 && 
			(gb = gi->second.find (ga->first)) != gi->second.end () &&
			relations[gb->second].isStacking ()) {
		      helix.push_back (make_pair (gi->first-ga->first, -1));
		      helix.push_back (make_pair (gb->first, gkp->first));
		      gi = graph.find (gb->first);
		      gk = gkp;
//  		      cout << "+ " << getResId(gi->first) << " " << getResId(gk->first) << endl;
		      found = true;
		      break;
		    }
		  }
		}
		if (!found) {
		  break;
		}
	      }
	    } else {
	      break;
	    }
	  }
	  
	  if (helix.size () > 2) {
	    Helix::iterator i;
	    for (i=helix.begin (); i!=helix.end (); ++i) {
	      if (i->first >= 0) {
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



void AnnotatedModel::findStrands ()
{
  int i;
  Graph< node, edge >::adjgraph::iterator gi, gk;
  Strand strand;

  for (gi=graph.begin (); gi!=graph.end (); ) {
    if (helix_mask[gi->first] == -1) {
      gk = gi;
      while (gk != graph.end () &&
	     sequence_mask[gi->first] == sequence_mask[gk->first] &&
	     helix_mask[gk->first] == -1) ++gk;
      gk--;
      strand.first = gi->first;
      strand.second = gk->first;
      strand.type = OTHER;
      for (i=gi->first; i<=gk->first; ++i) {
	strand_mask[i] = strands.size ();
      }
      strands.push_back (strand);
      gi = gk;
    }
    gi++;
  }
}


void AnnotatedModel::classifyStrands ()
{
  vector< Strand >::iterator j;

  for (j=strands.begin (); j!=strands.end (); ++j) {
    if (j->first == 0 || j->second == (int)graph.size () - 1 ||
	sequence_mask[j->first-1] != sequence_mask[j->second+1])
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
	  cerr << "Pairing not found for " << j->first-1 << endl;

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
	  cerr << "Pairing not found for " << j->first-1 << endl;

	if (helices[helix_mask[j->second+1]].front ().first == j->second+1)
	  k = helices[helix_mask[j->second+1]].front ().second;
	else if (helices[helix_mask[j->second+1]].front ().second == j->second+1)
	  k = helices[helix_mask[j->second+1]].front ().first;
	else if (helices[helix_mask[j->second+1]].back ().first == j->second+1)
	  k = helices[helix_mask[j->second+1]].back ().second;
	else if (helices[helix_mask[j->second+1]].back ().second == j->second+1)
	  k = helices[helix_mask[j->second+1]].back ().first;
	else 
	  cerr << "Pairing not found for " << j->first-1 << endl;

	if (sequence_mask[i] != sequence_mask[k]) {
	  j->type = OTHER;
	} else if (i==k+1 || i==k-1) {
	  j->type = BULGE;
	} else {
	  // Find the other part of the internal loop
	  vector< Strand >::iterator s;
	  for (s=strands.begin (); s!=strands.end (); ++s) {
	    if (s->first == k+1 && s->second == i-1) {
	      j->type = INTERNAL_LOOP;
	      j->ref = s-strands.begin ();
	      break;
	    } else {
	      j->type = OTHER;
	    }
	  }
	}
      }    
  }
}



void AnnotatedModel::findKissingHairpins ()
{
  Graph< node, edge >::adjgraph::iterator gi;
  Graph< node, edge >::adjlist::iterator gj;
  int check = 0;

  cout << "Tertiary base-pairs ---------------------------------------------" << endl;

  for (gi=graph.begin (); gi!=graph.end (); ++gi) {
    for (gj=gi->second.begin (); gj!=gi->second.end (); ++gj) {

      if (gi->first < gj->first && isPairing (gj->second)) {
	
	if (helix_mask [gi->first] != -1 &&
	    helix_mask [gi->first] == helix_mask [gj->first]) {
	  // Simple helical base-pair.
	  check++;
	} 
	
	if (helix_mask [gi->first] != -1 &&
	    helix_mask [gi->first] != helix_mask [gj->first]) {
	  cout << setw (20) << "inter helix: " 
	       << getResId (gi->first) << "-" 
	       << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	  check++;
	}

	if (strand_mask [gi->first] != -1 &&
	    helix_mask [gj->first] != -1) {
	  if (strands[strand_mask [gi->first]].type == OTHER) {
	    cout << setw (20) << "strand/helix: " 
		 <<  getResId (gi->first) << "-" 
		 << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	    check++;
	  } else if (strands[strand_mask [gi->first]].type == BULGE ||
		     strands[strand_mask [gi->first]].type == BULGE_OUT) {
	    cout << setw (20) << "bulge/helix: " 
		 <<  getResId (gi->first) << "-" 
		 << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	    check++;
	  } else if (strands[strand_mask [gi->first]].type == LOOP) {
	    cout << setw (20) << "loop/helix: " 
		 <<  getResId (gi->first) << "-" 
		 << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	    check++;
	  } else if (strands[strand_mask [gi->first]].type == INTERNAL_LOOP) {
	    cout << setw (20) << "internal loop/helix: " 
		 <<  getResId (gi->first) << "-" 
		 << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	    check++;
	  } else {
	    cout << "Other = A" << endl;
	  }
	}
	
	if (strand_mask [gi->first] != -1 &&
	    strand_mask [gj->first] != -1) {
	  
	  if (strand_mask [gi->first] == strand_mask [gj->first]) {
	    cout << setw (20) << "intraloop: " 
		 << getResId (gi->first) << "-"
		 << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	    check++;
	  } else {
	    int size = 0;
	    if (strands[strand_mask [gi->first]].type == OTHER) {
	      cout << "strand/";
	      size = 20-7;	    
	    } else if (strands[strand_mask [gi->first]].type == BULGE || 
		     strands[strand_mask [gi->first]].type == BULGE_OUT) {
	      cout << "bulge/";
	      size = 20-6;
	    } else if (strands[strand_mask [gi->first]].type == LOOP) {
	      cout << "loop/";
	      size = 20-5;
	    } else if (strands[strand_mask [gi->first]].type == INTERNAL_LOOP) {
	      cout << "internal loop/";
	      size = 20-14;
	    }
	    if (strands[strand_mask [gj->first]].type == OTHER) {
	      cout << setw (size) << "strand: " 
		   <<  getResId (gi->first) << "-" 
		   << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	      check++;
	    } else if (strands[strand_mask [gj->first]].type == BULGE ||
		       strands[strand_mask [gj->first]].type == BULGE_OUT) {
	      cout << setw (size) << "bulge: " 
		   <<  getResId (gi->first) << "-" 
		   << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	      check++;
	    } else if (strands[strand_mask [gj->first]].type == LOOP) {
	      cout << setw (size) << "loop: " 
		   <<  getResId (gi->first) << "-" 
		   << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	      check++;
	    } else if (strands[strand_mask [gj->first]].type == INTERNAL_LOOP) {
	      cout << setw (size) << "internal loop: " 
		   <<  getResId (gi->first) << "-" 
		   << getResId (gj->first) << endl;
	  tertiary_mask[gi->first] = 1;
	  tertiary_mask[gj->first] = 1;
	      check++;
	    } else {
	      cout << "Other = B" << endl;
	    }
	  }
	}
      }
    }
  }

  if (nb_pairings != check) 
    cout << "Missing interactions: " << check << "/" << nb_pairings << endl;

}


void AnnotatedModel::findPseudoknots ()
{
  vector< Helix >::iterator i, j;
  
  for (i=helices.begin (); i!=helices.end (); ++i) {
    if (strand_mask[i->front ().first] == 
	strand_mask[i->front ().second]) { 
      node a, b, c, d;
      a = i->front ().first;
      b = i->back ().first;
      c = i->back ().second;
      d = i->front ().second;
      for (j=i; j!=helices.end (); ++j) {
	if (i!=j && strand_mask[j->front ().first] == 
	    strand_mask[j->front ().second]) {
	  node ap, bp, cp, dp;
	  ap = j->front ().first;
	  bp = j->back ().first;
	  cp = j->back ().second;
	  dp = j->front ().second;

	  if (ap > b && bp < c && cp > d) {
	    cout << "Pseudoknot : " << "H" << helix_mask[a]+1 << " " 
		 << "H" << helix_mask[ap]+1 << endl;
	  }
	}
      }
    }
  }
}



// I/O  -----------------------------------------------------------------




void AnnotatedModel::dumpSequences () 
{
  
  int i, j;
  int currseq = -1;

  cout << "Sequences -------------------------------------------------------" << endl;

  for (i=0; i<(int)conformations.size (); ) {
    
    currseq = sequence_mask[i];

    cout << "Sequence " << currseq+1 << " (length = " 
	 << sequence_length[currseq] << "): " << endl;
    
    int pos = 0;
    while (pos < sequence_length[currseq]) 
      {
	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
	  if ((j)%50==0) cout << setw (6) << getResId (i+j) << " ";
	  cout << *getType (i+j);
	  if ((j+1)%10==0) cout << " ";
	}
	cout << endl;
	
	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
	  if ((j)%50==0) cout << setw (6) << " " << " " << flush;
	  cout << marks[i+j] << flush;
	  if ((j+1)%10==0) cout << " " << flush;
	}
	cout << endl;

//  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
//  	  if ((j)%50==0) cout << setw (6) << "seq" << " " << flush;
//  	  if (sequence_mask[i+j] == -1) cout << "-" << flush;
//  	  else cout << sequence_mask[i+j]+1 << flush;
//  	  if ((j+1)%10==0) cout << " " << flush;
//  	}
//  	cout << endl;

//  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
//  	  if ((j)%50==0) cout << setw (6) << "hel" << " " << flush;
//  	  if (helix_mask[i+j] == -1) cout << "-" << flush;
//  	  else cout << helix_mask[i+j]+1 << flush;
//  	  if ((j+1)%10==0) cout << " " << flush;
//  	}
//  	cout << endl;

//  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
//  	  if ((j)%50==0) cout << setw (6) << "str" << " " << flush;
//  	  if (strand_mask[i+j] == -1) cout << "-" << flush;
//  	  else cout << strand_mask[i+j]+1 << flush;
//  	  if ((j+1)%10==0) cout << " " << flush;
//  	}
//  	cout << endl;

	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
	  if ((j)%50==0) cout << setw (6) << "ter" << " " << flush;
	  if (tertiary_mask[i+j] == -1) cout << "-" << flush;
	  else cout << 'X' << flush;
	  if ((j+1)%10==0) cout << " " << flush;
	}

	cout << endl << endl;

	pos = j+1;
      }
    i += pos-1;
  }
}


void AnnotatedModel::dumpGraph () 
{
  Graph< node, edge >::adjgraph::iterator i;
  Graph< node, edge >::adjlist::iterator j;
  PropertySet::iterator k;

  cout << "Graph of relations ----------------------------------------------" << endl;

  for (i=graph.begin (); i!=graph.end (); ++i) {
    for (j=i->second.begin (); j!=i->second.end (); ++j) {
      cout << getResId (i->first) << "-" 
	   << getResId (j->first) << " : ";
      for (k=relations[j->second].getGlobalProp ().begin (); 
	   k!=relations[j->second].getGlobalProp ().end (); ++k) {
	cout << (const char*)**k << " " ;
      }
      cout << endl;
    }
  } 
}


void AnnotatedModel::dumpConformations () 
{
  Graph< node, edge >::adjgraph::iterator i;
  cout << "Residue conformations -------------------------------------------" << endl;
  for (i=graph.begin (); i!=graph.end (); ++i) {
    cout << getResId (i->first) << " : " 
	 << *getType (i->first) << " " 
	 << *conformations[i->first].getPucker () << " " 
	 << *conformations[i->first].getGlycosyl () << endl; 
  }
}


void AnnotatedModel::dumpPairs () 
{
  Graph< node, edge >::adjgraph::iterator i;
  Graph< node, edge >::adjlist::iterator j;
  PropertySet::iterator k;

  cout << "Base-pairs ------------------------------------------------------" << endl;

  for (i=graph.begin (); i!=graph.end (); ++i) {
    for (j=i->second.begin (); j!=i->second.end (); ++j) {
      if (i->first < j->first && 
	  isPairing (j->second)) {
	cout << getResId (i->first) << "-" 
	     << getResId (j->first) << " : ";
	cout << (const char*)(*getType (i->first)) << "-" << flush
	     << (const char*)(*getType (j->first)) << " " << flush;
	  
	cout << (const char*)(*relations[j->second].getRefProp ()) << "/" << flush
	     << (const char*)(*relations[j->second].getResProp ()) << " " << flush;
	
	for (k=relations[j->second].getGlobalProp ().begin (); 
	     k!=relations[j->second].getGlobalProp ().end (); ++k) {
	  cout << (const char*)**k << " " ;
	}
	cout << endl;
      }
    } 
  }
}


void AnnotatedModel::dumpStacks () 
{
  Graph< node, edge >::adjgraph::iterator i, k;
  Graph< node, edge >::adjlist::iterator j;
  PropertySet::iterator l;

  cout << "Adjacent relations ----------------------------------------------" << endl;
  for (i=graph.begin (), k=i, ++k; i!=graph.end () && k!=graph.end (); ++i, ++k) {
    if ((j = i->second.find (k->first)) != i->second.end () &&
	relations[j->second].isAdjacent ()) {
      cout << getResId (i->first) << "-" 
	   << getResId (j->first) << " : ";
      for (l=relations[j->second].getGlobalProp ().begin (); 
	   l!=relations[j->second].getGlobalProp ().end (); ++l) {
	if (*l != p_DIR_5p)
	  cout << (const char*)**l << " " ;
      }
      cout << endl;
    }
  }
  cout << endl;

  cout << "Non-Adjacent stackings ------------------------------------------" << endl;
  for (i=graph.begin (); i!=graph.end (); ++i) {
    for (j=i->second.begin (); j!=i->second.end (); ++j) {
      if (i->first < j->first && 
	  relations[j->second].isStacking () && !relations[j->second].isAdjacent ()) {
	cout << getResId (i->first) << "-" 
	     << getResId (j->first) << " : ";
	for (l=relations[j->second].getGlobalProp ().begin (); 
	     l!=relations[j->second].getGlobalProp ().end (); ++l) {
	  if (*l != p_reverse)
	    cout << (const char*)**l << " " ;
	}
	cout << endl;
      }
    }
  }
}


void AnnotatedModel::dumpHelices () 
{
  vector< Helix >::iterator i;
  Helix::iterator j;
    
  cout << "Helices ---------------------------------------------------------" << endl;

  for (i=helices.begin (); i!=helices.end (); ++i) {
    char type;
    int tmpsize = 0;
    if (conformations[i->begin ()->first].getPucker ()->is_C3p_endo () &&
	conformations[i->begin ()->first].getGlycosyl ()->is_anti ())
      type = 'A';
    else if (conformations[i->begin ()->first].getPucker ()->is_C2p_endo () &&
	     conformations[i->begin ()->first].getGlycosyl ()->is_anti ())
      type = 'B';
    else type = 'X';

    for (j=i->begin (); j!=i->end (); ++j) {
      if (j->first >= 0) {
	char tmptype = 'X';
	if (conformations[j->first].getPucker ()->is_C3p_endo () &&
	    conformations[j->first].getGlycosyl ()->is_anti ())
	  tmptype = 'A';
	else if (conformations[j->first].getPucker ()->is_C2p_endo () &&
		 conformations[j->first].getGlycosyl ()->is_anti ())
	  tmptype = 'B';
	else tmptype = 'X';
	
	if (type != tmptype) type = 'X';
	
	if (conformations[j->second].getPucker ()->is_C3p_endo () &&
	    conformations[j->second].getGlycosyl ()->is_anti ())
	  tmptype = 'A';
	else if (conformations[j->second].getPucker ()->is_C2p_endo () &&
		 conformations[j->second].getGlycosyl ()->is_anti ())
	  tmptype = 'B';
	else tmptype = 'X';
	
	if (type != tmptype) type = 'X';
	
	tmpsize++;
      }
    }

    

    cout << "H" << i-helices.begin ()+1 << ", length = " 
	 << tmpsize << ", " << flush;
    
    cout << "type = " << type << " : " << endl;
    
    int k;


    cout.setf (ios::right, ios::adjustfield);
    cout << setw (6) << getResId ((*i)[0].first) << "-" << flush;
    for (k=0; k<(int)i->size (); ++k) {
      if ((*i)[k].first >=0)
	cout << *getType ((*i)[k].first) << flush;
      else {
	if (- (*i)[k].first -1 == 0) {
	  cout << "-" << setw (max ((int)floor(log10 (-(*i)[k].first)), (int)floor(log10 (-(*i)[k].second)))+1) 
	       << setfill ('-') << "-" << "-" << flush;
	} else {
	  cout << "-" << setw (max ((int)floor(log10 (-(*i)[k].first)), (int)floor(log10 (-(*i)[k].second)))+1) 
	       << - (*i)[k].first -1 << "-" << flush;
	}
      }
    }
    cout << "-" << getResId ((*i)[i->size ()-1].first) << endl;
    
    cout.setf (ios::right, ios::adjustfield);
    cout << setw (6) << getResId ((*i)[0].second) << "-" << flush;
    for (k=0; k<(int)i->size (); ++k) {
      if ((*i)[k].first >=0)
	cout << *getType ((*i)[k].second) << flush;
      else {
	if (- (*i)[k].second -1 == 0) {
	  cout << "-" << setw (max ((int)floor(log10 (-(*i)[k].first)), (int)floor(log10 (-(*i)[k].second)))+1) 
	       << setfill ('-') << "-" << "-" << flush;	
	} else {  
	  cout << "-" << setw (max ((int)floor(log10 (-(*i)[k].first)), (int)floor(log10 (-(*i)[k].second)))+1) 
	       << - (*i)[k].second -1 << "-" << flush;
	}
      }
    }
    cout << "-" << getResId ((*i)[i->size ()-1].second) << setfill (' ') << endl;
  }
}


void AnnotatedModel::dumpTriples () 
{
  Graph< node, edge >::adjgraph::iterator gi;
  Graph< node, edge >::adjlist::iterator gj;
  int count;
  int nb = 1;
  set< node > nodes;
  set< node >::iterator k, l;
  vector< bool > treated;

  treated.resize (graph.size (), false);

  cout << "Base Triples ----------------------------------------------------" << endl;
    
  for (gi=graph.begin (); gi!=graph.end (); ++gi) {
    count = 0;
    nodes.clear ();
    if (!treated[gi->first]) {
      nodes.insert (gi->first);
      for (gj=gi->second.begin (); gj!=gi->second.end (); ++gj) {
	if (isPairing (gj->second)) {
	  count++;
	  nodes.insert (gj->first);
	}
      }
      if (count > 1) {
	cout.setf (ios::left, ios::adjustfield);
	cout << "T" << setw (5) << nb++ << " "; 
	for (k=nodes.begin (); k!=nodes.end (); ++k) {	
	  for (l=k; l!=nodes.end (); ++l) {
	    if (k!=l && isPairing (*k, *l)) {
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
}

void 
AnnotatedModel::dumpStrands () 
{
  int j;
  int k;

  cout << "Strands ---------------------------------------------------------" << endl;

  for (j=0; j<(int)strands.size (); ++j) {

    cout.setf (ios::left, ios::adjustfield);
    cout << "S" << setw (5) << j+1 << " " <<  flush;

    cout.setf (ios::left, ios::adjustfield);

    if (strands[j].type == OTHER) cout << setw (16) << "single strand:";
    else if (strands[j].type == BULGE) cout << setw (16) << "bulge:";
    else if (strands[j].type == BULGE_OUT) cout << setw (16) << "bulge out:";
    else if (strands[j].type == LOOP) cout << setw (16) << "loop:";
    else if (strands[j].type == INTERNAL_LOOP) cout << setw (16) << "internal loop:";

    cout << getResId (strands[j].first) << "-" << flush;
    for (k=strands[j].first; k<=strands[j].second; ++k) {
      cout << *getType (k) << flush;
    }
    cout << "-" << getResId (strands[j].second);

    if (strands[j].type == INTERNAL_LOOP) {
      cout << " -- " << flush;
      cout << getResId (strands[strands[j].ref].first) << "-" << flush;
      for (k=strands[strands[j].ref].first; k<=strands[strands[j].ref].second; ++k) {
	cout << *getType (k) << flush;
      }
      cout << "-" << getResId (strands[strands[j].ref].second);
    }
    cout << endl;
  }
}



void 
AnnotatedModel::dumpMcc (const char* pdbname)
{
  char str[256];
  cout << "//" << endl;
  cout << "// Annotation results ------------------------------------------------" << endl;
  cout << "//" << endl;
  gethostname (str, 255);
  cout << "// Author         : " << getenv ("LOGNAME") << '@' << str << endl;
  cout << "// Creation date  : " << get_current_time (str) << endl;
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
    
    for (i=0; i<(int)conformations.size (); ) {
      
      currseq = sequence_mask[i];
      
      cout << "sequence( r " << endl;
      
      int pos = 0;
      while (pos < sequence_length[currseq]) 
	{
	  for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
	    if (j==0) cout << setw (8) << getResId (i+j) << " ";
	    cout << *getType (i+j);
	    if ((j+1)%10==0) cout << " ";
	  }
	  cout << endl;
	  pos = j+1;
	}
      
      cout << ")" << endl;
      
      i += pos-1;
    }
  }

  cout << "//" << endl;
  cout << "// Nucleotide conformations ------------------------------------------" << endl;
  cout << "// The distance between atoms O3' and P of two residues must be       " << endl;
  cout << "// inferior to 20nm for them to be considered adjacent.               " << endl;
  cout << "//" << endl;
  {
    Graph< node, edge >::adjgraph::iterator i;

    cout << "residue(" << endl;
    for (i=graph.begin (); i!=graph.end (); ++i) {
      cout << setw(8) << getResId (i->first) << " { ";
      cout.setf (ios::left, ios::adjustfield);
      cout << setw(8) << *conformations[i->first].getPucker () << " && " 
	   << setw(4) << *conformations[i->first].getGlycosyl () << " }" << "  100%" << endl; 
      cout.setf (ios::right, ios::adjustfield);
    }
    cout << ")" << endl;
  }

  cout << "//" << endl;
  cout << "// Non-adjacent base-pairs and stackings ------------------------" << endl;
  cout << "//" << endl;
 
  {
    Graph< node, edge >::adjgraph::iterator i;
    Graph< node, edge >::adjlist::iterator j;
    PropertySet::iterator k;

    if (nb_pairings > 0) {
      cout << "pair(" << endl;
      for (i=graph.begin (); i!=graph.end (); ++i) {
	for (j=i->second.begin (); j!=i->second.end (); ++j) {
	  if (i->first < j->first && 
	      !isAdjacent (j->second)) {
	    
	    cout << setw(8) << getResId (i->first) 
		 << setw(8) << getResId (j->first) << " { ";
	    
	    bool useand = false;
	    if (isPairing (j->second)) {
	      cout << (const char*)(*relations[j->second].getRefProp ()) << "/" << flush
		   << (const char*)(*relations[j->second].getResProp ()) << " " << flush;
	      useand = true;
	    }
	    for (k=relations[j->second].getGlobalProp ().begin (); 
		 k!=relations[j->second].getGlobalProp ().end (); ++k) {
	      if ((**k).is_reverse () || 
		  (**k).is_cis () ||
		  (**k).is_trans () ||
		  (**k).is_stack ()) {      
		if (useand) cout << "&& ";
		cout << (const char*)**k << " " ;
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
    Graph< node, edge >::adjgraph::iterator i;
    Graph< node, edge >::adjlist::iterator j;
    PropertySet::iterator k;
  
    if (nb_connect > 0) {
      cout << "connect(" << endl;
      for (i=graph.begin (); i!=graph.end (); ++i) {
	for (j=i->second.begin (); j!=i->second.end (); ++j) {
	  if (i->first < j->first && 
	      isAdjacent (j->second)) {
	    
	    cout << setw(8) << getResId (i->first) 
		 << setw(8) << getResId (j->first) << " { ";
	    
	    if (isPairing (j->second)) {
	      cout << (const char*)(*relations[j->second].getRefProp ()) << "/" << flush
		   << (const char*)(*relations[j->second].getResProp ()) << " " << flush;
	      cout << "&& ";
	    }
	    
	    if (isStacking (j->second))
	      cout << "stack" << " ";
	    else 
	      cout << "!stack" << " ";
	    
	    for (k=relations[j->second].getGlobalProp ().begin (); 
		 k!=relations[j->second].getGlobalProp ().end (); ++k) {
	      if ((**k).is_reverse () || 
		  (**k).is_cis () ||
		  (**k).is_trans () ||
		  (**k).is_pairing ())
		cout << "&& " << (const char*)**k << " " ;
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

    Graph< node, edge >::adjgraph::iterator i;
    Graph< node, edge >::adjlist::iterator j;
    
    for (i=graph.begin (); i!=graph.end (); ++i) {
      for (j=i->second.begin (); j!=i->second.end (); ++j) {
	if (isPairing (j->second))
	  graph.setEdgeValue (j->second, 1);
	else
	  graph.setEdgeValue (j->second, 2);
      }
    }

    vector< edge > edges;
    vector< edge >::reverse_iterator k;

    edges = graph.minimum_spanning_tree ();

    list< list< node > > treated;
    list< list< node > >::iterator x;
    list< node >::iterator y;

    for (k=edges.rbegin (); k!=edges.rend (); ++k) {
      node a, b;
      a = graph.getNodes (*k).first;
      b = graph.getNodes (*k).second;

      bool done = false;
      for (x=treated.begin (); x!=treated.end (); ++x) {
	if (x->front () == b) {
	  x->push_front (a);
	  done = true;
	} else if (x->back () == a) {
	  x->push_back (b);
	  done = true;
	}
	if (done) break;
      }
      if (!done) {
	list< node > tmp;
	tmp.push_back (a);
	tmp.push_back (b);
	treated.push_front (tmp);
      }
    }
    
    for (x=treated.begin (); x!=treated.end (); ++x) {
      cout << "  ( ";
      for (y=x->begin (); y!=x->end (); ++y) {
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

  {
    
    cout << "global_cache = cache(" << endl;
    cout << "  global" << endl;
    cout << "  rmsd (1.0 align base_only  no_hydrogen)" << endl;
    cout << ")" << endl;
  }
  
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
