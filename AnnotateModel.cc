//                              -*- Mode: C++ -*- 
// AnnotatedModel.cc
// Copyright © 2001, 2002, 2003, 2004 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// Last Modified By : Martin Larose
// Last Modified On : Wed Dec  8 17:23:29 2004
// Update Count     : 208
// Status           : Unknown.
// 


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <list>
#include <utility>
#include <iterator>

#include "mccore/UndirectedGraph.h"
#include "mccore/GraphModel.h"
#include "mccore/Pdbstream.h"
#include "AnnotatedModel.h"


namespace annotate
{
  AnnotatedModel::AnnotatedModel (GraphModel &gfm) : gfm(gfm)
  {     
    nb_connect = 0;
    nb_pairings = 0;    
    GraphModel::edge_const_iterator edgeIt;
    for (edgeIt = gfm.edge_begin(); gfm.edge_end() != edgeIt; ++edgeIt, ++nb_connect) {
      if ((*edgeIt)->is(PropertyType::pPairing)) {
        marks[(*edgeIt)->getRef()] = 'b';
        marks[(*edgeIt)->getRes()] = 'b';
        nb_pairings++;
      }
    }
    nb_pairings /= 2;
    min_helix_size = 3;

    buildStrands();
    findHelices();
    findStrands();
    classifyStrands();
    gOut(0) << "Residue conformations -------------------------------------------" << endl;
    dumpConformations();
    dumpStacks();
    gOut(0) << "Base-pairs ------------------------------------------------------" << endl;
    findKissingHairpins();
    gOut(0) << "Triples ---------------------------------------------------------" << endl;
    dumpTriples();
    gOut(0) << "Helices ---------------------------------------------------------" << endl;
    dumpHelices();
    gOut(0) << "Strands ---------------------------------------------------------" << endl;
    dumpStrands();
    gOut(0) << "Various features ------------------------------------------------" << endl;
    findPseudoknots();
    gOut(0) << "Sequences -------------------------------------------------------" << endl;
    dumpSequences();
    gOut(0) << endl;
  }
  
  void 
  AnnotatedModel::findHelices ()
  {
    GraphModel::const_iterator gi, gk, gip, gkp;
    list< Residue * > neighbor;
    list< Residue * >::iterator nborIt;
    Helix helix;
  
    gi = gfm.begin ();
    while (gi!=gfm.end ()) {
      
      // Find a pairing involving gi
      neighbor = gfm.neighborhood ( (Residue *) &(*gi) );
      for (nborIt = neighbor.begin(); neighbor.end() != nborIt; ++nborIt) {
        if (isHelixPairing (gfm.getEdge ((Residue *)(&(*gi)), *nborIt)))
          break;
      }

      if (nborIt != neighbor.end ()) {
 
        if ((*gi).getResId() < (*nborIt)->getResId() 
            && marks.end() != marks.find(*nborIt) 
            && marks[*nborIt] != '(' && marks[*nborIt] != ')') {
          
         
         // Find counterpart of gi.
          gk = gfm.find ((*nborIt)->getResId());
      
          helix.push_back (make_pair (&(*gi), &(*gk)));
          
          // Extending with bulge detection
          while (true) {
             gip = gi;  ++gip;
            if (gip == gfm.end ()) break;
            if (gk == gfm.begin ()) break;
            gkp = gk;  --gkp;
        
            if (sequence_mask[&(*gip)] == sequence_mask[(&(*gi))] &&
                  sequence_mask[&(*gkp)] == sequence_mask[&(*gk)]) {
           
              neighbor = gfm.neighborhood ((Residue *) &(*gip));
              if ((nborIt = std::find (neighbor.begin (), neighbor.end (), (Residue *) &(*gkp))) != neighbor.end () &&
                  isHelixPairing (gfm.getEdge ((Residue *) &(*gip), (Residue *) &(*gkp)))) {
                  
                helix.push_back (make_pair ((Residue *) &(*gip), (Residue *) &(*gkp)));
                gi = gip;
                gk = gkp;
      
              } else {
                break;
              }
            } else {
              break;
            }
          }
      
          if ((int)helix.size () >= min_helix_size) {
            Helix::iterator i;
            for (i=helix.begin (); i!=helix.end (); ++i) {
 //             if (i->first >= 0) {
                helix_mask[i->first] = helices.size ();
                marks[i->first] = '(';
                helix_mask[i->second] = helices.size ();
                marks[i->second] = ')';
 //             }
            }
            helices.push_back (helix);
          } 
          helix.clear ();
        }
      }
      ++gi;
    }
  }

  /**
   * Output helices in text representation
   * Ex:
   * H0, length = 14
   *      A1-UUAUAUAUAUAUAA-A14
   *      B14-AAUAUAUAUAUAUU-B1
   */
  void
  AnnotatedModel::dumpHelices (void) const 
  {
    vector< Helix >::const_iterator i;
    
    for (i=helices.begin (); i!=helices.end (); ++i) {

      // Helix index and length
      gOut(0) << "H" << i - helices.begin () << ", length = " << i->size () << "" << endl;
    
      int k;
      // First strand
      gOut(0).setf (ios::right, ios::adjustfield);
      gOut(0) << setw (6) << ((*i)[0].first)->getResId() << "-";
      for (k=0; k<(int)i->size (); ++k) {
        if ((*i)[k].first >=0)
          gOut(0) << (((*i)[k].first)->getType()->isNucleicAcid () ? Pdbstream::stringifyResidueType (((*i)[k].first)->getType()) : "X");
        else 
        {
          // Code from old version of annotate.. ?
//            gOut(0) << "-" << setw (max ((int)floor(log10 ((double)-(*i)[k].first)), (int)floor(log10 ((double)-(*i)[k].second)))+1) 
//          if (- (*i)[k].first -1 == 0) {
//            gOut(0) << setfill ('-') << "-" << "-";
//          } else {
//           gOut(0)  << (*i)[k].first  << "-";
//          }
        }
      }
      gOut(0) << "-" << ((*i)[i->size ()-1].first)->getResId() << endl;

      // Second strand 
      gOut(0).setf (ios::right, ios::adjustfield);
      gOut(0) << setw (6) << ((*i)[0].second)->getResId() << "-";
    
      for (k=0; k<(int)i->size (); ++k) {
        if ((*i)[k].first >=0)
          gOut(0) << (((*i)[k].second)->getType()->isNucleicAcid () ? Pdbstream::stringifyResidueType (((*i)[k].second)->getType()) : "X");
        else {
          // Code from old version of annotate.. ?
//            gOut(0) << "-" << setw (max ((int)floor(log10 ((double)-(*i)[k].first)), (int)floor(log10 ((double)-(*i)[k].second)))+1) 
//          if (- (*i)[k].second -1 == 0) {
//              gOut(0) << setfill ('-') << "-" << "-";	
//          } else {  
//              gOut(0) << (*i)[k].second << "-";
//          }
        }
      }
      gOut(0) << "-" << ((*i)[i->size ()-1].second)->getResId() << setfill (' ') << endl;
    }
    gOut(0).setf (ios::left, ios::adjustfield);
  }
 
  void
  AnnotatedModel::findStrands (void)
  {

    GraphModel::const_iterator gi, gk, gl;
    map< const Residue *, int >::iterator gf;
    OStrand strand;

    for (gi=gfm.begin (); gi!=gfm.end (); )
    {
      if (helix_mask.find(&(*gi)) == helix_mask.end()) {
      
        gk = gi;
        while (gk != gfm.end () &&
            helix_mask.find(&(*gk)) == helix_mask.end() &&
            sequence_mask[&(*gi)] == sequence_mask[&(*gk)]
            )
          ++gk;
        gk--;
        strand.first = &(*gi);
        strand.second = &(*gk);
        strand.type = OTHER;
        for (gl = gi; gl<=gk; ++gl) { 
          strand_mask[&(*gl)] = strands.size ();
        }
        strands.push_back (strand);
        gi = gk;
      }
      gi++;
    }
  }
  
  void 
  AnnotatedModel::buildStrands(void)
  {
    list< GraphModel::edge_const_iterator > unsorted;
    GraphModel::edge_const_iterator edgeIt;

    sequences.clear ();
    for (edgeIt = gfm.edge_begin (); gfm.edge_end () != edgeIt; ++edgeIt)
    {
      if ((*edgeIt)->is (PropertyType::pAdjacent5p))
      {
        unsorted.push_back (edgeIt);
      }
    }
    while (! unsorted.empty ())
    {
      GraphModel::edge_const_iterator sj;
      list< GraphModel::edge_const_iterator > sorted;
      list< GraphModel::edge_const_iterator >::iterator it;
    
      sj = unsorted.front ();
      unsorted.pop_front ();
      sorted.push_back (sj);
      for (it = unsorted.begin (); unsorted.end () != it;)
      {
        if ((**it)->getRes ()->getResId () == (*sj)->getRef ()->getResId ())
        {
          sj = *it;
          sorted.push_front (sj);
          unsorted.erase (it);
          it = unsorted.begin ();
        }
        else
        {
          ++it;
        }
      }
      sj = sorted.back ();
      for (it = unsorted.begin (); unsorted.end () != it;)
      {
        if ((**it)->getRef ()->getResId () == (*sj)->getRes ()->getResId ())
        {
          sj = *it;;
          sorted.push_back (sj);
          unsorted.erase (it);
          it = unsorted.begin ();
        }
        else
        {
          ++it;
        }
      }
      sequences.push_back (Strand ());
      Strand &seq = sequences.back ();
    
      it = sorted.begin ();
      seq.push_back ((**it)->getRef ());
      while (sorted.end () != it)
      {
        seq.push_back ((**it)->getRes ());
        ++it;
      }
    }

    // Verify that the sequences have the same chain id and have consecutive
    // residue number.
    StrandSet::iterator seqIt;
    int k = 0;

    for (seqIt = sequences.begin (); sequences.end () != seqIt; ++seqIt, ++k)
    {
      Strand &resVec = *seqIt;
      Strand::iterator resVecIt;
      char chain;
      int number;
    int len = 1;
    
      resVecIt = resVec.begin ();
      chain = (*resVecIt)->getResId ().getChainId ();
      number = (*resVecIt)->getResId ().getResNo ();
      sequence_mask[*resVecIt] = k;
      ++resVecIt;
      while (resVec.end () != resVecIt)
      {
        if ((*resVecIt)->getResId ().getChainId () != chain
          || (*resVecIt)->getResId ().getResNo () != number + 1)
        {
          Strand rest;
      
          if (++resVec.begin () == resVecIt)
          {
            rest.insert (rest.end (), resVecIt, resVec.end ());
            resVec.erase (++resVecIt, resVec.end ());
          }
          else
          {
            rest.insert (rest.end (), resVecIt - 1, resVec.end ());
            resVec.erase (resVecIt, resVec.end ());
          }
          ++seqIt;
          if (1 < rest.size ())
          {
            seqIt = --sequences.insert (seqIt, rest);
          }
          break;
        }
        else
        {
          sequence_mask[*resVecIt] = k;
          number += 1;
          ++resVecIt;
        len++;
        }
      }
      sequence_length.push_back(len);
    }
  }

  void AnnotatedModel::classifyStrands ()  {
    vector< OStrand >::iterator j;
    GraphModel::const_iterator first, last, prev, next;
  
    for (j=strands.begin (); j!=strands.end (); ++j) {
  
      prev = gfm.find((*(j->first)).getResId());
      next = gfm.find((*(j->second)).getResId());
      --prev; ++next;
      first = gfm.begin();
      last = gfm.end();
      --last;
  
      if (&(*(j->first)) == &(*first) || &(*(j->second)) == &(*last) ||
        sequence_mask[&(*prev)] != sequence_mask[&(*next)])
      {
        j->type = OTHER;
      }
      else if (helix_mask[&(*prev)] == helix_mask[&(*next)])
      {
        const Residue *i = 0;
        // Find pairs of the extrimities.
        if (helices[helix_mask[&(*prev)]].front ().first == &(*prev))
          i = helices[helix_mask[&(*prev)]].front ().second;
        else if (helices[helix_mask[&(*prev)]].front ().second == &(*prev))
          i = helices[helix_mask[&(*prev)]].front ().first;
        else if (helices[helix_mask[&(*prev)]].back ().first == &(*prev))
          i = helices[helix_mask[&(*prev)]].back ().second;
        else if (helices[helix_mask[&(*prev)]].back ().second == &(*prev))
          i = helices[helix_mask[&(*prev)]].back ().first;
        else 
          cerr << "Pairing not found for (a) " << *prev << endl;
  
        if (i == &(*next))
          j->type = LOOP;
        else 
          j->type = BULGE_OUT;
      }
      else
      { 
        const Residue *i = 0, *k = 0;
    
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
      
        if (sequence_mask[i] != sequence_mask[k]) {
          j->type = OTHER;
        } else if (i==k+1 || i==k-1) {
          j->type = BULGE;
        } else {
          // Find the other part of the internal loop
          vector< OStrand >::iterator s;
          for (s=strands.begin (); s!=strands.end (); ++s) {
            if (s->first == k+1 && s->second == i-1) {
              j->type = INTERNAL_LOOP;
              j->ref = s-strands.begin ();  // FIXIT!
              break;
            } else {
              j->type = OTHER;
            }
          }
        }
      }    
    }
  }
  
  void AnnotatedModel::dumpStrands ()
  {
    int j;
    GraphModel::const_iterator k;
  
    for (j=0; j<(int)strands.size (); ++j) {
      
      gOut(0).setf (ios::left, ios::adjustfield);
      gOut(0) << "S" << setw (5) << j << " ";
      
      gOut(0).setf (ios::left, ios::adjustfield);
      
      if (strands[j].type == OTHER) gOut(0) << setw (16) << "single strand:";
      else if (strands[j].type == BULGE) gOut(0) << setw (16) << "bulge:";
      else if (strands[j].type == BULGE_OUT) gOut(0) << setw (16) << "bulge out:";
      else if (strands[j].type == LOOP) gOut(0) << setw (16) << "loop:";
      else if (strands[j].type == INTERNAL_LOOP) gOut(0) << setw (16) << "internal loop:";
      
      gOut(0) << (strands[j].first)->getResId() << "-";
      
      for (k = gfm.find((strands[j].first)->getResId()); k <= gfm.find((strands[j].second)->getResId()); ++k) {
        gOut(0) << ((k)->getType()->isNucleicAcid () ? Pdbstream::stringifyResidueType ((k)->getType()) : "X");
      }
      gOut(0) << "-" << (strands[j].second)->getResId();
  
      if (strands[j].type == INTERNAL_LOOP) {
        gOut(0) << " -- ";
      // FIXIT!    
  /*      gOut(0) << (strands[strands[j].ref].first)->getResId() << "-";
        for (k=strands[strands[j].ref].first; k<=strands[strands[j].ref].second; ++k) {
          gOut(0) << ((k)->getType()->isNucleicAcid () ? Pdbstream::stringifyResidueType ((k)->getType()) : "X");
        }
        gOut(0) << "-" << (strands[strands[j].ref].second)->getResId();
  */      
      }
      gOut(0) << endl;
    }
    gOut(0).setf (ios::left, ios::adjustfield);  
  }
  
  void AnnotatedModel::findKissingHairpins (void)
  {
    GraphModel::const_iterator gi;
    list< Residue * > neighbor;
    list< Residue * >::iterator nborIt;
    int check = 0;
    int nb_helical_bp = 0;
    set< const PropertyType* >::iterator k;
  
    gOut(0).setf (ios::left, ios::adjustfield);
  
    for (gi=gfm.begin (); gi!=gfm.end (); ++gi) {
      neighbor = gfm.neighborhood (const_cast<Residue *> (&(*gi)));  
      for (nborIt = neighbor.begin (); nborIt != neighbor.end (); ++nborIt) {      
        if (gi->getResId() < (*nborIt)->getResId() && isPairing (gfm.getEdge ((Residue *)(&(*gi)), *nborIt))) {
          
          if (helix_mask.find(&(*gi)) != helix_mask.end() &&
              helix_mask [&(*gi)] == helix_mask [*nborIt]) {
            // Simple helical base-pair.
            gOut(0) << setw (20) << "helix bp: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
            check++;
            nb_helical_bp++;
          } 
      
          if (helix_mask.find(&(*gi)) != helix_mask.end() &&
              helix_mask [&(*gi)] != helix_mask [*nborIt]) {
            gOut(0) << setw (20) << "inter helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
            tertiary_mask[&(*gi)] = 1;
            tertiary_mask[*nborIt] = 1;
            check++;
          }
  
          if (strand_mask.find(&(*gi)) != strand_mask.end() &&
              helix_mask.find(&(*gi)) != helix_mask.end()) {
            if (strands[strand_mask [&(*gi)]].type == OTHER) {
              gOut(0) << setw (20) << "strand/helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
            tertiary_mask[&(*gi)] = 1;
            tertiary_mask[*nborIt] = 1;
              check++;
            } else if (strands[strand_mask [&(*gi)]].type == BULGE ||
                strands[strand_mask [&(*gi)]].type == BULGE_OUT) {
              gOut(0) << setw (20) << "bulge/helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
            tertiary_mask[&(*gi)] = 1;
            tertiary_mask[*nborIt] = 1;
              check++;
            } else if (strands[strand_mask [&(*gi)]].type == LOOP) {
              gOut(0) << setw (20) << "loop/helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
            tertiary_mask[&(*gi)] = 1;
            tertiary_mask[*nborIt] = 1;
              check++;
            } else if (strands[strand_mask [&(*gi)]].type == INTERNAL_LOOP) {
              gOut(0) << setw (20) << "internal loop/helix: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
            tertiary_mask[&(*gi)] = 1;
            tertiary_mask[*nborIt] = 1;
              check++;
            } else {
              gOut(0) << "Other A: " << " : ";
            }
          }
  
          if (strand_mask.find(&(*gi)) != strand_mask.end() &&
              strand_mask.find(*nborIt) != strand_mask.end()) {
            
            if (strand_mask [&(*gi)] == strand_mask [*nborIt]) {
              gOut(0) << setw (20) << "intraloop: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
            tertiary_mask[&(*gi)] = 1;
            tertiary_mask[*nborIt] = 1;
              check++;
            } else {             
              int size = 0;
              if (strands[strand_mask [&(*gi)]].type == OTHER) {
                gOut(0) << "strand/";
                size = 20-7;	    
              } else if (strands[strand_mask [&(*gi)]].type == BULGE || 
                strands[strand_mask [&(*gi)]].type == BULGE_OUT) {
                gOut(0) << "bulge/";
                size = 20-6;
              } else if (strands[strand_mask [&(*gi)]].type == LOOP) {
                gOut(0) << "loop/";
                size = 20-5;
              } else if (strands[strand_mask [&(*gi)]].type == INTERNAL_LOOP) {
                gOut(0) << "internal loop/";
                size = 20-14;
              }              
              if (strands[strand_mask [*nborIt]].type == OTHER) {
                gOut(0) << setw (size) << "strand: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
                tertiary_mask[&(*gi)] = 1;
                tertiary_mask[*nborIt] = 1;
                check++;
               } else if (strands[strand_mask [*nborIt]].type == BULGE ||
                  strands[strand_mask [*nborIt]].type == BULGE_OUT) {
                gOut(0) << setw (size) << "bulge: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
                tertiary_mask[&(*gi)] = 1;
                tertiary_mask[*nborIt] = 1;
                check++;
              } else if (strands[strand_mask [*nborIt]].type == LOOP) {
                gOut(0) << setw (size) << "loop: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
                tertiary_mask[&(*gi)] = 1;
                tertiary_mask[*nborIt] = 1;
                check++;
              } else if (strands[strand_mask [*nborIt]].type == INTERNAL_LOOP) {
                gOut(0) << setw (size) << "internal loop: " << gi->getResId() << "-" << (*nborIt)->getResId() << " : ";
                tertiary_mask[&(*gi)] = 1;
                tertiary_mask[*nborIt] = 1;
                check++;
              } else {
                gOut(0) << "Other B:" << " : ";
              }
            }
          }
         
          gOut(0) << Pdbstream::stringifyResidueType (gi->getType()) << "-"
                  << Pdbstream::stringifyResidueType ((*nborIt)->getType()) << " ";
          
          Relation *e = gfm.getEdge ((Residue *) &(*gi), *nborIt);
                   
          if (isPairing(e))
            gOut(0) << e->getRefFace () << "/" << e->getResFace () << " ";
 
          const set< const PropertyType* > &labels = e->getLabels();
          copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (gOut(0), " "));          
          gOut(0) << endl;       
        }
      }
    }
  
    gOut(0) << "Number of base pairs = " << nb_pairings << endl;
    gOut(0) << "Number of helical base pairs = " << nb_helical_bp << endl;
  
    if (nb_pairings != check) 
      gOut(0) << "Missing interactions: " << check << "/" << nb_pairings << endl;
  }

  void AnnotatedModel::findPseudoknots (void)
  {
  vector< Helix >::iterator i, j;
  
  for (i=helices.begin (); i!=helices.end (); ++i) {
    if (strand_mask[i->front ().first] == 
	strand_mask[i->front ().second]) { 
      const Residue *a, *b, *c, *d;
      a = (*i).front ().first;
      b = (*i).back ().first;
      c = (*i).back ().second;
      d = (*i).front ().second;
      for (j=i; j!=helices.end (); ++j) {
	if (i!=j && strand_mask[j->front ().first] == 
	    strand_mask[j->front ().second]) {
	  const Residue *ap, *bp, *cp, *dp;
	  ap = (*j).front ().first;
	  bp = (*j).back ().first;
	  cp = (*j).back ().second;
	  dp = (*j).front ().second;
	  
	  if (b->getResId() < ap->getResId() && bp->getResId() < c->getResId() && d->getResId() < cp->getResId()) {
	    gOut(0) << "Pseudoknot : " << "H" << helix_mask[a]+1 << " " 
		 << "H" << helix_mask[ap]+1 << endl;
	  }
	}
      }
    }
  }
}
 
 void AnnotatedModel::dumpSequences (bool detailed) 
{
  GraphModel::iterator i;
  int j = 0;
  int currseq = -1;
  int size = 0;
  
  if (!detailed) {
//    for (i=0; i<(int)conformations.size (); ++i) {
//      if (i==0 || sequence_mask[i] != sequence_mask[i-1])
//	gOut(0) << " " << setw (8) << getResId (i) << " ";
//      gOut(0) << getType (i)->toString ();
//    }	
    return;
  } 
  
  for (i=gfm.begin (); i!=gfm.end (); ) {
    
    currseq = sequence_mask[(Residue *) &(*i)];
    size = sequence_length[currseq];

    gOut(0) << "Sequence " << currseq << " (length = " 
	 << size << "): " << endl;
    int pos = 0;
    while (pos < size) 
      {
	for (j=pos; (j<size && (j+1)%50!=0); ++j) {
	  if ((j)%50==0) {
	    gOut(0).setf (ios::right, ios::adjustfield);
	    gOut(0) << setw (3) << (*(i+j))->getResId().getChainId ();
	    gOut(0).setf (ios::left, ios::adjustfield);
	    gOut(0) << setw (5) << (*(i+j))->getResId().getResNo () << " ";
	  }
	  
	  gOut(0) << ((*(i+j))->getType()->isNucleicAcid () ? Pdbstream::stringifyResidueType ((*(i+j))->getType()) : "X");
	  
	  if ((j+1)%10==0) gOut(0) << " ";
	}
	gOut(0) << endl;
	
	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
	  if ((j)%50==0) gOut(0) << setw (8) << " " << " ";
 	  if (marks.find(*(i+j)) == marks.end()) gOut(0) << "-";
          else gOut(0) << marks[*(i+j)];
	  if ((j+1)%10==0) gOut(0) << " ";
	}
	gOut(0) << endl;

//  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
//  	  if ((j)%50==0)
//          gOut(0) << setw (6) << "seq" << " ";
//  	  if (sequence_mask[i+j] == -1)
//          gOut(0) << "-";
//  	  else gOut(0) << sequence_mask[i+j];
//  	  if ((j+1)%10==0)
//          gOut(0) << " ";
//  	}
//  	gOut(0) << endl;

 	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
 	  if ((j)%50==0) gOut(0) << setw (8) << "helix" << " ";
 	  if (helix_mask.find(*(i+j)) == helix_mask.end()) gOut(0) << "-";
 	  else gOut(0) << helix_mask[*(i+j)];
 	  if ((j+1)%10==0) gOut(0) << " ";
 	}
 	gOut(0) << endl;

//  	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
//  	  if ((j)%50==0) gOut(0) << setw (6) << "str" << " ";
//  	  if (strand_mask[i+j] == -1) gOut(0) << "-";
//  	  else gOut(0) << strand_mask[i+j];
//  	  if ((j+1)%10==0) gOut(0) << " ";
//  	}
//  	gOut(0) << endl;

// 	for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
// 	  if ((j)%50==0) gOut(0) << setw (6) << "ter" << " ";
// 	  if (tertiary_mask[i+j] == -1) gOut(0) << "-";
// 	  else gOut(0) << 'X';
// 	  if ((j+1)%10==0) gOut(0) << " ";
// 	}

// 	gOut(0) << endl << endl;

	pos = j+1;
      }
    i += pos-1;
  }
  gOut(0).setf (ios::left, ios::adjustfield);
}

void AnnotatedModel::dumpConformations (void) 
{
  GraphModel::const_iterator i;
  
  for (i=gfm.begin (); i!=gfm.end (); ++i) {
    gOut(0) << (*i).getResId () << " : " 
	 << (*i).getType() << " " 
	 << (*i).getPucker() << " " 
	 << (*i).getGlycosyl() << endl; 
  }
}

void AnnotatedModel::dumpPairs (void) 
{
  GraphModel::const_iterator i;
  GraphModel::edge_const_iterator edgeIt;
  list< Residue * > neighbor;
  list< Residue * >::iterator nborIt;
  set< const PropertyType* >::iterator k;

  for (i=gfm.begin (); i!=gfm.end (); ++i) {

    neighbor = gfm.neighborhood (( (Residue *) &(*i) ));

    for (nborIt=neighbor.begin (); nborIt!=neighbor.end (); ++nborIt) {
      if ((*i).getResId() < (*nborIt)->getResId() && isPairing ((Residue *) &(*i), *nborIt)) {
	gOut(0) << (*i).getResId() << "-" << (*nborIt)->getResId() << " : ";

	gOut(0) << Pdbstream::stringifyResidueType ((*i).getType()) << "-"
	     << Pdbstream::stringifyResidueType ((*nborIt)->getType()) << " ";
	  
        Relation *e = gfm.getEdge ((Residue *) &(*i), *nborIt);
        gOut(0) << e->getRefFace () << "/" << e->getResFace () << " ";

        const set< const PropertyType* > &labels = e->getLabels();
        copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (gOut(0), " "));          
	gOut(0) << endl;
      }
    } 
  }
}

void AnnotatedModel::dumpStacks (void) 
{
  GraphModel::const_iterator i, k;
  list< Residue * > neighbor;
  list< Residue * >::iterator nborIt;
  set< const PropertyType* >::iterator l;

  int nb_stacks = 0;
  int nb_helical_stacks = 0;
  int nb_adjacent_stacks = 0;

  gOut(0) << "Adjacent stackings ----------------------------------------------" << endl;

  for (i=gfm.begin (), k=i, ++k; i!=gfm.end () && k!=gfm.end (); ++i, ++k) {
    
    neighbor = gfm.neighborhood (( (Residue *) &(*i) ));

    if ((nborIt = std::find (neighbor.begin (), neighbor.end (), (Residue *) &(*k))) != neighbor.end () 
	&& gfm.getEdge ((Residue *) &(*i), *nborIt)->is (PropertyType::pAdjacent))
    {
      
      if (gfm.getEdge ((Residue *) &(*i), *nborIt)->is (PropertyType::pStack)) {
	nb_adjacent_stacks++;
	nb_stacks++;
	if (helix_mask.end() != helix_mask.find(&(*i))
            && helix_mask[&(*i)] == helix_mask[*nborIt])
	  nb_helical_stacks++;
      }
      
      gOut(0) << (*i).getResId () << "-" << (*nborIt)->getResId () << " : ";
      Relation *e = gfm.getEdge ((Residue *) &(*i), (Residue *) *nborIt);
      const set< const PropertyType* > &labels = e->getLabels();
      copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (gOut(0), " "));          
      gOut(0) << endl;
    }
  }
  gOut(0) << endl;

  gOut(0) << "Non-Adjacent stackings ------------------------------------------" << endl;
  for (i=gfm.begin (); i!=gfm.end (); ++i) {

    neighbor = gfm.neighborhood (( (Residue *) &(*i) ));
    
    for (nborIt=neighbor.begin (); nborIt!=neighbor.end (); ++nborIt) {
      if ((*i).getResId() < (*nborIt)->getResId() && 
	  gfm.getEdge ((Residue *) &(*i), *nborIt)->is (PropertyType::pStack) && 
	  !gfm.getEdge ((Residue *) &(*i), *nborIt)->is (PropertyType::pAdjacent)) {
	nb_stacks++;
	
      gOut(0) << (*i).getResId () << "-" << (*nborIt)->getResId () << " : ";
      Relation *e = gfm.getEdge ((Residue *) &(*i), (Residue *) *nborIt);
      const set< const PropertyType* > &labels = e->getLabels();
      copy (labels.begin (), labels.end (), ostream_iterator< const PropertyType* > (gOut(0), " "));          
      gOut(0) << endl;
      }
    }
  }

  gOut(0) << "Number of stackings = " << nb_stacks << endl
          << "Number of helical stackings = " << nb_helical_stacks << endl
          << "Number of adjacent stackings = " << nb_adjacent_stacks << endl
          << "Number of non adjacent stackings = " << nb_stacks - nb_adjacent_stacks << endl;
}

void AnnotatedModel::dumpTriples (void) 
{
/*
  GraphModel::const_iterator gi;
  list< Residue * > neighbor;
  list< Residue * >::iterator nborIt;

  int count;
  int nb = 1;
  set< node > nodes;
  set< node >::iterator k, l;
  vector< bool > treated;

  treated.resize (gfm.size (), false);

  for (gi=gfm.begin (); gi!=gfm.end (); ++gi) {
    count = 0;
    nodes.clear ();
    if (!treated[*gi]) {
      nodes.insert (*gi);
      
      neighbor = gfm.getNeighbors (*gi);

      for (nborIt=neighbor.begin (); nborIt!=neighbor.end (); ++nborIt) {
	if (isPairing (gfm.getEdge (*gi,*nborIt))) {
	  count++;
	  nodes.insert (*nborIt);
	}
      }
      if (count > 1) {
	gOut(0).setf (ios::left, ios::adjustfield);
	gOut(0) << "T" << setw (5) << nb++ << " "; 
	for (k=nodes.begin (); k!=nodes.end (); ++k) {	
	  for (l=k; l!=nodes.end (); ++l) {
	    if (k!=l && isPairing (*k, *l)) {
	      gOut(0) << getResId (*k) << "-"
		   << getResId (*l) << " ";
	      treated[*k] = true;
	      treated[*l] = true;
	    }
	  }
	}
	gOut(0) << endl;
      }
    }
  }
  gOut(0).setf (ios::left, ios::adjustfield);
*/
}

  ostream&
  AnnotatedModel::output (ostream &os) const
  {
    return os;
  }


  iPdbstream&
  AnnotatedModel::input (iPdbstream &ips)
  {
    return ips;
  }
  
  iBinstream&
  AnnotatedModel::input (iBinstream &is)
  {
    return is;
  }
  
  oBinstream&
  AnnotatedModel::output (oBinstream &os) const
  {
    return os;
  }
 
  iBinstream&
  operator>> (iBinstream &is, AnnotatedModel &model)
  {
    return model.input (is);
  }

  oBinstream&
  operator<< (oBinstream &os, const AnnotatedModel &model)
  {
    return model.output (os);
  }

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
  operator<< (ostream &os, const annotate::AnnotatedModel &am)
  {
    return am.output (os);
  }  
}
  
namespace annotate {

void 
AnnotatedModel::dumpMcc (const char* pdbname)
{
//   char str[256];
//   cout << "//" << endl;
//   cout << "// Annotation results ------------------------------------------------" << endl;
//   cout << "//" << endl;
//   gethostname (str, 255);
//   cout << "// Author         : " << getenv ("LOGNAME") << '@' << str << endl;
//   cout << "// Structure file : " << pdbname << endl;
//   cout << "// Structural annotation generated by " << PACKAGE << " "  << VERSION << endl;
//   cout << "// ";
//   cout << endl;
// 
//   cout << "// Sequences ---------------------------------------------------------" << endl;
//   cout << "// The distance between atoms O3' and P of two residues must be       " << endl;
//   cout << "// inferior to 20nm for them to be considered adjacent.               " << endl;
//   cout << "//" << endl;
//   {
//     int i, j;
//     int currseq = -1;
//     
//     for (i=0; i<(int)conformations.size (); ) {
//       
//       currseq = sequence_mask[i];
//       
//       cout << "sequence( r " << endl;
//       
//       int pos = 0;
//       while (pos < sequence_length[currseq]) 
// 	{
// 	  for (j=pos; (j<sequence_length[currseq] && (j+1)%50!=0); ++j) {
// 	    if (j==0) cout << setw (8) << getResId (i+j) << " ";
// 	    cout << getType (i+j)->toString ();
// 	    if ((j+1)%10==0) cout << " ";
// 	  }
// 	  cout << endl;
// 	  pos = j+1;
// 	}
//       
//       cout << ")" << endl;
//       
//       i += pos-1;
//     }
//   }
// 
//   cout << "//" << endl;
//   cout << "// Nucleotide conformations ------------------------------------------" << endl;
//   cout << "// The distance between atoms O3' and P of two residues must be       " << endl;
//   cout << "// inferior to 20nm for them to be considered adjacent.               " << endl;
//   cout << "//" << endl;
//   {
//     UndirectedGraph< node, edge >::iterator i;
// 
//     cout << "residue(" << endl;
//     for (i=graph.begin (); i!=graph.end (); ++i) {
//       cout << setw(8) << getResId (*i) << " { ";
//       cout.setf (ios::left, ios::adjustfield);
//       cout << setw(8) << *conformations[*i]->getPucker () << " && " 
// 	   << setw(4) << *conformations[*i]->getGlycosyl () << " }" << "  100%" << endl; 
//       cout.setf (ios::right, ios::adjustfield);
//     }
//     cout << ")" << endl;
//   }
// 
//   cout << "//" << endl;
//   cout << "// Non-adjacent base-pairs and stackings ------------------------" << endl;
//   cout << "//" << endl;
//  
//   {
//     UndirectedGraph< node, edge >::iterator i;
//     list< node >::iterator j;
//     list< node > neigh;
//     set< const PropertyType* >::iterator k;
//     
//     if (nb_pairings > 0) {
//       cout << "pair(" << endl;
//       for (i=graph.begin (); i!=graph.end (); ++i) {
// 	neigh = graph.getNeighbors (*i);
// 	for (j=neigh.begin (); j!=neigh.end (); ++j) {
// 	  edge e = graph.getEdge (*i, *j);
// 	  if (*i < *j && 
// 	      !isAdjacent (e)) {
// 	    
// 	    cout << setw(8) << getResId (*i) 
// 		 << setw(8) << getResId (*j) << " { ";
// 	    
// 	    bool useand = false;
// 	    if (isPairing (e) && !relations[e].is (PropertyType::parseType ("unclassified"))) {
// 	      cout << (const char*)(*relations[e].getRefFace ()) << "/" << flush
// 		   << (const char*)(*relations[e].getResFace ()) << " " << flush;
// 	      useand = true;
// 	    }
// 	    for (k=relations[e].getLabels ().begin (); 
// 		 k!=relations[e].getLabels ().end (); ++k) {
// 	      if ((*k) == PropertyType::pReverse || 
// 		  (*k) == PropertyType::pCis ||
// 		  (*k) == PropertyType::pTrans ||
// 		  (*k) == PropertyType::pStack) {      
// 		if (useand) cout << "&& ";
// 		cout << (const char*)**k << " " ;
// 		useand = true;
// 	      }
// 	    }
// 	    
// 	    cout << "}" << "  100%" << endl; 
// 	  }
// 	}
//       }
//       cout << ")" << endl;
//     }
//   }
//   
//   cout << "//" << endl;
//   cout << "// Adjacent relations -------------------------------------------" << endl;
//   cout << "//" << endl;
// 
//   {
//     UndirectedGraph< node, edge >::iterator i;
//     list< node >::iterator j;
//     list< node > neigh;
//     set< const PropertyType* >::iterator k;
//   
//     if (nb_connect > 0) {
//       cout << "connect(" << endl;
//       for (i=graph.begin (); i!=graph.end (); ++i) {
// 	neigh = graph.getNeighbors (*i);
// 	for (j=neigh.begin (); j!=neigh.end (); ++j) {
// 	  edge e = graph.getEdge (*i, *j);
// 	  if (*i < *j && 
// 	      isAdjacent (e)) {
// 	    
// 	    cout << setw(8) << getResId (*i) 
// 		 << setw(8) << getResId (*j) << " { ";
// 	    
// 	    if (isPairing (e)) {
// 	      cout << (const char*)(*relations[e].getRefFace ()) << "/" << flush
// 		   << (const char*)(*relations[e].getResFace ()) << " " << flush;
// 	      cout << "&& ";
// 	    }
// 	    
// 	    if (isStacking (e))
// 	      cout << "stack" << " ";
// 	    else 
// 	      cout << "!stack" << " ";
// 	    
// 	    for (k=relations[e].getLabels ().begin (); 
// 		 k!=relations[e].getLabels ().end (); ++k) {
// 	      if ((*k) == PropertyType::pReverse || 
// 		  (*k) == PropertyType::pCis ||
// 		  (*k) == PropertyType::pTrans ||
// 		  (*k) == PropertyType::pPairing) {
// 		cout << "&& " << (const char*)**k << " " ;
// 	      }
// 	    }
// 	    cout << "}" << "  100%" << endl; 
// 	  }
// 	} 
//       }
//       cout << ")" << endl;    
//     }
//   }
//   
//   cout << "//" << endl;
//   cout << "// Construction order -------------------------------------------" << endl;
//   cout << "// This section defines a spanning tree of minimal height that   " << endl;
//   cout << "// connects all residues of the molecule and states the order in " << endl;
//   cout << "// which they are placed in the modeling process."                 << endl;
//   cout << "//" << endl;
// 
//   // Here, we need to find connected components of the initial graph and build a backtrack
//   // for each of them.  The generated script is not guaranteed to work when an annotated
//   // PDB file contains many fragments of RNA...  What should we do???????
// 
//   {
//     cout << "global = backtrack(" << endl;
// 
//     UndirectedGraph< node, edge >::iterator i;
//     list< node >::iterator j;
//     list< node > neigh;
//     
//     for (i=graph.begin (); i!=graph.end (); ++i) {
//       neigh = graph.getNeighbors (*i);
//       for (j=neigh.begin (); j!=neigh.end (); ++j) {
// 	edge e = graph.getEdge (*i, *j);
// 	if (isPairing (e) && !relations[e].is (PropertyType::parseType ("unclassified")))
// 	  graph.setEdgeWeight (*i, *j, 1);
// 	else
// 	  graph.setEdgeWeight (*i, *j, 2);
//       }
//     }
// 
//     vector< pair< node, node > > edges;
//     vector< pair< node, node > >::reverse_iterator k;
// 
//     edges = graph.minimumSpanningTree ();
// 
//     list< list< node > > treated;
//     list< list< node > >::iterator x;
//     list< node >::iterator y;
// 
//     for (k=edges.rbegin (); k!=edges.rend (); ++k) {
//       node a, b;
//       a = *graph.find (k->first);
//       b = *graph.find (k->second);
// 
//       bool done = false;
//       for (x=treated.begin (); x!=treated.end (); ++x) {
// 	if (x->front () == b) {
// 	  x->push_front (a);
// 	  done = true;
// 	} else if (x->back () == a) {
// 	  x->push_back (b);
// 	  done = true;
// 	}
// 	if (done) break;
//       }
//       if (!done) {
// 	list< node > tmp;
// 	tmp.push_back (a);
// 	tmp.push_back (b);
// 	treated.push_front (tmp);
//       }
//     }
// 
//     // correction pass to keep a valid backtrack statement
// 
//     size_t placed_sz = 0;
//     set< int > placed;
// 
//     for (y = treated.begin ()->begin (); y != treated.begin ()->end (); ++y)
//       placed.insert (*y);
// 
//     if (treated.size () > 1)
//       {
// 	x = treated.begin ();
// 	x++;
// 	for (; x != treated.end (); ++x) 
// 	  {
// 	    placed_sz = placed.size ();
// 	    placed.insert (*x->begin ());
// 
// 	    if (placed.size () > placed_sz)
// 	      {
// 		// oups! reference residue not placed yet!
// 		// What if we just flip the sub-list over...
// 		placed_sz = placed.size ();
// 		placed.insert (x->back ());
// 
// 		if (placed.size () > placed_sz)
// 		  {
// 		    // hum...something is really wrong here!
// 		    cerr << "Fatal Error: unable to build a valid backtrack statement" << endl;
// 		    exit (EXIT_FAILURE);
// 		  }
// 
// 		// update placed residues set, then flip the sub-list over
// 		list< int > tmp;
// 		for (y = x->begin (); y != x->end (); ++y)
// 		  {
// 		    placed.insert (*y);
// 		    tmp.push_front (*y);
// 		  }
// 		x->clear ();
// 		for (y = tmp.begin (); y != tmp.end (); ++y)
// 		  x->push_back (*y);
// 	      }
// 	    else
// 	      {
// 		// just update placed residues set
// 		for (y = x->begin (); y != x->end (); ++y)
// 		  placed.insert (*y);
// 	      }
// 
// 	  }
//       }
//     
//     for (x=treated.begin (); x!=treated.end (); ++x) {
//       cout << "  ( ";
//       for (y=x->begin (); y!=x->end (); ++y) {
// 	cout << getResId (*y) << " ";
//       }
//       cout << ")" << endl;
//     }
//     
//     cout << ")" << endl;
//   }
//   
//   cout << "//" << endl;
//   cout << "// Molecule cache -----------------------------------------------" << endl;
//   cout << "// A cache is normally placed on top of the backtrack so         " << endl;
//   cout << "// generated models are kept only if they are different enough   " << endl;
//   cout << "// from previously generated models."                              << endl;
//   cout << "//" << endl;
// 
//   {
//     
//     cout << "global_cache = cache(" << endl;
//     cout << "  global" << endl;
//     cout << "  rmsd (1.0 align base_only  no_hydrogen)" << endl;
//     cout << ")" << endl;
//   }
//   
//   cout << "//" << endl;
//   cout << "// Constraints --------------------------------------------------" << endl;
//   cout << "//" << endl;
//   cout << "adjacency( " << endl
//        << "  global 1.0 5.0" << endl
//        << ")" << endl << endl;
//   cout << "res_clash(" << endl
//        << "  global" << endl
//        << "  fixed_distance 1.0" << endl
//        << "  all no_hydrogen" << endl
//        << ")" << endl;
//   cout << "//" << endl;
//   cout << "// Exploration type ---------------------------------------------" << endl;
//   cout << "//" << endl;
//   cout << "explore(" << endl
//        << "  global_cache" << endl
//        << "  file_pdb (\"global-%05d.pdb\" zipped)" << endl
//        << ")" << endl;
// 
//   cout << "//" << endl
//        << "// --------------------------------------------------------------" << endl
//        << "//" << endl;
}



// void AnnotatedModel::dumpCt (const char* pdbname)
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


// void AnnotatedModel::PDF_drawLoop2 (PDF *p, int li, int source_jct)
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



// void AnnotatedModel::PDF_drawHelix (PDF *p, int hi, int li_ref, int x, int y)
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



// void AnnotatedModel::PDF_drawLoop (PDF *p, int li, int hi_ref, int x, int y)
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



// void AnnotatedModel::dumpPDF (const char* pdfname)
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



// // void AnnotatedModel::PDF_drawGraph (PDF *p, node curr, node prev, int x, int y, int a)
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


// // void AnnotatedModel::dumpSimplePDF (const char* pdbname, const char* pdfname)
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

}
