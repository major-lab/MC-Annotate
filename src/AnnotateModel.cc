//    Copyright 2006 Laboratoire de Biologie Informatique et Theorique
// 
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
// 
//        http://www.apache.org/licenses/LICENSE-2.0
// 
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.


#include <config.h>

#include <algorithm>
#include <iterator>
#include <list>

#include "mccore/Binstream.h"
#include "mccore/Messagestream.h"
#include "mccore/Pdbstream.h"
#include "mccore/UndirectedGraph.h"
#include "mccore/stlio.h"

#include "AnnotateModel.h"


namespace annotate
{

  static const unsigned int MIN_HELIX_SIZE = 3;
  static const unsigned int PAIRING_MARK = 1;
  static const unsigned int LHELIX = 2;
  static const unsigned int RHELIX = 4;

  AbstractModel* AnnotateModelFM::createModel () const
  {
    return new AnnotateModel (residueSelection, environment, rFM);
  }

  AbstractModel* AnnotateModelFM::createModel (const AbstractModel &model) const
  {
    return new AnnotateModel (model, residueSelection, environment, rFM);
  }

  oBinstream& AnnotateModelFM::write (oBinstream& obs) const
  {
    return obs;
  }

  void AnnotateModel::annotate ()
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
  }

  void AnnotateModel::fillSeqBPStacks ()
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

  bool AnnotateModel::isHelixPairing (const Relation &r)
  {
    return(r.is (PropertyType::pPairing)
          && (r.is (PropertyType::pSaenger)
          || r.is (PropertyType::pOneHbond)));
  }

  void AnnotateModel::findHelices ()
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
  }

  void AnnotateModel::dumpHelices () const
  {
    /**
    * Output helices in text representation
    * Ex:
    * H0, length = 14
    *      A1-UUAUAUAUAUAUAA-A14
    *      B14-AAUAUAUAUAUAUU-B1
    */
    vector< Helix >::const_iterator i;

    for (i = helices.begin (); helices.end () != i; ++i)
    {
      Helix::const_iterator hIt;

      // Helix index and length
      gOut (0) << "H" << i - helices.begin () << ", length = " << i->size ()<< endl;

      // First strand
      gOut (0).setf (ios::right, ios::adjustfield);
      gOut (0) << setw (6) << internalGetVertex (i->front ().first)->getResId () << "-";

      for (hIt = i->begin (); i->end () != hIt; ++hIt)
      {
        const ResidueType *type = internalGetVertex (hIt->first)->getType ();
        gOut (0) << (type->isNucleicAcid () ? Pdbstream::stringifyResidueType (type) : "X");
      }
      --hIt;
      gOut (0) << "-" << internalGetVertex (hIt->first)->getResId () << endl;

      // Second strand
      gOut (0).setf (ios::right, ios::adjustfield);
      gOut (0) << setw (6) << internalGetVertex (i->front ().second)->getResId () << "-";

      for (hIt = i->begin (); i->end () != hIt; ++hIt)
      {
        const ResidueType *type = internalGetVertex (hIt->second)->getType ();
        gOut (0) << (type->isNucleicAcid () ? Pdbstream::stringifyResidueType (type) : "X");
      }
      --hIt;

      gOut (0) << "-" << internalGetVertex (hIt->second)->getResId() << endl;
    }
    gOut (0).setf (ios::left, ios::adjustfield);
  }

  void AnnotateModel::findStrands ()
  { }

  void AnnotateModel::buildStrands ()
  { }

  void AnnotateModel::findKissingHairpins ()
  { }

  void AnnotateModel::findPseudoknots ()
  { }

  void AnnotateModel::dumpSequences (bool detailed)
  { }

  void AnnotateModel::dumpConformations () const
  {
    const_iterator i;

    for (i = begin (); i != end (); ++i)
    {
      gOut(0) << i->getResId () << " : " << Pdbstream::stringifyResidueType (i->getType ());
      if (i->getType ()->isNucleicAcid ())
      {
        gOut (0) << " " << i->getPucker () << " " << i->getGlycosyl ();
      }
      gOut (0) << endl;
    }
  }

  void AnnotateModel::dumpNucleotides () const
  {
    const_iterator i;

    for (i = begin (); i != end (); ++i)
    {
       if (i->getType ()->isNucleicAcid ())
      {
        gOut(0) << i->getResId () << " " << Pdbstream::stringifyResidueType (i->getType ());
        gOut (0) << endl;
      }

    }
  }

  void AnnotateModel::dumpStacks () const
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

  void AnnotateModel::dumpPairs () const
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

  void AnnotateModel::dumpTriples ()
  { }

  ostream& AnnotateModel::output (ostream &os) const
  {
    gOut (0) << "Residue conformations -------------------------------------------" << endl;
    dumpConformations ();
    dumpStacks ();
    gOut (0) << "Base-pairs ------------------------------------------------------" << endl;
    dumpPairs ();
    return os;
  }

  ostream& AnnotateModel::outputNts(ostream &os) const
  {
    dumpNucleotides ();
    gOut (0) << ">Base pairs" << endl;
    dumpPairs ();
    gOut (0) << ">End";
    return os;
  }

  iPdbstream& AnnotateModel::input (iPdbstream &is)
  {
    return GraphModel::input (is);
  }

  iBinstream& AnnotateModel::input (iBinstream &is)
  {
    return GraphModel::input (is);
  }

}



namespace std
{
  ostream & operator<< (ostream &out, const annotate::Strand &t)
  {
    return t.output (out);
  }

  ostream& operator<< (ostream &os, const annotate::AnnotateModel &am)
  {
    return am.output (os);
  }
}