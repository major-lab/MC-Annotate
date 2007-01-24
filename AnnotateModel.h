//                    -*- Mode: C++; coding: UTF-8 -*- 
// AnnotateModel.h
// Copyright © 2001-07 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// $Revision$
// $Id$


#ifndef _annotate_AnnotateModel_h_
#define _annotate_AnnotateModel_h_

#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "mccore/GraphModel.h"
#include "mccore/Model.h"
#include "mccore/ModelFactoryMethod.h"
#include "mccore/PdbFileHeader.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"
#include "mccore/ResId.h"
#include "mccore/ResIdSet.h"
#include "mccore/Residue.h"
#include "mccore/ResidueType.h"

#include "AnnotateResidue.h"
#include "BaseLink.h"
#include "BasePair.h"
#include "BaseStack.h"
#include "Helix.h"
#include "Sequence.h"

using namespace mccore;
using namespace std;



namespace mccore
{
  class iBinstream;
  class iPdbstream;
}



namespace annotate
{
  class AnnotateModel;
  
  typedef int strandId;
  
  enum stype { BULGE_OUT, BULGE, INTERNAL_LOOP, LOOP, HELIX, OTHER };
 
  class Strand : public vector< const Residue* >
  {
  private:
    char consensusChainId;
  
  public:
     
    virtual ostream& output (ostream &os) const
    {
      Strand::const_iterator resIt;
        
      os << size () << " residues : ";
      for (resIt = begin (); end () != resIt; ++resIt)
        {
          os << (**resIt).getType();
          if (end () != resIt + 1)
	    {
	      os << "-";
	    }
        }
      return os;
    }
      
    virtual ~Strand () { }
    
  };

  class StrandSet : public vector< Strand >
  {
    /**
     * Map from integer to original ResId.  The ResId of the residues are
     * replaced with sequential ids.
     */
    map< unsigned int, const ResId * > int2ResIdMap;
      
  };

  /**
   * @short ModelFactoryMethod implementation for AnnotateModel class.
   *
   * This is the model factory method implementation for the AnnotateModel
   * class.
   *
   * @author Martin Larose (<a href="larosem@iro.umontreal.ca">larosem@iro.umontreal.ca</a>)
   * @version $Id$
   */
  class AnnotateModelFM : public ModelFactoryMethod
  {

    /**
     * A ResIdSet of residues to annotate.
     */
    ResIdSet residueSelection;

    /**
     * The number of relation layers around the residue selection to
     * annotate.
     */
    unsigned int environment;

  public:

    // LIFECYCLE ------------------------------------------------------------

    /**
     * Initializes the object.
     * @param rs a ResIdSet of residue ids to annotate.
     * @param env the number of relation layers around the residue selection to
     * annotate.
     * @param fm the residue factory method.
     */
    AnnotateModelFM (const ResIdSet &rs, unsigned int env)
      : ModelFactoryMethod (new AnnotateResidueFM ()),
	residueSelection (rs),
	environment (env)
    { }

    /**
     * Initializes the object with the right content.
     * @param right the object to copy.
     */
    AnnotateModelFM (const ModelFM &right) : ModelFactoryMethod (right) { }

    /**
     * Clones the object.
     * @return the copy of the object.
     */
    virtual ModelFactoryMethod* clone () const
    {
      return (ModelFactoryMethod*) new AnnotateModelFM (*this);
    }
  
    /**
     * Destroys the object.
     */
    virtual ~AnnotateModelFM () { }

    // OPERATORS ------------------------------------------------------------

    // ACCESS ---------------------------------------------------------------

    // METHODS --------------------------------------------------------------

    /**
     * Creates a new model of Model type.
     * @return the newly created empty model.
     */
    virtual AbstractModel* createModel () const;

    /**
     * Creates the model initialized with right.  This is like a copy
     * constructor.
     * @param right the model to copy.
     * @return the newly created copied model.
     */
    virtual AbstractModel* createModel (const AbstractModel &model) const;

    // I/O  -----------------------------------------------------------------

    /**
     * Writes the object to the output stream.
     * @param obs the output stream.
     * @return the written stream.
     */
    virtual oBinstream& write (oBinstream& obs) const;

  };
  

  /**
   * AnnotateModel
   * @author 
   * @version 
   */
  class AnnotateModel : public GraphModel
  {
    
  public:
    
    static const unsigned int MIN_HELIX_SIZE = 3;
    static const unsigned int PAIRING_MARK = 1;
    static const unsigned int LHELIX = 2;
    static const unsigned int RHELIX = 4;
    static const unsigned int weight_GC = 100;
    static const unsigned int weight_AU = 99;
    static const unsigned int weight_GU = 95;
    static const unsigned int weight_GA = 90;
    static const unsigned int weight_AC = 90;
    static const unsigned int weight_CU = 90;
    static const unsigned int weight_GG = 90;
    static const unsigned int weight_AA = 90;
    static const unsigned int weight_CC = 90;
    static const unsigned int weight_UU = 90;

  private:
    
    /**
     * The model name.
     */
    string name;

    /**
     * The model PdbFileHeader
     */
    PdbFileHeader fileHeader;

    /**
     * Collection of AnnotateModel sequences.
     */
    list< Sequence > sequences;
    list< Helix > helices;

//     StrandSet sequences;

//     int nb_pairings;

//     struct OStrand : public pair< label, label > {
//       stype type;
//       int ref;
//     };

    vector< BasePair > basepairs;
    vector< BaseStack > stacks;
    vector< BaseLink > links;
//     vector< OStrand > strands;
    
//     vector< int > sequence_length;

//     map< label, int > helix_mask;
//     map< label, int > strand_mask;
//     map< label, int > sequence_mask;
//     map< label, int > tertiary_mask;

    ResIdSet residueSelection;

    unsigned int environment;
        
  public:
    
    // LIFECYCLE ------------------------------------------------------------
    
    /**
     * Initializes the object.
     * @param rs a ResIdSet of residues to annotate.
     * @param env the number of relation layers around the residue selection to
     * annotate.
     */
    AnnotateModel (const ResIdSet &rs, unsigned int env)
      : GraphModel (new AnnotateResidueFM ()),
	residueSelection (rs),
	environment (env)
    { }
    
    /**
     * Initializes the object with the right's content (deep copy).
     * @param right the object to copy.
     * @param rs a ResIdSet of residues to annotate.
     * @param env the number of relation layers around the residue selection to
     * annotate.
     */
    AnnotateModel (const AbstractModel &right, const ResIdSet &rs, unsigned int env)
      : GraphModel (right, new AnnotateResidueFM ()),
	residueSelection (rs),
	environment (env)
    { }

    /**
     * Destroys the object.
     */
    virtual ~AnnotateModel () { }
    
    // OPERATORS ------------------------------------------------------------
    
    // ACCESS ---------------------------------------------------------------
    
    // METHODS --------------------------------------------------------------

  private:

    void s_secondaire (const vector< BasePair > &bps, set< BasePair > &stable);
    bool areHelix (const BasePair &bp1, const BasePair &bp2) const;

  public:
    
    /**
     * Builds the graph of relations, find strands and helices.
     */
    void annotate ();
    
  private :
    
    bool isHelixPairing (const Relation &r);

    bool isPairing (const Relation *r)
    {
      return r->is (PropertyType::pPairing);
    }

    bool isPairing (Residue *i, Residue *j)
    {
      return areConnected (i, j) && isPairing (getEdge (i, j));
    }
           
  public:
 
    void fillBPStacks ();
    void buildSequences ();
    void findHelices (const set< BasePair > &stable);
    void buildStrands();
    void findStrands ();
    void classifyStrands ();
    void findKissingHairpins ();
    void findPseudoknots ();

    // I/O  -----------------------------------------------------------------
  
    void dumpHelices () const;
    void dumpStrands ();
    void dumpSequences () const;
    void dumpPairs () const;
    void dumpConformations () const;
    void dumpTriples () ;
    void dumpStacks () const;

    /**
     * Ouputs the model to the stream.
     * @param os the output stream.
     * @return the used output stream.
     */
    virtual ostream& output (ostream &os) const;

    /**
     * Reads the model from a pdb input stream.
     * @param is the pdb data stream.
     * @return the consumed pdb stream.
     */
    virtual iPdbstream& input (iPdbstream &is);
  
//     /**
//      * Writes the model to a binary output stream.
//      * @param obs the binary data stream.
//      * @return the consumed binary stream.
//      */
//     virtual oBinstream& output (oBinstream &obs) const;

    /**
     * Reads the model from a binary input stream.
     * @param is the binary data stream.
     * @return the consumed binary stream.
     */
    virtual iBinstream& input (iBinstream &iss);

  };

//   /**
//    * Reads the AnnotateModel from a binary input stream.
//    * @param is the binary input stream.
//    * @param model the AnnotateModel.
//    * @return the consumed binary stream.
//    */
//   iBinstream& operator>> (iBinstream &is, AnnotateModel &model);

//   /**
//    * Writes the AnnotateModel to a binary output stream.
//    * @param os the binary input stream.
//    * @param model the AnnotateModel.
//    * @return the consumed binary stream.
//    */
//   oBinstream& operator<< (oBinstream &os, const AnnotateModel &model);

}



namespace std
{
  
  ostream &
  operator<< (ostream &out, const annotate::Strand &t);

  /**
   * Prints the AnnotateModel to the stream.
   * @param os the output stream.
   * @param am the AnnotateModel.
   * @return the output stream.
   */
  ostream& operator<< (ostream &os, const annotate::AnnotateModel &am);
  
}
  
#endif
