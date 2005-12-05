//                              -*- Mode: C++ -*- 
// AnnotateModel.h
// Copyright � 2001-05 Laboratoire de Biologie Informatique et Th�orique.
//                     Universit� de Montr�al
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// $Revision$
// $Id$


#ifndef _annotate_AnnotateModel_h_
#define _annotate_AnnotateModel_h_

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
// #include <pdflib.h>

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

using namespace std;
using namespace mccore;



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

  public:

    // LIFECYCLE ------------------------------------------------------------

    /**
     * Initializes the object.
     */
    AnnotateModelFM (const ResidueFactoryMethod *fm = 0)
      : ModelFactoryMethod (fm) { }

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
    /**
     * The model name.
     */
    string name;

    /**
     * The model PdbFileHeader
     */
    PdbFileHeader fileHeader;

//     GraphModel &gfm;
    StrandSet sequences;

    int nb_pairings;
    int min_helix_size;

    typedef vector< pair< label, label > > Helix;
    struct OStrand : public pair< label, label > {
      stype type;
      int ref;
    };

    vector< Helix > helices;
    vector< OStrand > strands;
  
    vector< int > sequence_length;
    map< label, char > marks;

    map< label, int > helix_mask;
    map< label, int > strand_mask;
    map< label, int > sequence_mask;
    map< label, int > tertiary_mask;
        
  public:
    
    // LIFECYCLE ------------------------------------------------------------
    
    /**
     * Initializes the object.
     * @param fm the residue factory methods that will instanciate new
     * residues (default is @ref ExtendedResidueFM).
     */
    AnnotateModel (const ResidueFactoryMethod *fm = 0)
      : GraphModel (fm)
    { }
    
    /**
     * Initializes the object with the right's content (deep copy).
     * @param right the object to copy.
     * @param fm the residue factory methods that will instanciate new
     * residues (default is @ref ExtendedResidueFM).
     */
    AnnotateModel (const AbstractModel &right, const ResidueFactoryMethod *fm = 0)
      : GraphModel (right, fm)
    { }

    /**
     * Initializes the object with the right's content (deep copy).
     * @param right the object to copy.
     * @param fm the residue factory methods that will instanciate new
     * residues (default is @ref ExtendedResidueFM).
     */
    AnnotateModel (const GraphModel &right, const ResidueFactoryMethod *fm = 0)
      : GraphModel (right, fm)
    { }

//     /**
//      * Initializes the object with the right's content (deep copy).
//      * @param right the object to copy.
//      * @param fm the residue factory methods that will instanciate new
//      * residues (default is @ref ExtendedResidueFM).
//      */
//     AnnotateModel (const AnnotateModel &gfm, const ResidueFactoryMethod *fm = 0);
    
    /**
     * Destroys the object.
     */
    virtual ~AnnotateModel () { }
    
    // OPERATORS ------------------------------------------------------------
    
    // ACCESS ---------------------------------------------------------------
    
    // METHODS --------------------------------------------------------------

    /**
     * Builds the graph of relations, find strands and helices.
     */
    void annotate ();
    
  private :
    
    bool isHelixPairing (const Relation &r);

    bool isPairing (const Relation *r)
    {
      return (r->is (PropertyType::pPairing));
    }

    bool isPairing (Residue *i, Residue *j) {
      if (areConnected (i, j))
	{
	  return getEdge (i, j)->is (PropertyType::pPairing);
	}
      return false;
    }
           
  public:
 
    void buildStrands();
    
    void findHelices (const set< pair< label, label > > &helixPairsCandidates);
    void dumpHelices () const;
    
    void findStrands ();
    void classifyStrands ();
    void dumpStrands ();
    void findKissingHairpins ();
    void findPseudoknots ();

    void dumpSequences (bool detailed = true) ;
    void dumpPairs () const;
    void dumpConformations () ;
    void dumpTriples () ;
    void dumpStacks () const;

    // I/O  -----------------------------------------------------------------
  
    void dumpMcc (const char* pdbname);
    
      
    /**
     * Ouputs the model to the stream.
     * @param os the output stream.
     * @return the used output stream.
     */
    virtual ostream& output (ostream &os) const;

    /**
     * Reads the model from a pdb input stream.
     * @param ips the pdb data stream.
     * @return the consumed pdb stream.
     */
    virtual iPdbstream& input (iPdbstream &ips);
  
    /**
     * Writes the model to a binary output stream.
     * @param obs the binary data stream.
     * @return the consumed binary stream.
     */
    virtual oBinstream& output (oBinstream &obs) const;

    /**
     * Reads the model from a binary input stream.
     * @param obs the binary data stream.
     * @return the consumed binary stream.
     */
    virtual iBinstream& input (iBinstream &ibs);

  };

  /**
   * Reads the AnnotateModel from a binary input stream.
   * @param is the binary input stream.
   * @param model the AnnotateModel.
   * @return the consumed binary stream.
   */
  iBinstream& operator>> (iBinstream &is, AnnotateModel &model);

  /**
   * Writes the AnnotateModel to a binary output stream.
   * @param os the binary input stream.
   * @param model the AnnotateModel.
   * @return the consumed binary stream.
   */
  oBinstream& operator<< (oBinstream &os, const AnnotateModel &model);

}

namespace std
{
  
  ostream &
  operator<< (ostream &out, const annotate::Strand &t);

  /**
   * Ouputs the residue to the stream.
   * @param os the output stream.
   * @param r the residue.
   * @return the used output stream.
   */
  ostream& operator<< (ostream &os, const annotate::AnnotateModel &am);
  
}
  
#endif
