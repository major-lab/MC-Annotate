//                              -*- Mode: C++ -*- 
// AnnotatedModel.h
// Copyright © 2001-04 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// $Revision$


#ifndef _AnnotatedModel_h_
#define _AnnotatedModel_h_

#include <map>
#include <string>
#include <vector>
// #include <pdflib.h>

#include "mccore/Model.h"
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
  class AnnotatedModel;
  class mccore::GraphModel;
  
  typedef int strandId;
  
  enum stype { BULGE_OUT, BULGE, INTERNAL_LOOP, LOOP, HELIX, OTHER };
 
  class Strand : public vector< const Residue * >
  {
  
    private:
      char consensusChainId;
  
    public:
     
      virtual ostream&
      output (ostream &os) const
      {   
        Strand::const_iterator resIt;
        
        os << size() << " residues : ";
        for (resIt = begin (); end () != resIt; ++resIt)
        {
          os << (**resIt).getType();
          if(end() != resIt+1)
            os << "-";
        }
        
        return os;
      }
      
      virtual ~Strand() { }
      
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
   * AnnotatedModel
   * @author 
   * @version 
   */
  class AnnotatedModel
  {
    /**
     * The model name.
     */
    string name;

    /**
     * The model PdbFileHeader
     */
    PdbFileHeader fileHeader;

    GraphModel &gfm;
    StrandSet sequences;

    int nb_pairings;
    int nb_connect;
    int min_helix_size;

    struct Helix : public vector< pair<const Residue *, const Residue * > > {
    };
    struct OStrand : public pair< const Residue *, const Residue * > {
      stype type;
      int ref;
    };

    vector< Helix > helices;
    vector< OStrand > strands;
  
    vector< int > sequence_length;
    map< const Residue *, char > marks;

    map< const Residue *, int > helix_mask;
    map< const Residue *, int > strand_mask;
    map< const Residue *, int > sequence_mask;
    map< const Residue *, int > tertiary_mask;
        
  public:
    
    // LIFECYCLE ------------------------------------------------------------
    
    /**
     * Initializes the object.
     */
    AnnotatedModel (GraphModel &gfm);
    
    /**
     * Destroys the object.
     */
    virtual ~AnnotatedModel () { }
    
    // OPERATORS ------------------------------------------------------------
    
    // ACCESS ---------------------------------------------------------------
    
  private :
    
    // METHODS --------------------------------------------------------------

    bool isHelixPairing (const Relation *r) {
      return (r->is (PropertyType::pPairing));
//     return (relations[e].is (PropertyType::pXX) ||
//          relations[e].is (PropertyType::pXIX) ||
//          relations[e].is (PropertyType::pXXVIII));
    }

    bool isPairing (const Relation *r) {
      return (r->is (PropertyType::pPairing));
//     return (relations[e].is (PropertyType::pXX) ||
//          relations[e].is (PropertyType::pXIX) ||
//          relations[e].is (PropertyType::pXXVIII));
    }

  bool isPairing (Residue *i, Residue *j) {
    if (gfm.areConnected (i, j))
      return gfm.getEdge (i, j)->is (PropertyType::pPairing);
    return false;
  }
           
  public:
 
    void buildStrands(void);
    
    void findHelices (void);
    void dumpHelices (void) const;
    
    void findStrands (void);
    void classifyStrands (void);
    void dumpStrands (void);
    void findKissingHairpins (void);
    void findPseudoknots (void);

    void dumpSequences (bool detailed = true) ;
    void dumpPairs (void) ;
    void dumpConformations (void) ;
    void dumpTriples (void) ;
    void dumpStacks (void) ;

    // I/O  -----------------------------------------------------------------
  
    void 
    dumpMcc (const char* pdbname);
    
      
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
   * Reads the AnnotatedModel from a binary input stream.
   * @param is the binary input stream.
   * @param model the AnnotatedModel.
   * @return the consumed binary stream.
   */
  iBinstream& operator>> (iBinstream &is, AnnotatedModel &model);

  /**
   * Writes the AnnotatedModel to a binary output stream.
   * @param os the binary input stream.
   * @param model the AnnotatedModel.
   * @return the consumed binary stream.
   */
  oBinstream& operator<< (oBinstream &os, const AnnotatedModel &model);

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
  ostream& operator<< (ostream &os, const annotate::AnnotatedModel &am);
  
}
  
#endif
