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
#include "mccore/UndirectedGraph.h"

using namespace std;
using namespace mccore;



namespace annotate
{
  class AnnotatedModel;
//   class Block;
  
  typedef int node;
  typedef int edge;
  typedef int strandid;
  
  enum stype { BULGE_OUT, BULGE, INTERNAL_LOOP, LOOP, HELIX, OTHER };
  
  /**
   * @short Helix
   * @author Patrick Gendron
   * -----------------------------------------------------------------------------
   */
  struct Helix : public vector< pair< node, node > > {
  };
  
  
  /**
   * @short Strand
   * @author Patrick Gendron
   * -----------------------------------------------------------------------------
   */
  struct Strand : public pair< node, node > {
    stype type;
    strandid ref;
  };
  

  
//   /**
//    * @short Vector
//    * @author Patrick Gendron
//    * -----------------------------------------------------------------------------
//  */
//   struct Vector
//   {
//     float x;
//     float y;
    
//     Vector () { }
//     Vector (float _x, float _y) : x (_x), y (_y) { }
    
//     Vector operator+ (const Vector& v) { return Vector (x+v.x, y+v.y); }
//     float operator* (const Vector& v) { return x*v.x + y*v.y; }
//     float operator| (const Vector& v) { return sqrt ((x-v.x)*(x-v.x) + (y-v.y)*(y-v.y)); }
//     Vector operator* (float f) { return Vector (f*x, f*y); }
//     Vector operator/ (float f) { return Vector (x/f, y/f); }
    
//     float norm () { return sqrt (x*x + y*y); }
    
//     friend ostream& operator<< (ostream &os, const Vector &obj) 
//     {
//       os << "<" << obj.x << ", " << obj.y << ">" << endl;
//       return os;
//     }
//   };
  
  
  /**
   * AnnotatedModel
   * @author Patrick Gendron (<a href="gendrop@iro.umontreal.ca">gendrop@iro.umontreal.ca</a>)
   * @version $Id$
   */
  class AnnotatedModel : public UndirectedGraph< Residue, Relation >
  {
    /**
     * The model name.
     */
    string name;

    /**
     * The model PdbFileHeader
     */
    PdbFileHeader fileHeader;

    /**
     * Map from integer to original ResId.  The ResId of the residues are
     * replaced with sequential ids.
     */
    map< unsigned int, ResId > int2ResIdMap;

//     vector< Model::iterator > conformations;
    
    vector< int > sequence_mask;
    vector< int > sequence_length;
    
    vector< int > helix_mask;
    vector< Helix > helices;
    
    vector< int > strand_mask;
    vector< Strand > strands;
    
    vector< int > tertiary_mask;
    
//     /**
//      * Relations
//      */
//     vector< Relation > relations;
    
//     UndirectedGraph< node, edge > graph;
    
    vector< char > marks;
    
    map< ResId, ResId > translation;
    
    int nb_pairings;
    int nb_connect;
    
    //   /**
    //    * Blocks...
    //    */
    //   vector< Block > h_blocks;
    //   vector< Block > l_blocks;
    //   vector< Block > s_blocks;
    
  public:
    
    // LIFECYCLE ------------------------------------------------------------
    
    /**
     * Initializes the object.
     */
    AnnotatedModel (string &name, PdbFileHeader &header, Model &m)
      : name (name), fileHeader (header), nb_pairings (0), nb_connect (0)
    {
      insert (m.begin (), m.end ());
    }
    
    /**
     * Destroys the object.
     */
    ~AnnotatedModel () { }
    
    // OPERATORS ------------------------------------------------------------
    
    // ACCESS ---------------------------------------------------------------
    
  private :
    
    ResId& getResId (node i)
    {
      return translation[conformations[i]->getResId ()]; 
    }
    
    ResId & getRefId (node i, node j)
    {
      return translation[relations[graph.getEdge(i, j)].getRef ()->getResId ()];
    }

    ResId & getResId (node i, node j)
    {
      return translation[relations[graph.getEdge(i, j)].getRes ()->getResId ()];
    }
    
    const ResidueType* getType (node i)
    {
      return conformations[i]->getType ();
    }

    bool isPairing (edge e)
    {
      return relations[e].is (PropertyType::pPairing);
    }
  
    bool isHelixPairing (edge e)
    {
      return (relations[e].is (PropertyType::pPairing));
    }
    
    bool isAdjacent (edge e)
    {
      return relations[e].is (PropertyType::pAdjacent);
    }
    
    bool isStacking (edge e)
    {
      return relations[e].is (PropertyType::pStack);
    }
    
    bool isPairing (node i, node j)
    {
      if (graph.areConnected (i, j))
	return relations[graph.getEdge (i, j)].is (PropertyType::pPairing);
      return false;
    }
    
    // METHODS --------------------------------------------------------------
    
  public:
    
    /**
     * Creates the annotation and determines the sequences present in the 
     * structure.
     */
    void annotate ();
    
    /**
     * Finds all helices.
     */
    void findHelices ();

    /**
     *
     */
    void findStrands ();
    
    /**
     *
     */
    void classifyStrands ();
    
    /**
     *
     */
    void findKissingHairpins ();
    
    /**
     * 
     */
    void findPseudoknots ();
    
    /**
     * Extracts...
     */
    ResIdSet extract (ResIdSet &seed, int size);
    
    // I/O  -----------------------------------------------------------------
    
    void dumpSequences (bool detailed=true);
    
    void dumpGraph ();
    void dumpConformations ();
    void dumpPairs ();
    
    void dumpHelices ();
    void dumpTriples ();
    void dumpStrands ();
    void dumpStacks ();
    
    void dumpMcc (const char* pdbname);
    
    //   void dumpCt (const char* pdbname);
    
    //   void PDF_drawLoop2 (PDF *p, int li, int source_jct);
    //   void PDF_drawHelix (PDF *p, int hi, int li_ref = -1, int x = -1, int y = -1);
    //   void PDF_drawLoop (PDF *p, int li, int hi_ref = -1, int x = -1, int y = -1);
    //   void dumpPDF (const char* pdfname);
    
    //   void PDF_drawGraph (PDF *p, node i, node prev, int x, int y, int a);
    //   void dumpSimplePDF (const char* pdbname, const char* pdfname);
};
  
  
  // /**
  //  * @short Block
  //  * @author Patrick Gendron
  //  * -----------------------------------------------------------------------------
  //  */
  // struct Block : public vector< vector< node > > {
  //   AnnotatedModel *amodel;
  //   vector< pair< node, node > > junctions;
  
  //   Block (AnnotatedModel* m) : amodel (m) {}
  
  //   ostream& output (ostream &os) const {
  //     vector< pair< node, node > >::const_iterator k;
  //     vector< vector< node > >::const_iterator i;
  //     vector< node >::const_iterator j;
  
  //     for (k=junctions.begin (); k!=junctions.end (); ++k) {
  //       os << "(" << amodel->getResId (k->first) << ", " 
  // 	   << amodel->getResId (k->second) << ")" << flush;
  //     }
  //     os << " " << flush;
  //     for (i=begin (); i!=end (); ++i) {
  //       os << "[";
  //       for (j=i->begin (); j!=i->end (); ++j) {
  // 	if (j!=i->begin ()) os << ", ";
  // 	os << amodel->getResId (*j);
  //       }
  //       os << "] ";
  //     }
  //     return os;
  //   }
  
  //   friend ostream& operator<< (ostream &os, const Block &b) {
  //     return b.output (os);
  //   }
  // };

#endif
