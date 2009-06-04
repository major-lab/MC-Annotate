//                              -*- Mode: C++ -*- 
// AnnotateModel.h
// Copyright © 2001-06 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// $Revision: 58 $
// $Id: AnnotateModel.h 58 2006-11-15 21:09:19Z larosem $


#ifndef _annotate_AnnotateModel_h_
#define _annotate_AnnotateModel_h_

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "mccore/GraphModel.h"
#include "mccore/Model.h"
#include "mccore/ModelFactoryMethod.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"
#include "mccore/ResId.h"
#include "mccore/ResIdSet.h"
#include "mccore/Residue.h"
#include "mccore/ResidueType.h"

#include "BaseLink.h"
#include "BasePair.h"
#include "BaseStack.h"
#include "Stem.h"
#include "Loop.h"
#include "Linker.h"

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

  /**
   * @short ModelFactoryMethod implementation for AnnotateModel class.
   *
   * This is the model factory method implementation for the AnnotateModel
   * class.
   *
   * @author Martin Larose (<a href="larosem@iro.umontreal.ca">larosem@iro.umontreal.ca</a>)
   * @version $Id: AnnotateModel.h 58 2006-11-15 21:09:19Z larosem $
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
    AnnotateModelFM (const ResIdSet &rs, unsigned int env, const ResidueFactoryMethod *fm = 0)
      : ModelFactoryMethod (fm),
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
    /**
     * The model name.
     */
    string name;

    std::vector< std::vector< const Residue * > > chains;

    vector< BasePair > basepairs;
    vector< BaseStack > stacks;
    vector< BaseLink > links;
    std::vector< Stem > stems;
    std::vector< Loop > loops;
    std::vector< Linker > linkers;

    vector< unsigned int > marks;

    ResIdSet residueSelection;

    unsigned int environment;
        
  public:
    
    // LIFECYCLE ------------------------------------------------------------
    
    /**
     * Initializes the object.
     * @param rs a ResIdSet of residues to annotate.
     * @param env the number of relation layers around the residue selection to
     * annotate.
     * @param fm the residue factory methods that will instanciate new
     * residues (default is @ref ExtendedResidueFM).
     */
    AnnotateModel (const ResIdSet &rs, unsigned int env, const ResidueFactoryMethod *fm = 0)
      : GraphModel (fm),
	residueSelection (rs),
	environment (env)
    { }
    
    /**
     * Initializes the object with the right's content (deep copy).
     * @param right the object to copy.
     * @param rs a ResIdSet of residues to annotate.
     * @param env the number of relation layers around the residue selection to
     * annotate.
     * @param fm the residue factory methods that will instanciate new
     * residues (default is @ref ExtendedResidueFM).
     */
    AnnotateModel (const AbstractModel &right, const ResIdSet &rs, unsigned int env, const ResidueFactoryMethod *fm = 0)
      : GraphModel (right, fm),
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

    /**
     * Builds the graph of relations, find strands and helices.
     */
    void annotate ();
    
  private :
  
  	std::set< BasePair > getWWBasePairs();
    
    bool isPairing (const Relation *r)
    {
      return r->is (PropertyType::pPairing);
    }

    bool isPairing (Residue *i, Residue *j)
    {
      return areConnected (i, j) && isPairing (getEdge (i, j));
    }
    
    void dumpPair (const BasePair& aBasePair) const;
    void dumpLinker(const Linker& aLinker) const;
    void dumpLoop (const Loop& aLoop) const;
    bool enclose(const BasePair& aBasePair, const Stem& aStem);
    std::vector<const Stem*> getEnclosedStems(const BasePair& aBasePair);
    std::map< const Residue*, const Stem* > getResidueStemAssociation(
    	unsigned int iChain) const;
    int getDirection(
    	const Stem& aStem,		// Starting stem
		const ResId& aResId);	// Residue from which we're starting
	int getDirection(const StemConnection& aConnections) const;
		
	Loop findLoop(const Stem& aStem);
	Linker findLinker(const StemConnection& aConnection) const;
	void findLinker(
		const Stem* apStem, 
		const Stem::enConnection& aeConnect,
  		std::set<Linker>& outLinkerSet);
	
	struct stResidueInfo
	{
		ResId resId;
		const mccore::Residue* pResidue;
		const Stem* pStem;
	};
	std::vector<stResidueInfo> mResidueInfos;
	void computeResidueInfos();
	std::vector<stResidueInfo>::const_iterator findResidueInfo(
		const ResId& aResId) const;
	mccore::ResId nextId(
		const Stem& aStem, 
		const StemConnection& aConnection) const;
	std::map<mccore::ResId, const Linker*> getResidueLinkerMap() const;
    Linker nextLinker(
		const Linker& aLinker,
		const std::map<mccore::ResId, const Linker*>& aResidueLinkerMap);
	void removeLinker(
		std::map<mccore::ResId, const Linker*>& aResidueLinkerMap,
		const Linker& aLinker);
  public:
 
 	void findStems ();
 	void dumpStems () const;
 	
 	void findLinkers ();
 	void dumpLinkers () const;
 	
	void findLoops();
 	void dumpLoops() const;
 	
 	void findChains();
 	void dumpChains () const;

    void fillSeqBPStacks ();

    void dumpPairs () const;
    void dumpConformations () const;
    void dumpStacks () const;

    // I/O  -----------------------------------------------------------------
  
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
  
    /**
     * Reads the model from a binary input stream.
     * @param is the binary data stream.
     * @return the consumed binary stream.
     */
    virtual iBinstream& input (iBinstream &iss);

  };
}



namespace std
{
   /**
   * Prints the AnnotateModel to the stream.
   * @param os the output stream.
   * @param am the AnnotateModel.
   * @return the output stream.
   */
  ostream& operator<< (ostream &os, const annotate::AnnotateModel &am);
  
}
  
#endif
