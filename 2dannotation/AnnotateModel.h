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
// Temporary until cycle annotation is fully moved into AnnotationCycles.
#include "mccore/Molecule.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"
#include "mccore/ResId.h"
#include "mccore/ResIdSet.h"
#include "mccore/Residue.h"
#include "mccore/ResidueType.h"

namespace mccore
{
  class iBinstream;
  class iPdbstream;
}

namespace annotate
{
  class AnnotateModel;
  class Annotation;

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
  class AnnotateModelFM : public mccore::ModelFactoryMethod
  {
    /**
     * A ResIdSet of residues to annotate.
     */
    mccore::ResIdSet residueSelection;

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
    AnnotateModelFM (
    	const mccore::ResIdSet &rs,
    	unsigned int env,
    	const mccore::ResidueFactoryMethod *fm = 0)
	: mccore::ModelFactoryMethod (fm),
		residueSelection (rs),
		environment (env)
    { }

    /**
     * Initializes the object with the right content.
     * @param right the object to copy.
     */
    AnnotateModelFM (const mccore::ModelFM &right) : mccore::ModelFactoryMethod (right) { }

    /**
     * Clones the object.
     * @return the copy of the object.
     */
    virtual mccore::ModelFactoryMethod* clone () const
    {
      return (mccore::ModelFactoryMethod*) new AnnotateModelFM (*this);
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
    virtual mccore::AbstractModel* createModel () const;

    /**
     * Creates the model initialized with right.  This is like a copy
     * constructor.
     * @param right the model to copy.
     * @return the newly created copied model.
     */
    virtual mccore::AbstractModel* createModel (const mccore::AbstractModel &model) const;

    // I/O  -----------------------------------------------------------------

    /**
     * Writes the object to the output stream.
     * @param obs the output stream.
     * @return the written stream.
     */
    virtual mccore::oBinstream& write (mccore::oBinstream& obs) const;

  };


  /**
   * AnnotateModel
   * @author
   * @version
   */
  class AnnotateModel : public mccore::GraphModel
  {
    /**
     * The model name.
     */
    std::string mName;

    /**
     * Model number in the pdb file
     */
    unsigned int muiId;

    /**
     * The relation mask used for the annotation
     */
     unsigned char mucRelationMask;

    // TODO : Cycle annotation is affecting self, const correctness work needs to be done.
    mccore::Molecule mCyclesMolecule;

    mccore::ResIdSet residueSelection;

    unsigned int environment;

    std::vector<Annotation *> annotations;
    std::set<std::string> mProvidedAnnotations;

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
    AnnotateModel (
    	const mccore::ResIdSet &rs,
    	unsigned int env,
    	const mccore::ResidueFactoryMethod *fm = 0);

    /**
     * Initializes the object with the right's content (deep copy).
     * @param right the object to copy.
     * @param rs a ResIdSet of residues to annotate.
     * @param env the number of relation layers around the residue selection to
     * annotate.
     * @param fm the residue factory methods that will instanciate new
     * residues (default is @ref ExtendedResidueFM).
     */
    AnnotateModel (
    	const mccore::AbstractModel &right,
    	const mccore::ResIdSet &rs,
    	unsigned int env,
    	const mccore::ResidueFactoryMethod *fm = 0);

    /**
     * Destroys the object.
     */
    virtual ~AnnotateModel ();

    // OPERATORS ------------------------------------------------------------

    // ACCESS ---------------------------------------------------------------
	const Annotation* getAnnotation(const std::string& astrAnnotName) const;

    template<class T>
    const T* getAnnotation() const
    {
    	return dynamic_cast<const T*>(getAnnotation(T::AnnotationName()));
    }

    const mccore::Molecule& getCyclesMolecule() const
    {
    	return mCyclesMolecule;
    }

    const std::string& name() const {return mName;}
    void name(const std::string& aName) {mName = aName;}

    // id of the model in the pdb file
    unsigned int id() const {return muiId;}
    void id(unsigned int auiId) {muiId = auiId;}

    unsigned char relationMask() const {return mucRelationMask;}

    // METHODS --------------------------------------------------------------

    /**
     * Builds the graph of relations, find strands and helices.
     */
    void annotate (unsigned char aspb =
    	mccore::Relation::adjacent_mask
    	|mccore::Relation::pairing_mask
    	|mccore::Relation::stacking_mask
    	|mccore::Relation::backbone_mask);

    /**
     * Computes the union of minimum cycle bases.
     */
    void computeUnionMinimumCycleBases();

    void addAnnotation(Annotation& aAnnotation);

  private :

    bool isPairing (const mccore::Relation *r)
    {
      return r->is (mccore::PropertyType::pPairing);
    }

    bool isPairing (mccore::Residue *i, mccore::Residue *j)
    {
      return areConnected (i, j) && isPairing (getEdge (i, j));
    }

  public:

    void dumpConformations () const;

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
    virtual mccore::iPdbstream& input (mccore::iPdbstream &is);

    /**
     * Reads the model from a binary input stream.
     * @param is the binary data stream.
     * @return the consumed binary stream.
     */
    virtual mccore::iBinstream& input (mccore::iBinstream &iss);

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
