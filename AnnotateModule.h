//                              -*- Mode: C++ -*- 
// AnnotateModule.h
// Copyright Å© 2004 Laboratoire de Biologie Informatique et ThÅÈorique
//                  UniversitÅÈ de MontrÅÈal.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Mon Dec  6 19:02:31 2004
// $Revision$


#ifndef _AnnotateModule_h_
#define _AnnotateModule_h_

#include <iostream>
#include <string>
#include <vector>

#include "mccore/ResIdSet.h"

#include "AnnotatedModel.h"
#include "Module.h"

using namespace mccore;
using namespace std;



namespace annotate
{
  /**
   * @author Martin Larose (<a href="mailto:larosem@iro.umontreal.ca">larosem@iro.umontreal.ca</a>)
   * @version $Id$
   */
  class AnnotateModule : public Module
  {
    /**
     * The residue selection.
     */
    ResIdSet selection;

    /**
     * The model collection.
     */
    vector< AnnotatedModel > models;

    // LIFECYCLE ------------------------------------------------------------
    
    /**
     * Initializes the module.
     */
    AnnotateModule () { }

    /**
     * Initializes the module with a copy.
     * @param right the module to copy.
     */
    AnnotateModule (const AnnotateModule &right)
      : selection (right.selection) { }
     
  public:

    /**
     * Initializes the module with the selection.
     * @param selection the residue selection.
     */
    AnnotateModule (const ResIdSet &selection)
      : selection (selection) { }

    /**
     * Destroys the module.
     */
    virtual ~AnnotateModule () { }

    // OPERATORS ------------------------------------------------------------
    
  private:

    /**
     * Copies the module right into this.
     * @param right the module to copy.
     * @return this.
     */
    AnnotateModule& operator= (const AnnotateModule &right)
    {
      if (this != &right)
	{
	  selection = right.selection;
	}
      return *this;
    }

    // ACCESS ---------------------------------------------------------------
    
    // METHODS --------------------------------------------------------------

    /**
     * Adds a model into the module.
     * @param model the model to insert.
     */
    virtual void add (Model &model) { }
    
  public:

    /**
     * Adds a model into the module.  Only the selected residues are kept.
     * @param name the name of the model.
     * @param header the header of the pdb model.
     * @param model the model to annotate.
     */
    void add (string &name, PdbFileHeader &header, Model &model);

    /**
     * Annotates every model in the module, finds helices and strands.
     */
    virtual void exec ();

    // I/O  -----------------------------------------------------------------

    /**
     * Writes the models into a MC-Sym script.
     * @param os the output stream.
     */
    void writeMCSym (ostream &os);

    /**
     * Writes the models into a RNAML file.
     * @param fileName the RNAML file name.
     */
    void writeRNAML (const string &fileName);
    
  };
}

#endif
