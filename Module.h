//                              -*- Mode: C++ -*- 
// Module.h
// Copyright Å© 2004 Laboratoire de Biologie Informatique et ThÅÈorique
//                  UniversitÅÈ de MontrÅÈal.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Tue Dec  7 13:34:00 2004
// $Revision$


#ifndef _Module_h_
#define _Module_h_

namespace mccore
{
  class Model;
}



namespace annotate
{
  /**
   * @author Martin Larose (<a href="mailto:larosem@iro.umontreal.ca">larosem@iro.umontreal.ca</a>)
   * @version $Id$
   */
  class Module
  {
  protected:

    // LIFECYCLE ------------------------------------------------------------
    
    /**
     * Initializes the module.
     */
    Module () { }

    /**
     * Initializes the module with a copy.
     * @param right the module to copy.
     */
    Module (const Module &right) { }
     
    /**
     * Destroys the module.
     */
    virtual ~Module () { }

    // OPERATORS ------------------------------------------------------------
    
  private:

    /**
     * Copies the module right into this.
     * @param right the module to copy.
     * @return this.
     */
    Module& operator= (const Module &right) { return *this; }

    // ACCESS ---------------------------------------------------------------
    
    // METHODS --------------------------------------------------------------
    
  public:

    /**
     * Adds a model into the module.
     * @param model the model.
     */
    virtual void add (Model &model) = 0;

    /**
     * Prepares the algo for calculation.
     */
    virtual void exec () = 0;

    // I/O  -----------------------------------------------------------------

  };
}

#endif
