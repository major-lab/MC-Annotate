//                              -*- Mode: C++ -*- 
// View.h
// Copyright © 2004 Laboratoire de Biologie Informatique et Théorique
//                  Université de Montréal.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Mon Dec  6 18:16:37 2004
// $Revision$
// $Id$


#ifndef _View_h_
#define _View_h_



namespace annotate
{
  class Module;


  
  /**
   * @author Martin Larose (<a href="mailto:larosem@iro.umontreal.ca">larosem@iro.umontreal.ca</a>)
   * @version $Id$
   */
  class View
  {
  protected:
    
    /**
     * The module object.
     */
    Module *module;

    // LIFECYCLE ------------------------------------------------------------
    
    /**
     * Initializes the object.
     */
    View () : module (0) { }
    
    /**
     * Initializes the object with the other's content.
     * @param other the object to copy.
     */
    View (const View &other) : module (other.module) { }
    
  public:
    
    /**
     * Initializes the frame with a module.
     */
    View (Module *module) : module (module) { }
    
    /**
     * Destroys the object.
     */
    virtual ~View () { }
    
    // OPERATORS ------------------------------------------------------------
    
  private:
    
    /**
     * Assigns the object with the other's content.
     * @param other the object to copy.
     * @return itself.
     */
    View& operator= (const View &other)
    {
      if (this != &other)
	{
	  this->module = other.module;
	}
      return *this;
    }
    
    // ACCESS ---------------------------------------------------------------
    
    // METHODS --------------------------------------------------------------
    
  public:
    
    virtual void show () const = 0;
    
    // I/O  -----------------------------------------------------------------
    
  };
  
}

#endif
