//                              -*- Mode: C++ -*- 
// TerminalView.h
// Copyright © 2004 Laboratoire de Biologie Informatique et Théorique
//                  Université de Montréal.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Mon Dec  6 18:11:47 2004
// $Revision$
// $Id$


#ifndef _TerminalView_h_
#define _TerminalView_h_

#include "View.h"

namespace mccore
{
  class Messagestream;
}

using namespace mccore;



namespace annotate
{
  class Module;


  
  /**
   * @author Martin Larose (<a href="mailto:larosem@iro.umontreal.ca">larosem@iro.umontreal.ca</a>)
   * @version $Id$
   */
  class TerminalView : public View
  {
    /**
     * The output message stream.
     */
    Messagestream *ms;
    
  public:

    // LIFECYCLE ------------------------------------------------------------

    /**
     * Initializes the object.
     */
    TerminalView (Module *module, Messagestream &ms)
      : View (module), ms (&ms) { }

  private:
    
    /**
     * Initializes the object with the right's content.
     * @param right the object to copy.
     */
    TerminalView (const TerminalView& right) : View (right), ms (right.ms) { }

    /**
     * Destructs the object.
     */
    virtual ~TerminalView () { }

    // OPERATORS ------------------------------------------------------------

  private:
    
    /**
     * Assigns the object with the right's content.
     * @param right the object to copy.
     * @return itself.
     */
    TerminalView& operator= (const TerminalView& right) { return *this; }

  public:

    // ACCESS ---------------------------------------------------------------

    // METHODS --------------------------------------------------------------

    virtual void show () const;
    
    // I/O  -----------------------------------------------------------------

  };
}

#endif
