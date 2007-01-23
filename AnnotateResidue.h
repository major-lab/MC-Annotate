//                    -*- Mode: C++; coding: UTF-8 -*- 
// AnnotateResidue.h
// Copyright © 2007 Laboratoire de Biologie Informatique et Théorique.
//                  Université de Montréal
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Mon Jan 22 17:24:34 2007
// $Revision$
// 
// This file is part of mccore.
// 
// mccore is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// 
// mccore is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with mccore; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#ifndef _mccore_AnnotateResidue_h_
#define _mccore_AnnotateResidue_h_

#include <vector>

#include "mccore/Atom.h"
#include "mccore/Residue.h"
#include "mccore/ResidueFactoryMethod.h"

using namespace std;



namespace mccore
{
  class oBinstream;
  class ResId;
  class ResidueType;


  
  /**
   * This class extends the Residue class from mccore. It adds a properties
   * member field to encode annotations.
   *
   * @author Martin Larose (<a href="larosem@iro.umontreal.ca">larosem@iro.umontreal.ca</a>)
   * @version $Id$
   */
  class AnnotateResidue : public Residue
  {
  public:

    static const unsigned int LPAREN = 1;
    static const unsigned int RPAREN = 2;
    static const unsigned int X = 4;
    static const unsigned int B = 8;

  protected:
    
    /**
     * Bit vector containing annotation properties.
     * 1 = (
     * 2 = )
     * 4 = X
     * 8 = b
     */
    unsigned int properties;

  public:

    /**
     * Initializes the AnnotateResidue.
     */
    AnnotateResidue () : Residue (), properties (0) { }

    /**
     * Initializes the AnnotateResidue with a type and id.
     * @param t the residue type.
     * @param i the residue id.
     */
    AnnotateResidue (const ResidueType *t, const ResId &i)
      : Residue (t, i),
	properties (0)
    { }
    
    /**
     * Initializes the AnnotateResidue with type, atom container and id.
     * @param type the residue type.
     * @param i the residue id.
     * @param vec the atom container.
     */
    AnnotateResidue (const ResidueType *t, const ResId &i, vector< Atom > &vec)
      : Residue (t, i, vec),
	properties (0)
    { }
    
    /**
     * Initializes this object's content with another's.
     * @param anres the other object from which to copy content.
     */
    AnnotateResidue (const AnnotateResidue &anres)
      : Residue (anres),
	properties (anres.properties)
    { }

    /**
     * Initializes this object's content with another's by resolving
     * its polymorphic type.
     * @param res the polymorphic object from which to copy content.
     */
    AnnotateResidue (const Residue &res)
      : Residue (res),
	properties (0)
    { }

    /**
     * Clones the residue.
     * @return the copy of the object.
     */
    virtual Residue* clone () const { return new AnnotateResidue (*this); }
    
    /**
     * AnnotateResidue destructor.
     */
    virtual ~AnnotateResidue () { }

    // OPERATORS ------------------------------------------------------------

    /**
     * Assigns this object's content with another's. Kept for compatibility
     * reason.
     * @param anres the other object from which to copy content.
     * @return itself.
     */
    AnnotateResidue& operator= (const AnnotateResidue &anres)
    {
      if (&anres != this)
	{
	  Residue::operator= ((const Residue&) anres);
	  properties = anres.properties;
	}
      return *this;
    }

    /**
     * Assigns this object's content with another's by resolving its
     * polymorphic type. Kept for compatibility reason.
     * @param res the polymorphic object from which to copy content.
     * @return itself.
     */
    AnnotateResidue& operator= (const Residue& res)
    {
      if (&res != this)
	{
	  Residue::operator= (res);
	  properties = 0;
	}
      return *this;
    }

    // ACCESS ---------------------------------------------------------------

    /**
     * Gets the annotation properties.
     * @return the annotation properties.
     */
    unsigned int getProperties () const { return properties; }

  };


  
  /**
   * This is the residue factory method implementation for the AnnotateResidue
   * class.
   *
   * @author Martin Larose (<a href="larosem@iro.umontreal.ca">larosem@iro.umontreal.ca</a>)
   * @version $Id$
   */
  class AnnotateResidueFM : public ResidueFactoryMethod
  {

  public:

    // LIFECYCLE ------------------------------------------------------------

    /**
     * Initializes the object.
     */
    AnnotateResidueFM () { }

    /**
     * Clones the object.
     * @return the copy of the object.
     */
    virtual ResidueFactoryMethod* clone () const
    {
      return new AnnotateResidueFM ();
    }
  
    /**
     * Destroys the object.
     */
    virtual ~AnnotateResidueFM () { }

    // OPERATORS ------------------------------------------------------------

    // ACCESS ---------------------------------------------------------------

    // METHODS --------------------------------------------------------------

    /**
     * Creates a new residue of Residue type.
     * @return the newly created empty residue.
     */
    virtual Residue* createResidue () const;

    /**
     * Creates a residue copy.
     * @param res the residue to copy from.
     * @return the newly created residue copy.
     */
    virtual Residue* createResidue (const Residue &res) const;

    // I/O  -----------------------------------------------------------------

    /**
     * Writes the object to the output stream.
     * @param obs the output stream.
     * @return the written stream.
     */
    virtual oBinstream& write (oBinstream &obs) const;

  };

}

#endif
