//                              -*- Mode: C++ -*- 
// MolecularUnit.h
// Copyright © 2003 Laboratoire de Biologie Informatique et Théorique
//                  Université de Montréal.
// Author           : Patrick Gendron
// Created On       : Thu Oct 23 13:20:34 2003
// Last Modified By : Patrick Gendron
// Last Modified On : Mon Oct 27 10:49:39 2003
// Update Count     : 9
// Status           : Unknown.
// 


#ifndef _MolecularUnit_h_
#define _MolecularUnit_h_

#include <iostream>

#include "mccore/Residue.h"
#include "mccore/Pdbstream.h"

using namespace std;
using namespace mccore;

/**
 * @short A sub model with a single chain ID
 *
 * Long Description
 *
 * @author Patrick Gendron
 */
class MolecularUnit : public map< const ResId, Residue* >
{

public:

  // LIFECYCLE ------------------------------------------------------------

  /**
   * Initializes the object.
   */
  MolecularUnit ();

  /**
   * Initializes the object with the other's content.
   * @param other the object to copy.
   */
  MolecularUnit (const MolecularUnit &other);

  /**
   * Destroys the object.
   */
  ~MolecularUnit ();

  // OPERATORS ------------------------------------------------------------

  /**
   * Assigns the object with the other's content.
   * @param other the object to copy.
   * @return itself.
   */
  MolecularUnit& operator= (const MolecularUnit &other);

  // ACCESS ---------------------------------------------------------------

  // METHODS --------------------------------------------------------------

  // I/O  -----------------------------------------------------------------

  iPdbstream& input (iPdbstream &ips);
  ostream& output (ostream& os) const;
};

iPdbstream& operator>> (iPdbstream &ips, MolecularUnit &obj);

ostream &operator<< (ostream &os, const MolecularUnit &obj);

#endif
