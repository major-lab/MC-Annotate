//                              -*- Mode: C++ -*- 
// MolecularComplex.h
// Copyright © 2002, 2003 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Thu Mar  7 13:06:59 2002
// Last Modified By : Patrick Gendron
// Last Modified On : Tue Oct 28 11:34:16 2003
// Update Count     : 20
// Status           : Unknown.
// 


#ifndef _MolecularComplex_h_
#define _MolecularComplex_h_

#include <iostream>
#include <map>

#include "mccore/ResId.h"
#include "mccore/Residue.h"
#include "mccore/Pdbstream.h"
#include "mccore/PdbFileHeader.h"
#include "mccore/Model.h"

#include "MolecularUnit.h"

using namespace std;
using namespace mccore;


/**
 * @short A Model and annotation information.
 * @author Patrick Gendron
 */
class MolecularModel : public map< const char, MolecularUnit* >
{

public:

  /**
   * Initializes the object.
   */
  MolecularModel ();
  ~MolecularModel ();


  void annotate ();


  iPdbstream& input (iPdbstream &ips);
};

iPdbstream& operator>> (iPdbstream &ips, MolecularModel &obj);




// A PDB, check the name!!!
class MolecularComplex : public map< int, MolecularModel* >
{
  /**
   * Pdb file header.
   */
  PdbFileHeader header;

public:

  /**
   * Initializes the object.
   */
  MolecularComplex ();
  ~MolecularComplex ();

  iPdbstream& input (iPdbstream &ips);
  ostream& output (ostream& os) const;
};

iPdbstream& operator>> (iPdbstream &ips, MolecularComplex &obj);

ostream &operator<< (ostream &os, const MolecularComplex &obj);

#endif
