//                              -*- Mode: C++ -*- 
// MolecularUnit.cc
// Copyright © 2003 Laboratoire de Biologie Informatique et Théorique
//                  Université de Montréal.
// Author           : Patrick Gendron
// Created On       : Thu Oct 23 13:20:34 2003
// Last Modified By : Patrick Gendron
// Last Modified On : Wed Dec 17 16:16:05 2003
// Update Count     : 5
// Status           : Unknown.
// 



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "MolecularUnit.h"



MolecularUnit::MolecularUnit ()
{
  

}

MolecularUnit::~MolecularUnit ()
{
}


ostream &
MolecularUnit::output (ostream &os) const
{
  MolecularUnit::const_iterator z;
  for (z=begin (); z!=end (); ++z) {
    os << z->first << z->second->getType () << " ";   
  }
  return os;
}

ostream &operator<< (ostream &os, const MolecularUnit &obj)
{
  return obj.output (os);
}
