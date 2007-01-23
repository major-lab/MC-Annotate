//                              -*- Mode: C++ -*- 
// MolecularComplex.cc
// Copyright © 2002, 2003 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Thu Mar  7 13:06:59 2002
// Last Modified By : Patrick Gendron
// Last Modified On : Tue Oct 28 11:46:49 2003
// Update Count     : 16
// Status           : Unknown.
// 



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "MolecularComplex.h"

#include "mccore/Residue.h"
#include "mccore/Model.h"




MolecularModel::MolecularModel ()
{
}


MolecularModel::~MolecularModel ()
{
}


void
MolecularModel::annotate () 
{

}



iPdbstream& MolecularModel::input (iPdbstream &ips)
{
  MolecularModel::iterator i;

  while (! (ips.eof ()))
    {
      Residue *res = new Residue;

      ips >> *res;

      if (res->size () != 0) {
	char chain = res->getResId ().getChainId ();

	i = find (chain);
	if (i!=end ())
	  (*(*i).second)[res->getResId ()] = res;
	else {
	  (*this)[chain] = new MolecularUnit;
	  (*(*this)[chain])[res->getResId ()] = res;
	}
      } else {
	delete res;
      }
      if (ips.eom ())
	break;
    }
  return ips;
}


iPdbstream& operator>> (iPdbstream &ips, MolecularModel &obj)
{
  return obj.input (ips);
}


MolecularComplex::MolecularComplex ()
{
}

MolecularComplex::~MolecularComplex ()
{
}

iPdbstream& MolecularComplex::input (iPdbstream &ips)
{
  int i = 1;

  header = ips.getHeader ();

  while (!ips.eof ())
    {
      MolecularModel *model = new MolecularModel;

      ips >> *model;
      
      if (!model->empty ())
	(*this)[i++] = model;
      else
	delete model;
    }
  return ips;
}


iPdbstream& operator>> (iPdbstream &ips, MolecularComplex &obj)
{
  return obj.input (ips);
}

ostream &
MolecularComplex::output (ostream &os) const
{
  MolecularComplex::const_iterator x;
  MolecularModel::iterator y;

  os << header << endl;
  for (x=begin (); x!=end (); ++x) {
    os << "model " << x->first << " : " << endl;
    for (y=x->second->begin (); y!=x->second->end (); ++y) {
      os << y->first << " -> " << *(y->second);
      
      os << endl;
    }
  }
  return os;
}

ostream &operator<< (ostream &os, const MolecularComplex &obj)
{
  return obj.output (os);
}
