//                    -*- Mode: C++; coding: UTF-8 -*- 
// AnnotateResidue.cc
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mccore/Binstream.h"
#include "mccore/ResId.h"
#include "mccore/ResidueType.h"

#include "AnnotateResidue.h"



namespace mccore
{
  
  Residue*
  AnnotateResidueFM::createResidue () const
  {
    return new AnnotateResidue ();
  }


  Residue*
  AnnotateResidueFM::createResidue (const Residue &res) const
  {
    return new AnnotateResidue (res);
  }


  oBinstream&
  AnnotateResidueFM::write (oBinstream &obs) const
  {
    return obs << 'R';
  }

}
