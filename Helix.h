//                    -*- Mode: C++; coding: UTF-8 -*- 
// Helix.h
// Copyright © 2006-07 Laboratoire de Biologie Informatique et Théorique
//                     Université de Montréal.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Wed Jul 12 13:57:45 2006
// $Revision$
// $Id$


#ifndef _annotate_Helix_h_
#define _annotate_Helix_h_

#include <vector>

#include "BasePair.h"

using namespace std;



namespace annotate
{

  class Helix : public vector< BasePair >
  {
    
    unsigned int id;

  public:

    // LIFECYCLE ------------------------------------------------------------

    Helix () : id (0) { }

    template< typename Iterator >
    Helix (Iterator start, Iterator end)
      : vector< BasePair > (start, end)
    { }

    ~Helix () { }

    // OPERATORS ------------------------------------------------------------

    // ACCESS ---------------------------------------------------------------

    unsigned int getId () const { return id; }

    void setId (unsigned int val) { id = val; }

    // METHODS --------------------------------------------------------------

    // I/O  -----------------------------------------------------------------
    
  };
  
}

#endif
