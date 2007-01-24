//                    -*- Mode: C++; coding: UTF-8 -*- 
// Sequence.h
// Copyright © 2006-07 Laboratoire de Biologie Informatique et Théorique
//                     Université de Montréal.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Mon Nov 20 15:05:57 2006
// $Revision$
// $Id$
// 


#ifndef _annotate_Sequence_h_
#define _annotate_Sequence_h_

#include <iostream>
#include <vector>

#include "mccore/GraphModel.h"

using namespace mccore;
using namespace std;



namespace annotate
{
  class Sequence : public vector< GraphModel::label >
  {
  public:

    Sequence () : vector< GraphModel::label > () { }

    virtual ~Sequence () { }

  };
}

#endif
