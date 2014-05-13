//                           mcannotate
// Copyright   Laboratoire de Biologie Informatique et Theorique.
//                     Universite de Montreal

#ifndef _annotate_BaseLink_h_
#define _annotate_BaseLink_h_

#include <utility>

#include "mccore/GraphModel.h"
#include "mccore/ResId.h"

using namespace mccore;
using namespace std;


namespace annotate
{
  class BaseLink : public pair< GraphModel::label, GraphModel::label >
  {
  public:
    ResId fResId;
    ResId rResId;

    // LIFECYCLE ------------------------------------------------------------
    BaseLink (GraphModel::label l,
              const ResId &fResId,
              GraphModel::label r,
              const ResId &rResId)
    : pair< GraphModel::label, GraphModel::label > (l, r), fResId (fResId), rResId (rResId)
    { }

    ~BaseLink () { }

    // OPERATORS ------------------------------------------------------------
    BaseLink& operator= (const BaseLink &right)
    {
      if (this != &right)
      {
        first = right.first;
        second = right.second;
        fResId = right.fResId;
        rResId = right.rResId;
      }
      return *this;
    }

    bool operator== (const BaseLink &right) const
    {
      return (&right == this || (first == right.first && second == right.second));
    }

    bool operator!= (const BaseLink &right) const
    {
      return ! operator== (right);
    }

    bool operator< (const BaseLink &right) const
    {
      return (&right != this
      && (fResId < right.fResId
      || (fResId == right.fResId && rResId < right.rResId)));
    }
  };
}

#endif
