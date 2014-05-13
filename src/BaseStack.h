//                           mcannotate
// Copyright   Laboratoire de Biologie Informatique et Theorique.
//                     Universite de Montreal

#ifndef _annotate_BaseStack_h_
#define _annotate_BaseStack_h_

#include <utility>

#include "mccore/GraphModel.h"
#include "mccore/ResId.h"

using namespace mccore;
using namespace std;


namespace annotate
{
  class BaseStack : public pair< GraphModel::label, GraphModel::label >
  {
    public:
    ResId fResId;
    ResId rResId;

    // LIFECYCLE ------------------------------------------------------------
    BaseStack (GraphModel::label l,
               const ResId &fResId,
               GraphModel::label r,
               const ResId &rResId)
    : pair< GraphModel::label, GraphModel::label > (l, r), fResId (fResId), rResId (rResId)
    { }

    ~BaseStack () { }

    // OPERATORS ------------------------------------------------------------
    BaseStack& operator= (const BaseStack &right)
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

    bool operator== (const BaseStack &right) const
    {
      return (&right == this
      || (first == right.first && second == right.second));
    }

    bool operator!= (const BaseStack &right) const
    {
      return ! operator== (right);
    }

    bool operator< (const BaseStack &right) const
    {
      return (&right != this
      && (fResId < right.fResId
      || (fResId == right.fResId
      && rResId < right.rResId)));
    }
  };
}

#endif
