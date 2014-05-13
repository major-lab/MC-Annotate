//    Copyright 2006 Laboratoire de Biologie Informatique et Theorique
// 
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
// 
//        http://www.apache.org/licenses/LICENSE-2.0
// 
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

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
