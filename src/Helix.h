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

      template< typename Iterator > Helix (Iterator start, Iterator end)
      : vector< BasePair > (start, end)
      { }

      ~Helix () { }

      // ACCESS ---------------------------------------------------------------
      unsigned int getId () const { return id; }

      void setId (unsigned int val) { id = val; }
  };
}

#endif
