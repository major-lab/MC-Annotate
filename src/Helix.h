//                           mcannotate
// Copyright   Laboratoire de Biologie Informatique et Theorique.
//                     Universite de Montreal

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
