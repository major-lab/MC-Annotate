//                              -*- Mode: C++ -*- 
// AnnotateModule.cc
// Copyright © 2004 Laboratoire de Biologie Informatique et Théorique.
//                  Université de Montréal
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Fri Dec  3 16:31:38 2004
// $Revision$
// $Id$


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>

#include "mccore/GraphModel.h"
#include "mccore/PdbFileHeader.h"
#include "mccore/ResId.h"
#include "mccore/ResIdSet.h"
#include "mccore/Residue.h"

#include "AnnotateModule.h"



namespace annotate
{

  void
  AnnotateModule::add (string &name, PdbFileHeader &header, GraphModel &model)
  {
    ResIdSet modelSet;
    ResIdSet tmp;
    GraphModel::iterator mit;
    ResIdSet::iterator sit;
    
    gOut (3) << "selection = " << selection << endl;
    gOut (3) << " (" << model.size () << " residues) " << flush;

    model.sort ();
    for (mit = model.begin (); model.end () != mit; ++mit)
      {
	modelSet.insert (mit->getResId ());
      }
    set_difference (modelSet.begin (), modelSet.end (),
		    selection.begin (), selection.end (),
		    inserter (tmp, tmp.begin ()));
    for (sit = tmp.begin (); tmp.end () != sit; ++sit)
      {
	model.erase (model.find (*sit));
      }
    if (! model.empty ())
      {
	models.push_back (AnnotatedModel (name, header, model));
      }
  }


  void
  AnnotateModule::exec ()
  {
    vector< AnnotatedModel >::iterator it;

    for (it = models.begin (); models.end () != it; ++it)
      {
	gOut (3) << "Annotating " << it->getName () << endl;
	it->annotate ();

	gOut (3) << "Finding helices..." << flush;
	it->findHelices ();
	gOut (3) << "done." << endl;
	gOut (3) << "Finding strands..." << flush;
	it->findStrands ();
	it->classifyStrands ();
	gOut (3) << "done." << endl;
      }
  }


  void
  AnnotateModule::writeMCSym (ostream &os)
  {
    vector< AnnotatedModel >::iterator it;

    for (it = models.begin (); models.end () != it; ++it)
      {
	it->dumpMcc (os);
	os << endl;
      }
  }

  
  void
  AnnotateModule::writeRNAML (const string &fileName)
  {
  }
}
