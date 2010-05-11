//                              -*- Mode: C++ -*-
// BasePair.cc
// Copyright Â© 2006 Laboratoire de Biologie Informatique et ThÃ©orique
//                  UniversitÃ© de MontrÃ©al.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Wed Aug  2 17:43:21 2006
// $Revision: 59 $
// $Id: BasePair.cc 59 2006-11-15 21:25:50Z larosem $
//


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "BasePair.h"

namespace annotate {

BasePair::face_vector BasePair::reverseFaces(const face_vector& aFaces)
{
	BasePair::face_vector faces;
	if(0 < aFaces.size())
	{
		faces.resize(aFaces.size());
		face_vector::const_iterator it;
		for(it = aFaces.begin(); it != aFaces.end(); ++it)
		{
			faces.push_back(std::pair<const mccore::PropertyType*, const mccore::PropertyType*>(it->second, it->first));
		}
	}
	return faces;
}

}; // namespace annotate
