/*
 * LabeledResId.cc
 *
 *  Created on: Sep 6, 2009
 *      Author: blanchmf
 */

#include "LabeledResId.h"

namespace annotate
{

	bool LabeledResId::isValid() const
	{
		// Insure that the instance has been initialized
		return *this != LabeledResId();
	}
}
