/*
 * GetModelRangeFunctor.h
 *
 *  Created on: Jan 5, 2010
 *      Author: blanchmf
 */

#ifndef GETMODELRANGEFUNCTOR_H_
#define GETMODELRANGEFUNCTOR_H_

#include "ModelInfo.h"

template <class T>
class GetModelRangeFunctor
{
public:
	typedef typename std::set< T >::const_iterator const_iterator;
	typedef std::pair<const_iterator, const_iterator> const_range;

	const_range operator ()(const std::set<T>& aSet, const annotate::ModelInfo& aModelInfo) const
	{
		const_range modelRange(aSet.end(), aSet.end());
		const_iterator it = aSet.begin();
		while(it != aSet.end() && it->getModelInfo() < aModelInfo)
		{
			++ it;
		}
		modelRange.first = it;
		while(it != aSet.end() && !(aModelInfo < it->getModelInfo()))
		{
			++ it;
		}
		modelRange.second = it;
		return modelRange;
	}
};

#endif /* GETMODELRANGEFUNCTOR_H_ */
