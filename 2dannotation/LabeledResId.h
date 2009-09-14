/*
 * LabeledResId.h
 *
 *  Created on: Sep 6, 2009
 *      Author: blanchmf
 */

#ifndef _annotate_labeledResId_H_
#define _annotate_labeledResId_H_

#include "mccore/GraphModel.h"
#include "mccore/ResId.h"

namespace annotate
{
	class LabeledResId : public mccore::ResId
	{
	public:
		/**
		 * Initializes the object.
		 */
		LabeledResId (int n = -1) : mccore::ResId(n) {}

		/**
		 * Initializes the object.
		 */
		LabeledResId (char c, int n, char i = ' ') : mccore::ResId(c, n, i) {}

		/**
		 * Initializes the object with the other's content.
		 * @param other the object to copy.
		 */
		LabeledResId (const ResId &other) : mccore::ResId(other) {}

		/**
		 * Initializes the object with the other's content and a label.
		 * @param other the object to copy.
		 * @param label the label to apply
		 */
		LabeledResId (const ResId &other, const mccore::GraphModel::label &aLabel)
		: mccore::ResId(other)
		, mLabel(aLabel)
		{}

		/**
		 * Initializes the structure with a text representation.
		 * @param str the text representation.
		 */
		LabeledResId (const char *str) : mccore::ResId(str) {}

		// ACCESS --------------------------------------------------------------
		const mccore::GraphModel::label label() const {return mLabel;}

		// METHODS -------------------------------------------------------------

		/**
		 * Verify that the resid has been initialized with values different than
		 * the default.
		 */
		bool isValid() const;

	private:
		mccore::GraphModel::label mLabel;
	};
}

#endif /* _annotate_labeledResId_H_ */
