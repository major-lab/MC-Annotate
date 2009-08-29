#ifndef _annotate_SecondaryStructure_H_
#define _annotate_SecondaryStructure_H_

#include <string>

#include "BaseInteraction.h"

namespace annotate
{
	class SecondaryStructure
	{
	public:

		// LIFECYCLE ------------------------------------------------------------
		SecondaryStructure (const std::string& aName = "") { mName = aName;  }
	    virtual ~SecondaryStructure () { }

		// OPERATORS ------------------------------------------------------------

		// ACCESS ---------------------------------------------------------------
		const std::string& name() const {return mName;}
		void name(std::string& aName) {mName = aName;}

		// METHODS --------------------------------------------------------------
		virtual bool isAdjacent(const SecondaryStructure& aStruct) const = 0;

		// TODO : Verify if this couldn't be changed to an operator ==
		/**
		 * @brief Check if given parameter and current structure are the same
		 */
		virtual bool isSame(const SecondaryStructure& aStruct) const = 0;

		/**
		 * @brief getBaseInteractions forming the secondary structure.
		 * @details BaseInteractions are unspecialized interactions between
		 * residues.  This effectively return only one interaction between any
		 * pair of interacting residues.  The interactions are unqualified, ( no
		 * pairing, stacking, adjacency, etc... ).
		 */
		virtual std::set<BaseInteraction> getBaseInteractions() const = 0;
	protected:
		std::string mName;

	private:
	};
}

#endif /*_annotate_SecondaryStructure_H_*/
