#ifndef _annotate_StemConnection_H_
#define _annotate_StemConnection_H_

#include "BasePair.h"
#include "Stem.h"

#include <vector>

namespace annotate
{
	class StemConnection
	{
	public:
		// LIFECYCLE ------------------------------------------------------------
		StemConnection();
		StemConnection(const Stem& aStem, const Stem::enConnection& aeConnect);

		// ACCESS ---------------------------------------------------------------
		const SecondaryStructure* getStructure() const {return mpStructure;}
		const Stem::enConnection getConnection() const { return meConnection; }

		// METHODS --------------------------------------------------------------
		mccore::ResId nextId() const throw(mccore::NoSuchElementException);
		BasePair getPair() const throw(mccore::NoSuchElementException);
		mccore::ResId getResidue() const;
		bool isValid() const;

		/**
		 * @brief Verify the two connections refers to a both residues of the
		 * same end of a stem.
    	 * @return true if they form the pair at the end of a stem.
    	 */
		// bool connects(const StemConnection& aConnection) const;

		// OPERATOR -------------------------------------------------------------
		bool operator== (const StemConnection &other) const;
	private:
		const SecondaryStructure* mpStructure;
		Stem::enConnection meConnection;
	};
}

#endif /*_annotate_StemConnection_H_*/
