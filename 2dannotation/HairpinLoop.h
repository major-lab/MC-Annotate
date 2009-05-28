#ifndef _annotate_HairpinLoop_H_
#define _annotate_HairpinLoop_H_

#include "Stem.h"

namespace annotate
{
	class HairpinLoop
	{
	public:
		HairpinLoop();
		HairpinLoop(const std::vector<const Residue *>& aResidues);
		~HairpinLoop();
		
		const std::vector< const Residue* >& getResidues() const;
		
	private:
		std::vector< const Residue * > mResidues;
	};
}

#endif /*_annotate_HairpinLoop_H_*/
