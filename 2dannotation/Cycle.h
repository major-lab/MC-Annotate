#ifndef _annotate_Cycle_H_
#define _annotate_Cycle_H_

#include <mccore/GraphModel.h>

#include "BaseInteraction.h"

namespace annotate
{	
	class Cycle
	{
	public:	
		// LIFECYCLE ------------------------------------------------------------
		Cycle(const mccore::GraphModel& aModel);
		virtual ~Cycle();
		
		// ACCESS ---------------------------------------------------------------
		const mccore::GraphModel& getModel() const;
		
		// OPERATORS ------------------------------------------------------------
		
		// METHODS --------------------------------------------------------------
		void order();
		
	private:
		mccore::GraphModel mModel;
		std::list<const BaseInteraction*> mInteractions;
		std::list<const mccore::Residue*> mResidues;
		std::list<unsigned int> mTopology;
		void clear();
	};
}

#endif /*_annotate_Linker_H_*/