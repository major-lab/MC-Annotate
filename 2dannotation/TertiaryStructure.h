#ifndef _annotate_TertiaryStructure_H_
#define _annotate_TertiaryStructure_H_

#include <mccore/GraphModel.h>

#include "Cycle.h"

namespace annotate
{
	class TertiaryStructure
	{
	public:		
		// LIFECYCLE ------------------------------------------------------------
		TertiaryStructure();
		virtual ~TertiaryStructure();
		
		// ACCESS ---------------------------------------------------------------
		
		const std::string& name() const {return mName;}
		void name(const std::string& aName) {mName = aName;}
		
		const std::list<Cycle>& getCycles() const {return mCycles;}
		const mccore::GraphModel& getModel() const {return mModel;}
		const std::set<mccore::ResId>& getResidues() const {return mResidues;}
		
		// OPERATORS ------------------------------------------------------------
		
		// METHODS --------------------------------------------------------------
		void update(const AnnotateModel& aModel);
		void addCycle(const Cycle& aCycle);
		
		
				
	private:
		std::string mName;
		
		std::set<mccore::ResId> mResidues;
		std::list<Cycle> mCycles;
		mccore::GraphModel mModel;
		
		void clear();
	};
}

#endif /*_annotate_TertiaryStructure_H_*/