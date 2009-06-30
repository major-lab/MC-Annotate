#ifndef _annotate_Cycle_H_
#define _annotate_Cycle_H_

#include <mccore/GraphModel.h>

#include "BaseInteraction.h"

#include "AnnotationInteractions.h"

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
		std::list<std::string> getInteractionLabels() const;
		
	private:
		
		typedef std::list<const BaseInteraction*> interactions_list;
		typedef interactions_list::const_iterator interactions_list_iterator;
		typedef std::list<unsigned int> strand_list;
		
		mccore::GraphModel mModel;
		AnnotationInteractions mInteractionsAnnotation;
		std::list<const BaseInteraction*> mInteractions;
		std::list<const mccore::Residue*> mResidues;
		
		void clear();
	};
	
	
}

#endif /*_annotate_Linker_H_*/