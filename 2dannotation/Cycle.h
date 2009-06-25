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
		const std::list<unsigned int>& getTopology() const;
		
		// OPERATORS ------------------------------------------------------------
		
		// METHODS --------------------------------------------------------------
		void order();
		
	private:
	
		class Topology
		{
		public:
			typedef std::list<const BaseInteraction*>::const_iterator const_interaction_iterator;
			Topology() {}
			
			Topology(const const_interaction_iterator &it, bool abReverse) 
				: mItStart(it), mbReverse(abReverse) 
			{}
			~Topology() {mStrands.clear();}
			
			// ACCESS ---------------------------------------------------------------
			bool empty() const {return mStrands.empty();}
			const_interaction_iterator getFirstIterator() const {return mItStart;}
			const std::list<unsigned int>& getStrands() const {return mStrands;}
		
			// OPERATORS ------------------------------------------------------------
			bool operator <(const Topology& aTopo) const;
			
			// METHODS --------------------------------------------------------------
			void addStrandLength(unsigned int auiLength) 
			{
				mStrands.push_back(auiLength);
			}
			
			
		private:
			const_interaction_iterator mItStart;
			bool mbReverse;
			std::list<unsigned int> mStrands;
		};
		
		typedef std::list<const BaseInteraction*> interactions_list;
		typedef interactions_list::const_iterator interactions_list_iterator;
		typedef std::list<unsigned int> strand_list;
		mccore::GraphModel mModel;
		AnnotationInteractions mInteractionsAnnotation;
		std::list<const BaseInteraction*> mInteractions;
		std::list<const mccore::Residue*> mResidues;
		Topology mTopology;
		void clear();
		
		Topology computeTopology() const;
		void completeTopology(Topology& aCandidate) const;
		void completeReverseTopology(Topology& aCandidate) const;
		unsigned int advanceTopoIterator(interactions_list_iterator& it) const;
		unsigned int rewindTopoIterator(interactions_list_iterator& it) const;
	};
	
	
}

#endif /*_annotate_Linker_H_*/