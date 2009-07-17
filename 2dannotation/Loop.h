#ifndef _annotate_Loop_H_
#define _annotate_Loop_H_

#include "SecondaryStructure.h"
#include "Linker.h"
#include <vector>

namespace annotate
{	
	class Loop : public SecondaryStructure
	{
	public:		
		Loop();
		Loop(const std::vector<Linker>& aLinkers);
		virtual ~Loop();
		
		const std::vector<Linker>& getLinkers() const;
		void reverse();
		void append(const Loop& aLoop);
				
		void clear();
		
		bool contains(const mccore::ResId& aResId) const;
		
		bool operator ==(const Loop& other) const;
		bool isAdjacent(const SecondaryStructure& aStruct) const;
		
		std::string describe() const;
		
	protected:
		std::vector<Linker> mLinkers;
	};
}

#endif /*_annotate_Loop_H_*/
