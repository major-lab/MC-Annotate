#ifndef _annotate_Linker_H_
#define _annotate_Linker_H_

#include "SecondaryStructure.h"
#include "Stem.h"
#include <vector>

namespace annotate
{	
	class Linker : public SecondaryStructure
	{
	public:	
		Linker();
		
		Linker(
			const std::vector<mccore::ResId>& aResidues, 
			const StemConnection& aStart,
			const StemConnection& aEnd);
		~Linker();
		
		const std::vector<mccore::ResId >& getResidues() const;
			
		void clear();
		
		bool isEmpty() const;
		
		bool operator== (const Linker& other) const;
		bool operator!= (const Linker& other) const;
		bool operator< (const Linker& other) const;
		
		void order();
		void reverse();
		
		const StemConnection& getStart() const {return mStart;}
		const StemConnection& getEnd() const {return mEnd;}
		
	protected:
		std::vector<mccore::ResId> mResidues;
		StemConnection mStart;
		StemConnection mEnd;
	};
}

#endif /*_annotate_Linker_H_*/
