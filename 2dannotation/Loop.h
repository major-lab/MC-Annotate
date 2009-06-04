#ifndef _annotate_Loop_H_
#define _annotate_Loop_H_

#include "Linker.h"
#include <vector>

namespace annotate
{	
	class Loop
	{
	public:		
		Loop();
		Loop(const std::vector<Linker>& aLinkers);
		~Loop();
		
		const std::vector<Linker>& getLinkers() const;
				
		void clear();
		
	protected:
		std::vector<Linker> mLinkers;
	};
}

#endif /*_annotate_Loop_H_*/
