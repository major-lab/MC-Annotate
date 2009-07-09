#ifndef _annotate_SecondaryStructure_H_
#define _annotate_SecondaryStructure_H_

#include <string>

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
	protected:
		std::string mName;
	
	private:
	};
}

#endif /*_annotate_SecondaryStructure_H_*/
