#ifndef _annotate_SecondaryStructure_H_
#define _annotate_SecondaryStructure_H_

namespace annotate
{
	class SecondaryStructure
	{		
	public:
				
		// LIFECYCLE ------------------------------------------------------------
		SecondaryStructure () {  }
	    virtual ~SecondaryStructure () { }

		// OPERATORS ------------------------------------------------------------

		// ACCESS ---------------------------------------------------------------

		// METHODS --------------------------------------------------------------
		virtual bool isAdjacent(const SecondaryStructure& aStruct) const = 0;
	private:
	};
}

#endif /*_annotate_SecondaryStructure_H_*/
