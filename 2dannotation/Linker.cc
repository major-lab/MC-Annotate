#include "Linker.h"

namespace annotate 
{
	Linker::Linker()
	{
	}
		
	Linker::Linker(
		const std::vector<mccore::ResId>& aResidues, 
		const StemConnection& aStart,
		const StemConnection& aEnd)
	{
		mResidues = aResidues;
		mStart = aStart;
		mEnd = aEnd;
	}
	
	Linker::~Linker()
	{
		clear();
	}
		
	const std::vector<mccore::ResId>& Linker::getResidues() const
	{
		return mResidues;
	}
	
	void Linker::clear()
	{
		mResidues.clear();
	}
	
	bool Linker::isEmpty() const
	{
		bool bIsEmpty = true;
		
		if((mStart.isValid() && mEnd.isValid()) || 0 < mResidues.size())
		{
			bIsEmpty = false;
		}
		return bIsEmpty;
	}
	
	bool Linker::operator== (const Linker &other) const
	{
		bool bEqual = false;
		
		if((mStart == other.mStart && mEnd == other.mEnd) // same
			|| (mStart == other.mEnd && mEnd == other.mStart)) // reversed
		{
			// Same order
			bEqual = true;
		}
		return bEqual;
	}
	
	bool Linker::operator!= (const Linker &other) const
	{
		bool bEqual = operator==(other);
		
		return !bEqual;
	}
	
	bool Linker::operator< (const Linker &other) const
    {
    	bool bIsSmaller = false;
    	
    	if(!isEmpty() && !other.isEmpty())
    	{
    		ResId thisId;
    		if(mStart.isValid())
	    	{
	    		thisId = mStart.getResidue();
	    	}
	    	else
	    	{
	    		thisId = mResidues.front();
	    	}
	    	
	    	ResId otherId;
	    	if(other.mStart.isValid())
	    	{
	    		otherId = other.mStart.getResidue();
	    	}
	    	else
	    	{
	    		otherId = other.mResidues.front();
	    	}
	    	bIsSmaller = thisId < otherId;
    	}
    	
      return bIsSmaller;
    }
	
	void Linker::order()
	{
		if(0 < mResidues.size())
		{
			std::sort(mResidues.begin(), mResidues.end());
		}
		
		if(mStart.isValid() && mEnd.isValid())
		{
			mccore::ResId startResId = mStart.getResidue();
			mccore::ResId endResId = mEnd.getResidue();
			if(endResId < startResId)
			{
				std::swap(mStart, mEnd);
			}
		}else if(mStart.isValid() && 0 < mResidues.size())
		{
			mccore::ResId startResId = mStart.getResidue();
			if(startResId < mResidues.back())
			{
				std::swap(mStart, mEnd);
			}			
		}else if(mEnd.isValid() && 0 < mResidues.size())
		{
			mccore::ResId endResId = mEnd.getResidue();
			if(endResId < mResidues.back())
			{
				std::swap(mStart, mEnd);
			}
		}
		else
		{
			// TODO : Add an exception here, the linker is invalid
		}
	}
	
	void Linker::reverse()
	{
		std::reverse(mResidues.begin(), mResidues.end());
		std::swap(mStart, mEnd);
	}
};