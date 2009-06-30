#ifndef _annotate_AlgorithmExtra_h_
#define _annotate_AlgorithmExtra_h_

#include <set>
#include <algorithm>

namespace annotate
{
	//----------------------------------------------------------------------
	template < class T >
	std::set< T > SetIntersection( std::set< T > & set1, std::set< T > & set2 )
	{
		std::set< T > setInter;
		std::insert_iterator< std::set< T > > iter( setInter, setInter.begin() );
	
		std::set_intersection( set1.begin(), set1.end(),
	          set2.begin(), set2.end(),
	          iter );
	
		return setInter;
	};
	
	template < class T >
	std::set< T > SetIntersection( const std::set< T > & set1, const std::set< T > & set2 )
	{
		std::set< T > setInter;
		std::insert_iterator< std::set< T > > iter( setInter, setInter.begin() );
	
		std::set_intersection( set1.begin(), set1.end(),
	          set2.begin(), set2.end(),
	          iter );
	
		return setInter;
	};
		
	//----------------------------------------------------------------------
	template < class T >
	std::set< T > SetDifference( std::set< T > & set1, std::set< T > & set2 )
	{
		std::set< T > setDiff;
		std::insert_iterator< std::set< T > > iter( setDiff, setDiff.begin() );
	
		std::set_difference( set1.begin(), set1.end(),
	        set2.begin(), set2.end(),
	        iter );
	
		return setDiff;
	};
	
	template <class InputIterator1, class InputIterator2>
  	bool set_intersects ( InputIterator1 first1, InputIterator1 last1,
                          InputIterator2 first2, InputIterator2 last2)
	{
		bool bIntersects = false;
		while (first1!=last1 && first2!=last2)
		{
		    if ((*first1) < (*first2))
		    {
		    	 ++first1;
		    }
		    else if ((*first2) < (*first1))
		    {
		    	 ++first2;
		    }
		    else 
		    { 
		    	bIntersects = true;
		    	break;
		    }
		}
		return bIntersects;
	}
	
	template <class InputIterator1, class InputIterator2, class Compare>
  	bool set_intersects ( InputIterator1 first1, InputIterator1 last1,
                          InputIterator2 first2, InputIterator2 last2,
                          Compare comp)
	{
		bool bIntersects = false;
		while (first1!=last1 && first2!=last2)
		{
		    if (comp(*first1, *first2))
		    {
		    	 ++first1;
		    }
		    else if (comp(*first2, *first1))
		    {
		    	 ++first2;
		    }
		    else 
		    { 
		    	bIntersects = true;
		    	break;
		    }
		}
		return bIntersects;
	}
}

#endif /*_annotate_AlgorithmExtra_h_*/
