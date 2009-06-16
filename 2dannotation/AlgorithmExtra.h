#ifndef _annotate_AlgorithmExtra_h_
#define _annotate_AlgorithmExtra_h_

namespace annotate
{
	//----------------------------------------------------------------------
	template < class T >
	set< T > SetIntersection( set< T > & set1, set< T > & set2 )
	{
	 set< T > setInter;
	 insert_iterator< set< T > > iter( setInter, setInter.begin() );
	
	 set_intersection( set1.begin(), set1.end(),
	          set2.begin(), set2.end(),
	          iter );
	
	 return setInter;
	};
	
	//----------------------------------------------------------------------
	template < class T >
	set< T > SetDifference( set< T > & set1, set< T > & set2 )
	{
	 set< T > setDiff;
	 insert_iterator< set< T > > iter( setDiff, setDiff.begin() );
	
	 set_difference( set1.begin(), set1.end(),
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
		    if (*first1<*first2)
		    {
		    	 ++first1;
		    }
		    else if (*first2<*first1)
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
