#ifndef _annotate_AlgorithmExtra_h_
#define _annotate_AlgorithmExtra_h_

#include <algorithm>
#include <set>
#include <memory>

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
	std::set< T > SetUnion( const std::set< T > & set1, const std::set< T > & set2 )
	{
		std::set< T > setUnion;
		std::insert_iterator< std::set< T > > iter( setUnion, setUnion.begin() );

		std::set_union( set1.begin(), set1.end(),
	          set2.begin(), set2.end(),
	          iter );

		return setUnion;
	};

	//----------------------------------------------------------------------
	template < class T >
	std::set< T > SetSymmetricDifference( const std::set< T > & set1, const std::set< T > & set2 )
	{
		std::set< T > setSymDiff;
		std::insert_iterator< std::set< T > > iter( setSymDiff, setSymDiff.begin() );

		std::set_symmetric_difference( set1.begin(), set1.end(),
	          set2.begin(), set2.end(),
	          iter );

		return setSymDiff;
	};


	//----------------------------------------------------------------------
	template < class T >
	std::set< T > SetDifference( const std::set< T > & set1, const std::set< T > & set2 )
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

	//----------------------------------------------------------------------
	template <class T> struct less_ptr : std::binary_function <T,T,bool>
	{
  		bool operator() (const T* x, const T* y) const
    	{
    		return *x < *y;
    	}
	};

	template <class T> struct less_auto_ptr : std::binary_function <T,T,bool>
	{
  		bool operator() (const std::auto_ptr<T>& x, const std::auto_ptr<T>& y) const
    	{
    		return (*x) < (*y);
    	}
	};

	template <class T1, class T2, class T3> class tuple3
	{
	public:
		typedef T1 first_type;
		typedef T2 second_type;
		typedef T3 third_type;

		T1 first;
		T2 second;
		T3 third;

		tuple3() : first(T1()), second(T2()), third(T3()) {}
		tuple3(const T1& x, const T2& y, const T3& z) : first(x), second(y), third(z) {}
		template <class U, class V, class W>
			tuple3(const tuple3<U,V,W> &t)
			: 	first(t.first),
				second(t.second),
				third(t.third) {}
		bool operator ==(const tuple3& t) const
		{
			return ((first == t.first) && (second == t.second) && (third == t.third));
		}
		bool operator !=(const tuple3& t) const
		{
			return !((first == t.first) && (second == t.second) && (third == t.third));
		}

		bool operator <(const tuple3& t) const
		{
			return ((first < t.first)
					|| ((first == t.first)
						&& ((second < t.second) || ((second == t.second) && third < t.third))
						)
					);
		}
		bool operator >(const tuple3& t) const
		{
			return ((t.first < first)
					|| ((t.first == first)
						&& ((t.second < second) || ((t.second == second) && t.third < third))
						)
					);
		}
		bool operator <=(const tuple3& t) const
		{
			return ((first < t.first)
					|| ((first == t.first)
						&& ((second < t.second)
								|| (((second == t.second) && third < t.third))
									|| (third == t.third))
						)
					);
		}
		bool operator >=(const tuple3& t) const
		{
			return ((t.first < first)
					|| ((t.first == first)
						&& ((t.second < second)
								|| (((t.second == second) && t.third < third))
									|| (third == t.third))
						)
					);
		}
	};
}

#endif /*_annotate_AlgorithmExtra_h_*/
