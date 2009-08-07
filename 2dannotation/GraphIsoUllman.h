#ifndef _graphtest_GraphIsoUllman_H_
#define _graphtest_GraphIsoUllman_H_

#include "Matrix.h"

namespace graphtest
{	
	//==================================================================
	class CGraphIsomorphism
	{
	public:
		class IsoContext
		{
		public:
			virtual bool potentialCheck(unsigned int i, unsigned int j) = 0;
			virtual bool processMatch(std::vector<unsigned int>& match) = 0;
			virtual bool isomorphismCheck(const std::vector<unsigned int>& match) = 0;
		};
		
		typedef bool (*pfnPotentialCheck)(unsigned int i, unsigned int j, void* apContext);
		
		//----- life cycle -----
		CGraphIsomorphism();
			
		void update(
			unsigned int auiM1NbNodes,
			unsigned int auiM2NbNodes,
			IsoContext* context);
	private:
	
		//----- methods -----
		bool checkIsomorphism(unsigned int I);
		bool alreadyAssigned(unsigned int I, unsigned int j) const;
	
	
		//----- attributes -----
		std::vector<unsigned int> m_M;
		unsigned int m_iAlpha;
		unsigned int m_iBeta;
		IsoContext *m_pContext;
		Matrix< short > m_M0;
	};
}	// namespace graphtest

#endif /*_graphtest_GraphIsoUllman_H_*/
