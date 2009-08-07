#include "mccore/Messagestream.h"

#include "GraphIsoUllman.h"

namespace graphtest
{
	//------------------------------------------------------------
	CGraphIsomorphism::CGraphIsomorphism()
	{
	}
	
	//------------------------------------------------------------
	void CGraphIsomorphism::update(
		unsigned int auiM1NbNodes,
		unsigned int auiM2NbNodes,
		IsoContext* apContext)
	{
		m_M0.clear();
				
		// save parameter
		m_pContext = apContext;
	
		// get dims of matrices
		m_iAlpha = auiM1NbNodes;
		m_iBeta  = auiM2NbNodes;
		
		if(m_iAlpha <= m_iBeta)
		{
			// set M; the vector of permutations
			//----------------------------------------
			m_M.clear();
			m_M.resize(m_iAlpha, 0);	
		
			// set and fill M0
			//----------------------------------------
			unsigned int i, j;
			m_M0.resize(m_iAlpha, m_iBeta);
			for ( i = 0; i < m_iAlpha; ++i )
			{
				for ( j = 0; j < m_iBeta; ++j )
				{
					if(m_pContext->potentialCheck(i, j))
					{
						m_M0.set(i, j, 1);
					}
					else
					{
						m_M0.set(i, j, 1);
					}
				}
			}
	
			// optimize M0 matrix
			
			// do the work
			// ---------------------------------------
			checkIsomorphism( 0 );
		
		
			// delete temp matrices
			m_M.clear();
		}
	}

	//------------------------------------------------------------
	bool CGraphIsomorphism::alreadyAssigned(unsigned int I, unsigned int j) const
	{
		bool bAssigned = false;
		for( unsigned int i = 0; i < I && !bAssigned; ++ i)
		{
			bAssigned = (j == m_M[i]);
		}
		return bAssigned;
	}
	
	//------------------------------------------------------------
	bool CGraphIsomorphism::checkIsomorphism( unsigned int I )
	{
		bool bContinue = true;
		if(I == m_iAlpha)
		{
			if(	m_pContext->isomorphismCheck(m_M))
			{
				bContinue = m_pContext->processMatch(m_M);
			}
		}
		else
		{
			for (unsigned int j = 0; j < m_iBeta && bContinue; ++j )
			{
				if ( m_M0.get( I , j ) != 0 )
				{
					if(!alreadyAssigned(I, j))
					{
						m_M[I] = j;
						bContinue = checkIsomorphism( I+1 );
					}
				}
			}
		}
		return bContinue;
	}
	
} // namespace graphtest