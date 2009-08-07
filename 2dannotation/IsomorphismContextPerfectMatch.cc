#include "IsomorphismContextPerfectMatch.h"

#include "mccore/Messagestream.h"
#include "mccore/PropertyType.h"

namespace graphtest
{
	PerfectMatchContext::PerfectMatchContext(
		const mccore::GraphModel& aG1,
		const mccore::GraphModel& aG2)
	{
		mG1 = aG1;
		mG2 = aG2;
		
		unsigned int i;
		mccore::GraphModel::const_iterator it = aG1.begin();
		for(it = aG1.begin(), i = 0; it != aG1.end(); ++ it, ++ i)
		{
			mccore::GraphModel::label label = aG1.getVertexLabel (const_cast< mccore::Residue* > (&(*it)));
			mMapping1.push_back(label);
			std::pair<mccore::GraphModel::label, unsigned int> reverse(label, i);
			mReverseMap1.insert(reverse);
		}
		for(it = aG2.begin(); it != aG2.end(); ++ it)
		{
			mccore::GraphModel::label label = aG2.getVertexLabel (const_cast< mccore::Residue* > (&(*it)));
			mMapping2.push_back(label);
		}
	}
	PerfectMatchContext::~PerfectMatchContext()
	{
		mMatches.clear();
		mReverseMap1.clear();
		mMapping1.clear();
		mMapping2.clear();
	}

	bool PerfectMatchContext::potentialCheck(unsigned int i, unsigned int j)
	{
		const mccore::Residue* pRes1 = mG1.internalGetVertex(mMapping1[i]);
		const mccore::Residue* pRes2 = mG2.internalGetVertex(mMapping2[j]);
		bool bPotential = (pRes1->getType() == pRes2->getType());
		if(bPotential)
		{
			unsigned int uiNbStacks1 = 0;
			unsigned int uiNbPairs1 = 0;
			unsigned int uiNbLinks1 = 0;
			unsigned int uiNbStacks2 = 0;
			unsigned int uiNbPairs2 = 0;
			unsigned int uiNbLinks2 = 0;
			
			std::list<mccore::GraphModel::label>::const_iterator it;
			std::list<mccore::GraphModel::label> neighbors1;
			neighbors1 = mG1.internalNeighborhood(mMapping1[i]);
			for(it = neighbors1.begin(); it != neighbors1.end(); ++ it)
			{
				try
				{
					mccore::Relation* pRelation = mG1.internalGetEdge(mMapping1[i], *it);
					if (pRelation->isPairing ())
					{
						uiNbPairs1 ++;
					}
					if (pRelation->is (mccore::PropertyType::pAdjacent5p))
					{
						uiNbLinks1 ++;
					}
					if (pRelation->isStacking ())
					{
						uiNbStacks1 ++;
					}
				}catch(mccore::NoSuchElementException& e)
				{
					// No relation, no count, no problem
				}
			}
			std::list<mccore::GraphModel::label> neighbors2;
			neighbors2 = mG2.internalNeighborhood(mMapping2[j]);
			for(it = neighbors2.begin(); it != neighbors2.end(); ++ it)
			{
				try
				{
					mccore::Relation* pRelation = mG2.internalGetEdge(mMapping2[j], *it);
					if (pRelation->isPairing ())
					{
						uiNbPairs2 ++;
					}
					if (pRelation->is (mccore::PropertyType::pAdjacent5p))
					{
						uiNbLinks2 ++;
					}
					if (pRelation->isStacking ())
					{
						uiNbStacks2 ++;
					}
				}catch(mccore::NoSuchElementException& e)
				{
					// No relation, no count, no problem
				}
			}
			bPotential = ((uiNbPairs1 == uiNbPairs2) 
				&& (uiNbStacks1 == uiNbStacks2) 
				&& (uiNbLinks1 == uiNbLinks2));			
		}
		return bPotential;
	}
	
	bool PerfectMatchContext::processMatch(std::vector<unsigned int>& match)
	{
		mMatches.push_back(match);
		return false; // We only care for one match for our test
	}
	
	bool PerfectMatchContext::isomorphismCheck(const std::vector<unsigned int>& match)
	{
		bool isomorphism = true;
		
		for(unsigned int i = 0; i < match.size() && (isomorphism); ++i)
		{
			std::list<mccore::GraphModel::label>::const_iterator it;
			std::list<mccore::GraphModel::label> neighbors1 = mG1.internalNeighborhood(mMapping1[i]);
			for(it = neighbors1.begin(); it != neighbors1.end() && (isomorphism); ++ it)
			{
				mccore::Relation* pRelation = NULL;
				mccore::Relation* pRelation2 = NULL;
				try
				{
					pRelation = mG1.internalGetEdge(mMapping1[i], *it);
					unsigned int j = mReverseMap1[*it];
					unsigned int matchI = match[i];
					unsigned int matchJ = match[j];
					unsigned int mappingI = mMapping2[matchI];
					unsigned int mappingJ = mMapping2[matchJ];
					pRelation2 = mG2.internalGetEdge(mappingI,mappingJ );
				}
				catch(mccore::NoSuchElementException& e)
				{
					isomorphism = false;
				}
				if (isomorphism && pRelation->isPairing ())
				{
					isomorphism = pRelation2->isPairing();
				}
				if (isomorphism && (pRelation->is (mccore::PropertyType::pAdjacent5p)))
				{
					isomorphism = (pRelation2->is (mccore::PropertyType::pAdjacent5p));
				}
				if (isomorphism && pRelation->isStacking ())
				{
					isomorphism = pRelation2->isStacking ();
				}
			}
		}
		return isomorphism;
	}
}