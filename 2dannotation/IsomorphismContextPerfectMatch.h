#ifndef _graphtest_IsomorphismContextPerfectMatch_H_
#define _graphtest_IsomorphismContextPerfectMatch_H_

#include "GraphIsoUllman.h"

#include "mccore/GraphModel.h"

namespace graphtest
{
	
	// bool potentialCheck_PerfectMatch(unsigned int i, unsigned int j, void* apContext);

	class PerfectMatchContext : public graphtest::CGraphIsomorphism::IsoContext
	{
	public:
		PerfectMatchContext(
			const mccore::GraphModel& aG1, 
			const mccore::GraphModel& aG2);
		~PerfectMatchContext();
		
		// ACCESS ---------------------------------------------------------------
		const std::vector<std::vector<unsigned int> >& getMatches() {return mMatches;}
		/*
		const mccore::GraphModel& getGraph1() const {return mG1;}
		const mccore::GraphModel& getGraph2() const {return mG2;}
		const std::vector<mccore::GraphModel::label> getMapping1() const {return mMapping1;}
		const std::vector<mccore::GraphModel::label> getMapping2() const {return mMapping2;}
		*/
					
		// METHODS --------------------------------------------------------------
		virtual bool potentialCheck(unsigned int i, unsigned int j);
		virtual bool processMatch(std::vector<unsigned int>& match);
		virtual bool isomorphismCheck(const std::vector<unsigned int>& match);
	private:
		mccore::GraphModel mG1;
		mccore::GraphModel mG2;
		std::vector<std::vector<unsigned int> > mMatches;
		
		std::vector<mccore::GraphModel::label> mMapping1;
		std::vector<mccore::GraphModel::label> mMapping2;
		
		std::map<mccore::GraphModel::label, unsigned int> mReverseMap1;
	};	
}

#endif // _graphtest_IsomorphismContextPerfectMatch_H_
