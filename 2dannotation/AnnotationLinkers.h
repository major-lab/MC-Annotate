#ifndef _annotate_AnnotationLinkers_H_
#define _annotate_AnnotationLinkers_H_

#include <vector>
#include <string>

#include "Annotation.h"
#include "Linker.h"

namespace annotate
{
	class Stem;
	class AnnotationLinkers : public Annotation
	{
	public:
		AnnotationLinkers();
		virtual ~AnnotationLinkers();
		
		virtual void update(AnnotateModel& aModel);		
		virtual std::string output() const;
		
		const std::vector< Linker >& getLinkers() const;
		
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
	private:
		static std::string mstrAnnotationName;
		struct stResidueInfo
		{
			mccore::ResId resId;
			const Stem* pStem;
		};
		std::vector< std::vector<stResidueInfo> > mResidueInfos;
		std::vector< Linker > mLinkers;
		virtual void clear();
		
	  	void allResiduesLinker(
	  		const std::vector<stResidueInfo>& chainInfo, 
			std::set<Linker>& linkers) const;
	  	
		void computeResidueInfos(const AnnotateModel& aModel);
		void updateLinkers(std::set<Linker>& linkers) const;
		void updateChainLinkers(
			const std::vector<stResidueInfo>& chainInfo, 
			std::set<Linker>& linkers) const;
		Linker createLinker(
	  		const std::vector<mccore::ResId>& aResidues, 
	  		const Stem* apStem1, 
	  		const mccore::ResId& aResId1, 
	  		const Stem* apStem2, 
	  		const mccore::ResId& aResId2) const;
			
		std::string outputLinker(const Linker& aLinker) const;
	};
	
}
#endif /*_annotate_AnnotationLinkers_H_*/
