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
		
		virtual void update(const AnnotateModel& aModel);		
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
		std::vector<stResidueInfo> mResidueInfos;
		std::vector< Linker > mLinkers;
		virtual void clear();
		
		Linker findLinker(const StemConnection& aConnection) const;
		void findLinker(
	  		const Stem* apStem, 
	  		const Stem::enConnection& aeConnect, 
	  		std::set<Linker>& outLinkerSet) const;
	  	Linker allResiduesLinker() const;
	  	
		void computeResidueInfos(const AnnotateModel& aModel);
		std::vector<stResidueInfo>::const_iterator findResidueInfo(
			const mccore::ResId& aResId) const;
	};
	
}
#endif /*_annotate_AnnotationLinkers_H_*/
