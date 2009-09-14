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
		typedef std::pair<Linker, const Stem*> linker_info;
		AnnotationLinkers();
		virtual ~AnnotationLinkers();

		// METHODS -------------------------------------------------------------
		virtual void update(AnnotateModel& aModel);
		virtual std::string output() const;

		// ACCESS --------------------------------------------------------------
		const std::vector< linker_info >& getLinkers() const { return mLinkers;	}
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
		std::vector< linker_info > mLinkers;
		std::set<mccore::ResId> mUnconnected;
		virtual void clear();

	  	void allResiduesLinker(const std::vector<stResidueInfo>& chainInfo);

		void computeResidueInfos(const AnnotateModel& aModel);
		void removeUnconnectedLinkers();
		bool isConnected(const Linker& aLinker) const;
		void updateLinkers();
		void updateChainLinkers(const std::vector<stResidueInfo>& chainInfo);
		Linker createLinker(
			const std::vector<LabeledResId>& aResidues,
			const SecondaryStructure* apStartStruct,
			const SecondaryStructure* apEndStruct) const;

		std::list<std::list<unsigned int> > getLinkerRanges(
				const std::vector<stResidueInfo>& chainInfo) const;

		std::list<AnnotationLinkers::linker_info>  getLinker(
			const std::vector<stResidueInfo>& aChainInfo,
			const std::list<unsigned int>& aRange);

		bool shouldConnect(
			const linker_info& aFirst,
			const linker_info& aSecond) const;
		void connectPseudoLinkers();
		void connectLinkersToLinkers();

		std::string outputLinker(const Linker& aLinker) const;
	};

}
#endif /*_annotate_AnnotationLinkers_H_*/
