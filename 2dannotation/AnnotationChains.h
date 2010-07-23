#ifndef _annotate_AnnotationChains_H_
#define _annotate_AnnotationChains_H_

#include "Annotation.h"
#include <mccore/ResId.h>
#include <list>
#include <map>

namespace annotate
{
	class AnnotationChains : public Annotation
	{
	public:
		typedef std::list<mccore::ResId> chain_content;
		typedef std::map<unsigned char, chain_content > chain_map;

		AnnotationChains();
		virtual ~AnnotationChains();

		virtual void update(AnnotateModel& aModel);
		virtual std::string output() const;

		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}

		// METHODS -------------------------------------------------------------
		const std::list<mccore::ResId>& getChain(unsigned char aucChainId) const;

		// ACCESSOR ------------------------------------------------------------
		const chain_map& chains() const {return mChains;}

	private:
		typedef std::pair<unsigned char, chain_content > chain_entry;

		static std::string mstrAnnotationName;
		chain_map mChains;

		virtual void clear();
	};
}
#endif /*_annotate_AnnotationChains_H_*/
