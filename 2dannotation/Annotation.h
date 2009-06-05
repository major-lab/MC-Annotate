#ifndef _annotation_Annotation_H_
#define _annotation_Annotation_H_

#include <set>
#include <string>

namespace annotate
{
	class AnnotateModel;
	
	class Annotation
	{
	public:
		Annotation() {}
		virtual ~Annotation() {}
		
		// PROVIDED
		const std::set<std::string >& requires() const;
		void addRequirement(const std::string& aRequirement);
		
		// MUST IMPLEMENT
		virtual void update(const AnnotateModel& aModel) = 0;		
		virtual std::string output() const = 0;
		virtual const std::string provides() const = 0;
		
	protected:
		virtual void clear() = 0;
		std::set< std::string > mRequirements;
	};
}

#endif /*_annotation_Annotation_H_*/
