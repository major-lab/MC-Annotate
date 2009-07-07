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
		
		// MUST IMPLEMENT
		virtual void update(const AnnotateModel& aModel) = 0;		
		virtual std::string output() const = 0;
		virtual const std::string& annotationName() = 0;
				
	protected:
		virtual void clear() = 0;
		std::set< std::string > mRequirements;
		
		template <class T> void addRequirement()
		{
			mRequirements.insert(T::AnnotationName());
		}
	};
}

#endif /*_annotation_Annotation_H_*/
