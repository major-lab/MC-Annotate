/*
 * AnnotatedModelSQLExporter.h
 *
 *  Created on: Feb 16, 2010
 *      Author: blanchmf
 */

#ifndef _annotate_AnnotateModelSQLExporter_h_
#define _annotate_AnnotateModelSQLExporter_h_

#include <string>
#include <sstream>

namespace annotate
{
class AnnotateModel;

class AnnotatedModelSQLExporter
{
public:
	AnnotatedModelSQLExporter();

	void addModel(const AnnotateModel& aModel);

	std::string toSQL() const;
private:
	std::ostringstream mOss;

	unsigned int muiModelId;
	unsigned int muiResidueId;

	std::string databaseSchema() const;
	std::string conformationsToSQL(const annotate::AnnotateModel& aModel) const;
};

}; // namespace annotate

#endif /* _annotate_AnnotateModelSQLExporter_h_ */
