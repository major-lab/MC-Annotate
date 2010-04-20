/*
 * AnnotatedModelSQLExporter.cc
 *
 *  Created on: Feb 16, 2010
 *      Author: blanchmf
 */

#include "AnnotatedModelSQLExporter.h"

#include "AnnotateModel.h"

#include "mccore/Pdbstream.h"

namespace annotate
{

AnnotatedModelSQLExporter::AnnotatedModelSQLExporter()
{
	muiModelId = 0;
	mOss << databaseSchema();
}

std::string AnnotatedModelSQLExporter::toSQL() const
{
	return mOss.str();
}

void AnnotatedModelSQLExporter::addModel(const AnnotateModel& aModel)
{
	std::ostringstream oss;

	muiModelId ++;

	// Model identification
	oss << "INSERT INTO models VALUES (" << muiModelId;
	oss << ", \"" << aModel.name();
	oss << "\", " << aModel.id() << ");" << std::endl;

	// Output the conformations
	oss << conformationsToSQL(aModel);

	mOss << oss.str();
}

std::string AnnotatedModelSQLExporter::conformationsToSQL(
	const annotate::AnnotateModel& aModel) const
{
	std::ostringstream oss;
	annotate::AnnotateModel::const_iterator i;

	for (i = aModel.begin (); i != aModel.end (); ++i)
	{
		oss << "INSERT INTO residues VALUES (" << muiModelId << ", ";
		oss << "\"" << i->getResId().getChainId() << "\", ";
		oss << "\"" << i->getResId().getResNo() << "\", ";
		oss << "\"" << i->getResId().getInsertionCode() << "\", ";
		oss << "\"" << mccore::Pdbstream::stringifyResidueType (i->getType ()) << "\", ";
		if (i->getType ()->isNucleicAcid ())
		{
			oss << "\"" << i->getPucker() << "\", ";
			oss << "\"" << i->getGlycosyl() << "\");";
		}
		else
		{
			oss << "NULL, NULL);";
		}
		oss << std::endl;
	}
	return oss.str();
}

std::string AnnotatedModelSQLExporter::databaseSchema() const
{
	std::ostringstream oss;
	// Database
	oss << "CREATE DATABASE MCAnnotation;" << std::endl;
	oss << "USE MCAnnotation;" << std::endl;

	// Models
	oss << "CREATE TABLE models (" << std::endl;
	oss << "\tid INTEGER UNSIGNED NOT NULL AUTO_INCREMENT, " << std::endl;
	oss << "\tpdbfile VARCHAR(4), " << std::endl;
	oss << "\tmodel INTEGER, " << std::endl;
	oss << "\tPRIMARY KEY(id)" << std::endl;
	oss << ");" << std::endl;

	// Residues
	oss << "CREATE TABLE residues (" << std::endl;
	oss << "\tid INTEGER UNSIGNED NOT NULL AUTO_INCREMENT, " << std::endl;
	oss << "\tmodelId INTEGER UNSIGNED NOT NULL REFERENCES models(id)," << std::endl;
	oss << "\tchain VARCHAR(20) NOT NULL, " << std::endl;
	oss << "\tresNo INTEGER NOT NULL, " << std::endl;
	oss << "\tinsertionCode SMALL INT, " << std::endl;
	oss << "\ttype VARCHAR(4), " << std::endl;
	oss << "\tpucker VARCHAR(8), " << std::endl;
	oss << "\tglycosyl VARCHAR(4), " << std::endl;
	oss << "\tPRIMARY KEY(id)" << std::endl;
	oss << ");" << std::endl;

	return oss.str();
}

};
