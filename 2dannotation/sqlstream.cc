/*
 * sqlstream.cc
 *
 *  Created on: Mar 8, 2010
 *      Author: blanchmf
 */

#include "sqlstream.h"

#include "AnnotateModel.h"
#include "AnnotationInteractions.h"
#include <mccore/Pdbstream.h>
#include <mccore/Residue.h>

namespace annotate {

oSQLStream::oSQLStream()
{
	*(std::ostream*)this << "CREATE DATABASE MCAnnotation;" << std::endl;
	*(std::ostream*)this << "USE MCAnnotation;" << std::endl;

	// Models
	*(std::ostream*)this << "CREATE TABLE models (" << std::endl;
	*(std::ostream*)this << "\tid INTEGER UNSIGNED NOT NULL, " << std::endl;
	*(std::ostream*)this << "\tpdbfile VARCHAR(4), " << std::endl;
	*(std::ostream*)this << "\tmodel INTEGER, " << std::endl;
	*(std::ostream*)this << "\tPRIMARY KEY(id)" << std::endl;
	*(std::ostream*)this << ");" << std::endl;

	// Residues
	*(std::ostream*)this << "CREATE TABLE residues (" << std::endl;
	*(std::ostream*)this << "\tid INTEGER UNSIGNED NOT NULL, " << std::endl;
	*(std::ostream*)this << "\tmodelId INTEGER UNSIGNED NOT NULL REFERENCES models(id)," << std::endl;
	*(std::ostream*)this << "\tchain VARCHAR(20) NOT NULL, " << std::endl;
	*(std::ostream*)this << "\tresNo INTEGER NOT NULL, " << std::endl;
	*(std::ostream*)this << "\tinsertionCode SMALL INT, " << std::endl;
	*(std::ostream*)this << "\ttype VARCHAR(4), " << std::endl;
	*(std::ostream*)this << "\tpucker VARCHAR(8), " << std::endl;
	*(std::ostream*)this << "\tglycosyl VARCHAR(4), " << std::endl;
	*(std::ostream*)this << "\tPRIMARY KEY(id)" << std::endl;
	*(std::ostream*)this << ");" << std::endl;

	// Interactions
	*(std::ostream*)this << "CREATE TABLE interactions (" << std::endl;
	*(std::ostream*)this << "\tid INTEGER UNSIGNED NOT NULL, " << std::endl;
	*(std::ostream*)this << "\tNuc1 INTEGER UNSIGNED NOT NULL REFERENCES residues(id), " << std::endl;
	*(std::ostream*)this << "\tNuc2 INTEGER UNSIGNED NOT NULL REFERENCES residues(id), " << std::endl;
	*(std::ostream*)this << "\ttype ENUM('stack', 'pairing'), " << std::endl;
	*(std::ostream*)this << "\tPRIMARY KEY(id)" << std::endl;
	*(std::ostream*)this << ");" << std::endl;

	// Stacks
	*(std::ostream*)this << "CREATE TABLE stacks (" << std::endl;
	*(std::ostream*)this << "\tid INTEGER UNSIGNED NOT NULL, " << std::endl;
	*(std::ostream*)this << "\tinteraction INTEGER UNSIGNED NOT NULL REFERENCES interactions(id), " << std::endl;
	*(std::ostream*)this << "\tadjacency ENUM('5p', '3p'), " << std::endl;
	*(std::ostream*)this << "\ttype ENUM('upward', 'downward', 'inward', 'outward'), " << std::endl;
	*(std::ostream*)this << "\tPRIMARY KEY(id)" << std::endl;
	*(std::ostream*)this << ");" << std::endl;
}

oSQLStream& oSQLStream::operator<< (const mccore::Residue& aResidue)
{
	muiResidueCount ++;
	mResIdMap.insert(model_res_id_pair(model_res_pair(muiModelCount, aResidue), muiResidueCount));
	*(std::ostream*)this << "INSERT INTO residues VALUES (" << muiResidueCount << ", ";
	*(std::ostream*)this << muiModelCount << ", ";
	*(std::ostream*)this << "\"" << aResidue.getResId().getChainId() << "\", ";
	*(std::ostream*)this << "\"" << aResidue.getResId().getResNo() << "\", ";
	*(std::ostream*)this << "\"" << aResidue.getResId().getInsertionCode() << "\", ";
	*(std::ostream*)this << "\"" << mccore::Pdbstream::stringifyResidueType (aResidue.getType ()) << "\", ";
	if (aResidue.getType ()->isNucleicAcid ())
	{
		*(std::ostream*)this << "\"" << aResidue.getPucker() << "\", ";
		*(std::ostream*)this << "\"" << aResidue.getGlycosyl() << "\");";
	}
	else
	{
		*(std::ostream*)this << "NULL, NULL);";
	}
	*(std::ostream*)this << std::endl;

	return *this;
}

oSQLStream& oSQLStream::operator<< (const AnnotateModel& aModel)
{
	muiModelCount ++;

	// Model identification
	*(std::ostream*)this << "INSERT INTO models VALUES (" << muiModelCount;
	*(std::ostream*)this << ", \"" << aModel.name();
	*(std::ostream*)this << "\", " << aModel.id() << ");" << std::endl;

	// Output the conformations
	AnnotateModel::const_iterator itRes;
	for (itRes = aModel.begin(); itRes != aModel.end (); ++itRes)
	{
		*this << *itRes;
		*(std::ostream*)this << std::endl;
	}

	return *this;
}

oSQLStream& oSQLStream::operator<< (const AnnotationInteractions& aInteractions)
{
	aInteractions.getSplitStacksAdjacencies();

	return *this;
}

oSQLStream& oSQLStream::flush()
{
	mResIdMap.clear();
	((std::ostream*)this)->flush();
	return *this;
}

}; // annotate

