/*
 * sqlstream.h
 *
 *  Created on: Mar 8, 2010
 *      Author: blanchmf
 */

#ifndef _annotate_SQLStream_h_
#define _annotate_SQLStream_h_

#include <iostream>
#include <map>

namespace mccore {
	class Residue;
};

namespace annotate {

class AnnotateModel;
class AnnotationInteractions;

class oSQLStream : public std::ostream
{
public:
	// LIFECYCLE -----------------------------------------------------------

	/**
	 * Initializes the stream.
	 */
	oSQLStream();

	/**
	 * Initializes the stream with a predefined stream buffer.
	 * @param sb the stream buffer.
	 */
	oSQLStream (std::streambuf* sb);

	oSQLStream (std::ostream &os);

	/**
	* Destroys the stream.
	*/
	virtual ~oSQLStream ();

	// OPERATORS -----------------------------------------------------------

	/**
	 * Casts the SQL stream to a ostream.
	 * @return the ostream.
	 */
	operator std::ostream* () { return dynamic_cast<std::ostream*>(this); }

    // I/O -----------------------------------------------------------------

	/**
	 * Writes an annotated model to the SQL stream.
	 * @param aModel the residue to write.
	 * @return itself.
	 */
	oSQLStream& operator<< (const AnnotateModel& aModel);

    /**
     * Writes a residue to the SQL stream.
     * @param aResidue the residue to write.
     * @return itself.
     */
	oSQLStream& operator<< (const mccore::Residue& aResidue);

	/**
	 * Synchronizes the buffer associated with the stream to its controlled
	 * output sequence. This methods also clears the mapping kept for objects in
	 * the stream. It is good practice to call this method between models when
	 * we're certain the information will no longer be required.
	 * @return itself.
	 */
	oSQLStream& flush();

private:
	// Number of models in the database
	unsigned int muiModelCount;
	unsigned int muiResidueCount;
	unsigned int muiInteractionCount;

	typedef std::pair<unsigned int, mccore::Residue> model_res_pair;
	typedef std::pair<model_res_pair, unsigned int > model_res_id_pair;
	std::map<model_res_pair, unsigned int> mResIdMap;

	oSQLStream& operator<< (const AnnotationInteractions& aInteractions);
};

}; // namespace annotate

#endif /* _annotate_SQLStream_h_ */
