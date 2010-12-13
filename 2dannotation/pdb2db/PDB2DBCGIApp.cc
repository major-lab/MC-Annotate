/*
 * PDB2DBCGIApp.cc
 *
 *  Created on: Jul 21, 2010
 *      Author: blanchmf
 */

#include "PDB2DBCGIApp.h"

#include "cgicc/Cgicc.h"
#include "cgicc/HTTPHTMLHeader.h"
#include "cgicc/HTMLClasses.h"

#include <sstream>
#include <string>

PDB2DotBracketParamsCGI::PDB2DotBracketParamsCGI() : PDB2DotBracketParams()
{
}

int main(int argc, char* argv[])
{
	PDB2DotBracketParamsCGI params;
	try
	{
		cgicc::Cgicc cgi;

		// Send HTTP header
		std::cout << cgicc::HTTPHTMLHeader() << std::endl;

		std::cout << "<!-- IE Quirks Mode -->" << std::endl;

		// Set up the HTML document
		std::cout << cgicc::html() << cgicc::head(cgicc::title("PDB to Dot Bracket")) << std::endl;
		std::cout << cgicc::body() << std::endl;

		std::cout << "<pre>" << std::endl;

		cgicc::form_iterator itElem = cgi.getElement("nbcombinedlayers");
		if(itElem != cgi.getElements().end())
		{
			params.muiCombinedLayers = itElem->getIntegerValue();
		}

		itElem = cgi.getElement("nbsplitlayers");
		if(itElem != cgi.getElements().end())
		{
			params.muiSplitLayers = itElem->getIntegerValue();
		}

		itElem = cgi.getElement("completegaps");
		if(itElem != cgi.getElements().end())
		{
			params.mbCompleteGaps = (itElem->getIntegerValue() == 1);
		}


		// Print out the submitted files
		cgicc::file_iterator referencepdb = cgi.getFile("referencepdb");
		if(referencepdb != cgi.getFiles().end())
		{
			// Write to a work stream
			std::stringbuf stringBuffer;
			std::ostream pdbStream(&stringBuffer);
			referencepdb->writeToStream(pdbStream);
			std::string strFileName = referencepdb->getFilename();
			PDB2DotBracket theApp(
				params,
				strFileName,
				pdbStream);
		}

		std::cout << "</pre>" << std::endl;

		// Close the HTML document
		std::cout << cgicc::body() << cgicc::html();
	}
	catch(std::exception& e)
	{
		// handle any errors - omitted for brevity
	}
	return EXIT_SUCCESS;
}
