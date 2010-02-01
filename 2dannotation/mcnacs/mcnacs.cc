//                              -*- Mode: C++ -*-
// mcnacs.cc
// Copyright © 2009 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Jul 29 10:28:00 2009

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mcnacs.h"
#include "NonAdjacentCycleStatsApp.h"
#include <iostream>

// PROTOTYPES ------------------------------------------------------------------

// FUNCTORS --------------------------------------------------------------------

// PROTOTYPES ------------------------------------------------------------------

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
int main (int argc, char *argv[])
{
	NonAdjacentCycleStatsApp theApp(argc, argv);

	// Get the cycles with their interactions
	std::set<NACycleInfo> connectedCycles = theApp.connectedCycles();

	// Remove the cycles not interacting with other NCMs
	theApp.filterOutOpenConnections();

	// Remove the cycles not covering entirely the relation
	theApp.filterOutPartialCoverage();

	// Remove the cycles who has subcycles in the same interactions
	// ( e.g. NGNRAN vs GNRA, keep GNRA only )
	theApp.filterOutEnclosingCycles();

	// Split the multiple adjacent cycles into multiple connections
	theApp.splitAdjacency();

	// Remove the cycles who connects to multibranch
	theApp.filterOutMultibranchCycles();

	// Remove the cycles who connects to multibranch
	theApp.filterOutLooseCycles();

	// Compile the statistics on interacting cycles
	theApp.compileStatistics();

	// Compute the interaction stats of the secondary structure cycles
	theApp.computeCyclesStats();

	// Output the statistics
	std::cout << theApp.toString();

	return EXIT_SUCCESS;
}
