#include "rnaselect2.h"

#include "AnnotateModel.h"

#include "mccore/Binstream.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/ResidueFactoryMethod.h"

#include <cassert>
#include <sstream>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
using namespace std;

// PROTOTYPE -------------------------------------------------------------------
mccore::Molecule loadMolecule(const string &filename, bool abBinary);
std::string getSequence(const mccore::Molecule& aMolecule);


// g++ -Wall -O5 rnaselect.C -o rnaselect.exe

//#define SHOW_INFO
#define SOURCE_BREAK  ":;( )$,-'"



//------------------------------------------------------------
class CDateData
{
public:
  CDateData( string & strDate )
  {
    if( !m_mapMonth2Int.size() )
      InitMap();

    m_iDay   = atoi( strDate.substr( 0, 2 ).c_str() );
    m_iMonth = m_mapMonth2Int[ strDate.substr( 3, 3 ).c_str() ];
    m_iYear  = atoi( strDate.substr( 7, 2 ).c_str() );
    if( m_iYear < 50 )
      m_iYear += 2000;
    else
      m_iYear += 1900;
  }

  int Compare( CDateData & pOther )
  {
    if( m_iYear > pOther.m_iYear )
      return +1;
    if( m_iYear == pOther.m_iYear )
      {
	if( m_iMonth > pOther.m_iMonth )
	  return +1;
	if( m_iMonth == pOther.m_iMonth )
	  {
	    if( m_iDay > pOther.m_iDay )
	      return +1;
	    if( m_iDay == pOther.m_iDay )
	      return 0;
	  }
      }

    return -1;
  }

  friend ostream & operator << ( ostream & out, CDateData & pData );

private:
  int m_iYear;
  int m_iMonth;
  int m_iDay;

  void InitMap( void )
  {
    m_mapMonth2Int[ "JAN" ] =  1;
    m_mapMonth2Int[ "FEB" ] =  2;
    m_mapMonth2Int[ "MAR" ] =  3;
    m_mapMonth2Int[ "APR" ] =  4;
    m_mapMonth2Int[ "MAY" ] =  5;
    m_mapMonth2Int[ "JUN" ] =  6;
    m_mapMonth2Int[ "JUL" ] =  7;
    m_mapMonth2Int[ "AUG" ] =  8;
    m_mapMonth2Int[ "SEP" ] =  9;
    m_mapMonth2Int[ "OCT" ] = 10;
    m_mapMonth2Int[ "NOV" ] = 11;
    m_mapMonth2Int[ "DEC" ] = 12;
  }

  static map< string, int > m_mapMonth2Int;
};

ostream & operator << ( ostream & out, CDateData & pData )
{
  int iYear = 0;
  if( pData.m_iYear >= 2000 )
    iYear = ( pData.m_iYear - 2000 );
  else
    iYear = ( pData.m_iYear - 1900 );

  if( iYear < 10 )
    out << "0";
  out << iYear;

  if( pData.m_iMonth < 10 )
    out << "0";
  out << pData.m_iMonth;

  if( pData.m_iDay < 10 )
    out << "0";
  out << pData.m_iDay;

  return out;
}

map< string, int > CDateData::m_mapMonth2Int;



//------------------------------------------------------------
class CDataNode
{
public:
  CDataNode( CDateData * pDate,
	     set< string > * pSource,
	     set< string > * pCompnd,
	     string & strSequence,
	     char * szBase,
	     char * szFile )
    : m_iColor( m_iNextColor++ ),
      m_pDate( pDate ),
      m_pSource( pSource ),
      m_pCompnd( pCompnd ),
      m_strSequence( strSequence ),
      m_strBase( szBase ),
      m_strFile( szFile )
  {
  }

  int             m_iColor;
  CDateData     * m_pDate;
  set< string > * m_pSource;
  set< string > * m_pCompnd;
  string          m_strSequence;
  string          m_strBase;
  string          m_strFile;

  static int m_iNextColor;
};

int CDataNode::m_iNextColor = 0;

map< string, CDataNode * > mapFile2Data;



//------------------------------------------------------------
#define MATCHSCORE     +10
#define MISMATCHSCORE  -7
#define GAPSCORE       -4
float SmithWaterman( const char * szSeq1, const char * szSeq2 )
{
  /* initialization */
  int i, j, dist, tmp;
  int m = strlen( szSeq1 );
  int n = strlen( szSeq2 );

  int ** distance = new int * [ m+1 ];
  for( i = 0; i <= m; ++i )
    distance[i] = new int[ n+1 ];

  for(i=0;i<=m;i++) distance[i][0]=0; /* do this for i=0,1,...,m */
  for(j=0;j<=n;j++) distance[0][j]=0;

  int minDist = 0;
  int minI = 0;
  int minJ = 0;

  for(i=1;i<=m;i++)
    {    /* note: we begin at i=1 ! */
      for(j=1;j<=n;j++)
	{
	  dist=0; /* distance to node (i,j) from virtual start node */

	  if(szSeq1[i-1]==szSeq2[j-1])
	    tmp = distance[i-1][j-1]-MATCHSCORE;
	  else
	    tmp = distance[i-1][j-1]-MISMATCHSCORE;
	  if(tmp<dist)
	    {
	      dist=tmp;
	    }

	  tmp = distance[i-1][j]-GAPSCORE;
	  if(tmp<dist)
	    {
	      dist=tmp;
	    }

	  tmp = distance[i][j-1]-GAPSCORE;
	  if(tmp<dist)
	    {
	      dist=tmp;
	    }

	  distance[i][j]=dist;

	  if(dist<minDist)
	    { /* keep track of where the minimum score is */
	      minDist=dist;
	      minI=i;
	      minJ=j;
	    }
	}
    }


  for( i = 0; i <= m; ++i )
    {
      /*
      for( j = 0; j <= n; ++j )
	{
	  char szBuf[10];
	  sprintf( szBuf, "%3d ", distance[i][j] );
	  cerr << szBuf;
	}
      cerr << "\n";
      */

      delete [] distance[i];
    }
  delete [] distance;


  return 100.0 * (float)minDist / ( -MATCHSCORE * max(m,n) );
}



//------------------------------------------------------------
void _splitpath( char *name,
		 char *drive, char *path, char *base, char *ext )
{
  char *s, *p;

  p = name;
  s = strchr(p, ':');
  if ( s != NULL )
    {
      if (drive)
	{
	  *s = '\0';
	  strcpy(drive, p);
	  *s = ':';
	}
      p = s+1;
      if (!p)
	return;
    }
  else
    if (drive)
      *drive = '\0';

  s = strrchr(p, '/');
  if ( s != NULL)
    {
      if (path)
	{
	  char c;

	  c = *(s+1);
	  *(s+1) = '\0';
	  strcpy(path, p);
	  *(s+1) = c;
	}
      p = s+1;
      if (!p)
	return;
    }
  else
    if (path)
      *path = '\0';

  s = strchr(p, '.');
  if ( s != NULL)
    {
      if (base)
	{
	  *s = '\0';
	  strcpy(base, p);
	  *s = '.';
	}
      p = s+1;
      if (!p)
	return;
    }
  else
    if (base)
      *base = '\0';

  if (ext)
    strcpy(ext, p);
}



//----------------------------------------------------------------------
bool IsNumber( char * szToken )
{
  int iMax = strlen( szToken );
  bool bIsNum = true;
  for( int i = 0; bIsNum && ( i < iMax ); ++i )
    bIsNum = ( bIsNum && isdigit( (int)szToken[ i ] ) );

  return bIsNum;
}



//----------------------------------------------------------------------
void FillSet( string & strBuffer,
	      set< string > & setString,
	      const char * szBreak )
{
  char * szBuffer = strdup( strBuffer.c_str() );
  char * szStart  = strstr( szBuffer, ":" );
  if( NULL == szStart )
    szStart = szBuffer + 10;

  char * szToken = strtok( szStart, SOURCE_BREAK );
  while( NULL != szToken )
    {
      if( ( strlen( szToken ) > 1 ) &&
	  ( strstr( szToken, "*" ) == NULL ) &&
	  !( strstr( szToken, "OF"      ) == szToken ) &&
	  !( strstr( szToken, "RNA"     ) == szToken ) &&
	  !( strstr( szToken, "DNA"     ) == szToken ) &&
	  !( strstr( szToken, "AND"     ) == szToken ) &&
	  !( strstr( szToken, "YES"     ) == szToken ) &&
	  !( strstr( szToken, "WITH"    ) == szToken ) &&
	  !( strstr( szToken, "PROTEIN" ) == szToken ) &&
	  !IsNumber( szToken )
	  )
	{
	  setString.insert( szToken );
	}

      szToken = strtok( NULL, SOURCE_BREAK );
    }

  free( szBuffer );
}



//----------------------------------------------------------------------
void FillSequence( string & strBuffer, string & strSequence )
{
  char * szBuffer = strdup( strBuffer.c_str() );

  char * szToken = strtok( szBuffer+18, " " );
  while( NULL != szToken )
    {
      if( strlen( szToken ) == 1 )
	{
	  if( ( *szToken == 'A' ) ||
	      ( *szToken == 'C' ) ||
	      ( *szToken == 'G' ) ||
	      ( *szToken == 'U' ) )
	    strSequence += szToken;
	}

      szToken = strtok( NULL, " " );
    }

  free( szBuffer );
}



//----------------------------------------------------------------------
void LoadFile( char * szFile )
{
	static int iFileNo = 0;

	// open the file
	ifstream hFile( szFile );
	if( !hFile )
	{
		std::cerr << "Failed to open PDB file: `" << szFile << "`\n";
		return;
	}


	char szBase[256];
	_splitpath( szFile, NULL, NULL, szBase, NULL );
	std::cerr << szBase << "\t";
	if( ( ++iFileNo % 10 ) == 0 )
		cerr << "\n";


	CDateData     * pDate   = NULL;
	set< string > * pCompnd = new set< string >;
	set< string > * pSource = new set< string >;


	// walk the file line by line
	string strBuffer;
	string strSequence;
	bool bRevDate = false;
	while( getline( hFile, strBuffer ) )
	{
		// get latest revision date
		if( !bRevDate && ( strBuffer.find( "REVDAT" ) == 0 ) )
		{
			bRevDate = true;

			string strDate = strBuffer.substr( 13, 9 );
			pDate = new CDateData( strDate );
		}
		else if( ( strBuffer.find( "SOURCE"  ) == 0 ) && ( strBuffer.find( "MOL_ID:" ) == string::npos ) )
		{
			// split the source line into pieces
			FillSet( strBuffer, *pSource, SOURCE_BREAK );
		}
		else if( ( strBuffer.find( "COMPND"  ) == 0 ) && ( strBuffer.find( "MOL_ID:" ) == string::npos ) && ( strBuffer.find( "CHAIN:"  ) == string::npos ) )
		{
			// split the source line into pieces
			FillSet( strBuffer, *pCompnd, SOURCE_BREAK );
		}
		else if( strBuffer.find( "SEQRES"  ) == 0 )
		{
			FillSequence( strBuffer, strSequence );
		}
    }

	// close the file
	hFile.close();


	// Load the sequence from the residues
	mccore::Molecule mol = loadMolecule(szFile, false);
	strSequence = getSequence(mol);

  // insert the data under proper file
	mapFile2Data[ szBase ] = new CDataNode( pDate, pSource, pCompnd, strSequence, szBase, szFile );
}



//----------------------------------------------------------------------
set< string > SetIntersection( set< string > & set1,
			       set< string > & set2 )
{
  set< string > setInter;
  insert_iterator< set< string > > iter( setInter, setInter.begin() );

  set_intersection( set1.begin(), set1.end(),
		    set2.begin(), set2.end(),
		    iter );

  return setInter;
}



//----------------------------------------------------------------------
int GetSetSimilarity( set< string > & set1,
		      set< string > & set2 )
{
  set< string > setInter = SetIntersection( set1, set2 );
  int iMaxSize = max( set1.size(), set2.size() );
  return( ( iMaxSize > 0 ) ? ( 100 * setInter.size() ) / iMaxSize : 100 );
}



//----------------------------------------------------------------------
int GetSimilarity( CDataNode * pData1, CDataNode * pData2 )
{
  if( NULL == pData1 )
    return 0;

  if( NULL == pData2 )
    return 0;
/*
  if( NULL == pData1->m_pSource )
    return 0;

  if( NULL == pData2->m_pSource )
    return 0;

  if( NULL == pData1->m_pCompnd )
    return 0;

  if( NULL == pData2->m_pCompnd )
    return 0;

  if( NULL == pData1->m_pDate )
    return 0;

  if( NULL == pData2->m_pDate )
    return 0;
*/
  int iSimil = 0;

  /*
  int iSimil = GetSetSimilarity
    ( *( pData1->m_pCompnd ),
      *( pData2->m_pCompnd ) );

  if( iSimil >= 50 )
  */
    {
      /* NOT USING SOURCE
      iSimil = GetSetSimilarity
	( *( pData1->m_pSource ),
	  *( pData2->m_pSource ) );
      */
  /*
      if( iSimil >= 50 )
  */
	{
          float fSW = SmithWaterman
	    ( pData1->m_strSequence.c_str(),
	      pData2->m_strSequence.c_str() );

	  if( fSW > 95.0 )
	    iSimil = 100;
	  else
	    iSimil = 0;
	}
    }

    return iSimil;
}

//----------------------------------------------------------------------
void ComputeGroups( void )
{
	cerr << "Computing Groups (one moment please):\n";
	map< string, CDataNode * >::iterator pos1;
	map< string, CDataNode * >::iterator pos2;
	for( pos1 = mapFile2Data.begin(); pos1 != mapFile2Data.end(); ++pos1 )
	{
		pos2 = pos1;
		++pos2;
		for( ; pos2 != mapFile2Data.end(); ++pos2 )
		{
			if( GetSimilarity( pos1->second, pos2->second ) >= 50 )
			{
				int iMinColor = min( pos1->second->m_iColor,
				pos2->second->m_iColor );
				pos1->second->m_iColor = iMinColor;
				pos2->second->m_iColor = iMinColor;
			}
		}
	}
}

//----------------------------------------------------------------------
void ShowGroups( void )
{
	int iNewColor = 1;
	vector< string > vecRemove;
	map< string, CDataNode * >::iterator pos1;
	for( int iColor = 0; iColor < CDataNode::m_iNextColor; ++iColor )
	{
		bool bShowColor = false;
		vector< CDataNode * > vecSequences;

		for( pos1 = mapFile2Data.begin(); pos1 != mapFile2Data.end(); ++pos1 )
		{
			if( iColor == pos1->second->m_iColor )
			{
				if( !bShowColor )
				{
					bShowColor = true;
					cerr << "\nGroup: " << iNewColor++ << "\n";
				}

				cerr << "\t " << pos1->first;

				vecSequences.push_back( pos1->second );
			}
		}
	}

	cerr << "\nGroups: " << (iNewColor-1)
		<< " on Total: " << CDataNode::m_iNextColor << "\n";
}

mccore::Molecule loadMolecule(const string &filename, bool abBinary)
{
	mccore::ResidueFM rFM;
	annotate::AnnotateModelFM aFM (mccore::ResIdSet(), 0, &rFM);
	mccore::Molecule molecule(&aFM);

	if (abBinary)
	{
		mccore::izfBinstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			std::ostringstream oss;
			oss << PACKAGE << ": cannot open binary file '" << filename << "'.";
			throw mccore::FileNotFoundException(oss.str(), __FILE__, __LINE__);
		}
		in >> molecule;
		in.close ();
	}
	else
	{
		mccore::izfPdbstream in;

		in.open (filename.c_str ());
		if (in.fail ())
		{
			std::ostringstream oss;
			oss << PACKAGE << ": cannot open pdb file '" << filename << "'.";
			throw mccore::FileNotFoundException(oss.str(), __FILE__, __LINE__);
		}
		in >> molecule;
		in.close ();
	}
	return molecule;
}

std::string getSequence(const mccore::Molecule& aMolecule)
{
	assert(1 == aMolecule.size());
	std::ostringstream oss;
	mccore::Molecule::const_iterator molIt = aMolecule.begin();
	mccore::GraphModel::const_iterator it;
	for(it = molIt->begin(); it != molIt->end(); ++ it)
	{
		oss << mccore::Pdbstream::stringifyResidueType (it->getType ());
	}
	return oss.str();
}

//----------------------------------------------------------------------
int main( int argc, char * argv[] )
{
  // check args
  if( argc < 2 )
    {
      cerr << "Usage:\n\t" << argv[0] << " /here/and/there/*.pdb\n";
      return -1;
    }


  // load the PDB file:
  //   get molecule
  //   get source
  //   get latest revision date
  cerr << "Loading Files:\n";
  for( int i = 1; i < argc; ++i )
    LoadFile( argv[i] );
  cerr << "\n\n";


  // if( mol1 != mol2 )
  //  then different
  // else
  // {
  //   // same molecule!
  //   if( source1 != source2 )
  //    then different
  //   else
  //   {
  //     // same source also!
  //     pick latest date
  //     if same then pick longest sequence
  //   }
  // }
  ComputeGroups();
  ShowGroups();

  return 0;
}
