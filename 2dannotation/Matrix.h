//                              -*- Mode: C++ -*- 
// Matrix.h
// Copyright © 2005-06 Laboratoire de Biologie Informatique et Théorique
//                     Université de Montréal
// Author           : Martin Larose
// Created On       : Wed Mar  9 16:10:51 2005
// $Revision: 88 $
// $Id: Matrix.h 88 2006-08-30 17:57:29Z larosem $
// 


#ifndef _mcsearch_Matrix_h_
#define _mcsearch_Matrix_h_

#include <vector>

#include <ostream>

namespace graphtest
{
  template< class Type >
  class Matrix : public std::vector< Type >
  {
    /**
     * The number of rows.
     */
    unsigned int rows;

    /**
     * The number of columns.
     */
    unsigned int columns;

  public:

    // LIFECYCLE ------------------------------------------------------------

    /**
     * Initializes the Matrix.
     */
    Matrix () : rows (0), columns (0) { }

    /**
     * Initializes the Matrix with row and column sizes.
     * @param rs the number of rows.
     * @param cs the number of columns.
     */
    Matrix (unsigned int rs, unsigned int cs)
      : std::vector< Type > (rs * cs),
	rows (rs),
	columns (cs)
    { }

    /**
     * Destroys the Matrix.
     */
    virtual ~Matrix () { }

    // OPERATORS ------------------------------------------------------------

    // ACCESS ---------------------------------------------------------------

    /**
     * Gets the value at position row column.
     * @param row the row index.
     * @param col the column index.
     * @return the value at row column indexes.
     */
    Type& get (unsigned int row, unsigned int col)
    {
      return std::vector< Type >::operator[] (row * columns + col);
    }

    /**
     * Gets the value at position row column.
     * @param row the row index.
     * @param col the column index.
     * @return the value at row column indexes.
     */
    const Type& get (unsigned int row, unsigned int col) const
    {
      return std::vector< Type >::operator[] (row * columns + col);
    }

    /**
     * Resizes the matrix for rows x cols size.
     * @param rows the numbers of rows.
     * @param cols the numbers of columns.
     */
    void resize (unsigned int rows, unsigned int cols)
    {
      this->rows = rows;
      columns = cols;
      std::vector< Type >::resize (rows * cols);
    }
    
    /**
     * Sets the value at position row column.
     * @param row the row index.
     * @param col the column index.
     * @param val the value to set.
     */
    void set (unsigned int row, unsigned int col, const Type &val)
    {
      std::vector< Type >::operator[] (row * columns + col) = val;
    }

    /**
     * Gets row size.
     * @return unsigned int
     */
    unsigned int getRowSize () const { return rows; }

    /**
     * Gets column size.
     * @return unsigned int
     */
    unsigned int getColumnSize () const { return columns; }
    
    // METHODS --------------------------------------------------------------

    // I/O  -----------------------------------------------------------------

    /**
     * Writes the Matrix to the output stream.
     * @param os the output stream.
     * @return the output stream.
     */
    virtual std::ostream& write (std::ostream &os) const
    {
      unsigned int i;

      os << "[Matrix]" << std::endl;
      for (i = 0; i < rows; ++i)
	{
	  unsigned int j;
	        
	  for (j = 0; j < columns; ++j)
	    {
	      os << std::vector< Type  >::operator[] (i * columns + j) << " ";
	    }
	  os << std::endl;
	}
      return os;
    }
    
  };

}



namespace std
{

  /**
   * Writes the Matrix to the output stream.
   * @param os the output stream.
   * @param obj the Matrix to write.
   * @return the output stream.
   */
  template< class Type >
  std::ostream& operator<< (std::ostream &os, const graphtest::Matrix< Type > &obj)
  {
    return obj.write (os);
  }
  
}

#endif

