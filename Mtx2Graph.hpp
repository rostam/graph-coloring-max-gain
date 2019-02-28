#ifndef MTX2BIPGRAPH_HPP
#define MTX2BIPGRAPH_HPP

#include <iostream>
#include "mmio.h"
#include <cstdio>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#endif

/** 
 * \struct matrix_market
 * \brief Convert graph of input file (MM format) to internal 
 * data representation
 *
 * This procedure opens a input file. This file contains a matrix in
 * MM-format. This matrix is converted into the interal data structure.
 */
struct matrix_market
{
  matrix_market(const char* filename);
  boost::numeric::ublas::matrix<int> to_ublas_matrix();
  ~matrix_market();
  bool write_to_file(char* filename);

  inline unsigned int nrows() const { return M; }
  inline unsigned int issym() const { return mm_is_symmetric(matcode); }

  int M;
  int N;
  int nz;
  unsigned int *I;
  unsigned int *J;
  MM_typecode matcode;
};
