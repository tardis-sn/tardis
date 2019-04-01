/*
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 */

/**\file Matrix3.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-18
 *\date 2012-08-2 Updated by rhl to be more generic. less code that does more!
 * TODO: Remove all 'inline' keywords as it is only a suggestion to the compiler
 * anyways, and it will ignore it or add it when it thinks its necessary.
 */

#ifndef MOAB_MATRIX3_HPP
#define MOAB_MATRIX3_HPP

#include <iostream>
#include <iosfwd>
#include <limits>
#include <cmath>
#include <cassert>

#include "moab/Util.hpp"
#include "moab/Types.hpp"
#include "moab/CartVect.hpp"

namespace moab {

namespace Matrix{
	template< typename Matrix>
	Matrix inverse( const Matrix & d, const double i){
		Matrix m( d);
		 m( 0) = i * (d(4) * d(8) - d(5) * d(7));
	         m( 1) = i * (d(2) * d(7) - d(8) * d(1));
	         m( 2) = i * (d(1) * d(5) - d(4) * d(2));
	         m( 3) = i * (d(5) * d(6) - d(8) * d(3));
	         m( 4) = i * (d(0) * d(8) - d(6) * d(2));
	         m( 5) = i * (d(2) * d(3) - d(5) * d(0));
	         m( 6) = i * (d(3) * d(7) - d(6) * d(4));
	         m( 7) = i * (d(1) * d(6) - d(7) * d(0));
	         m( 8) = i * (d(0) * d(4) - d(3) * d(1));
		return m;
	}

	template< typename Matrix>
	inline bool positive_definite( const Matrix & d, 
				       double& det ){
	        double subdet6 = d(1)*d(5)-d(2)*d(4);
	        double subdet7 = d(2)*d(3)-d(0)*d(5);
	        double subdet8 = d(0)*d(4)-d(1)*d(3);
	        det = d(6)*subdet6 + d(7)*subdet7 + d(8)*subdet8;
	        return d(0) > 0 && subdet8 > 0 && det > 0;
	}
	template< typename Matrix>
	inline Matrix transpose( const Matrix & d){
	      return Matrix( d(0), d(3), d(6),
	                     d(1), d(4), d(7),
	                     d(2), d(5), d(8) );
	}
	template< typename Matrix>
	inline Matrix mmult3( const Matrix& a, const Matrix& b ) {
	  return Matrix( a(0,0) * b(0,0) + a(0,1) * b(1,0) + a(0,2) * b(2,0),
	                 a(0,0) * b(0,1) + a(0,1) * b(1,1) + a(0,2) * b(2,1),
	                 a(0,0) * b(0,2) + a(0,1) * b(1,2) + a(0,2) * b(2,2),
	                 a(1,0) * b(0,0) + a(1,1) * b(1,0) + a(1,2) * b(2,0),
	                 a(1,0) * b(0,1) + a(1,1) * b(1,1) + a(1,2) * b(2,1),
	                 a(1,0) * b(0,2) + a(1,1) * b(1,2) + a(1,2) * b(2,2),
	                 a(2,0) * b(0,0) + a(2,1) * b(1,0) + a(2,2) * b(2,0),
	                 a(2,0) * b(0,1) + a(2,1) * b(1,1) + a(2,2) * b(2,1),
	                 a(2,0) * b(0,2) + a(2,1) * b(1,2) + a(2,2) * b(2,2) );
	}

	template< typename Vector, typename Matrix>
	inline Matrix outer_product( const Vector & u,
	                              const Vector & v,
				      Matrix & m ) {
	  	m = Matrix( u[0] * v[0], u[0] * v[1], u[0] * v[2],
	                    u[1] * v[0], u[1] * v[1], u[1] * v[2],
	                    u[2] * v[0], u[2] * v[1], u[2] * v[2] );
		return m;
	}
	template< typename Matrix>
	inline double determinant3( const Matrix & d){
		return (d(0) * d(4) * d(8) 
		     + d(1) * d(5) * d(6)
		     + d(2) * d(3) * d(7)
		     - d(0) * d(5) * d(7)
		     - d(1) * d(3) * d(8)
		     - d(2) * d(4) * d(6)); 
	}
	
	template< typename Matrix>
	inline const Matrix inverse( const Matrix & d){
		const double det = 1.0/determinant3( d);
		return inverse( d, det);
	}
	
	template< typename Vector, typename Matrix>
	inline Vector vector_matrix( const Vector& v, const Matrix& m ) {
	  return Vector( v[0] * m(0,0) + v[1] * m(1,0) + v[2] * m(2,0),
	                 v[0] * m(0,1) + v[1] * m(1,1) + v[2] * m(2,1),
	                 v[0] * m(0,2) + v[1] * m(1,2) + v[2] * m(2,2) );
	}
	
	template< typename Vector, typename Matrix>
	inline Vector matrix_vector( const Matrix& m, const Vector& v ){
	   Vector res = v;
	   res[ 0] = v[0] * m(0,0) + v[1] * m(0,1) + v[2] * m(0,2);
	   res[ 1] = v[0] * m(1,0) + v[1] * m(1,1) + v[2] * m(1,2);
	   res[ 2] = v[0] * m(2,0) + v[1] * m(2,1) + v[2] * m(2,2);
	   return res;
	} 

	// moved from EigenDecomp.hpp

	// Jacobi iteration for the solution of eigenvectors/eigenvalues of a nxn
	// real symmetric matrix. Square nxn matrix a; size of matrix in n;
	// output eigenvalues in w; and output eigenvectors in v. Resulting
	// eigenvalues/vectors are sorted in decreasing order; eigenvectors are
	// normalized.
	// TODO: Remove this
	#define VTK_ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	        a[k][l]=h+s*(g-h*tau)

	//TODO: Refactor this method into subroutines
	//use a namespace { }  with no name to
	//contain subroutines so that the compiler
	//automatically inlines them.

	template< typename Matrix, typename Vector>
	ErrorCode EigenDecomp( const Matrix & _a,
	                       double w[3],
	                       Vector o[3] ) {
	  Vector v[3];
	  const int MAX_ROTATIONS = 20;
	  const double one_ninth = 1./9;
	  int i, j, k, iq, ip, numPos;
	  double tresh, theta, tau, t, sm, s, h, g, c, tmp;
	  double b[3], z[3];
	  Matrix a( _a);

	  // initialize
	  for (ip=0; ip<3; ip++) {
	    for (iq=0; iq<3; iq++){
	      v[ip][iq] = 0.0;
	    }
	    v[ip][ip] = 1.0;
	  }
	  for (ip=0; ip<3; ip++) {
	    b[ip] = w[ip] = a[ip][ip];
	    z[ip] = 0.0;
	  }

	  // begin rotation sequence
	  for (i=0; i<MAX_ROTATIONS; i++){
	    sm = 0.0;
	    for (ip=0; ip<2; ip++){
	      for (iq=ip+1; iq<3; iq++){ sm += fabs(a[ip][iq]); }
	    }

	    if ( sm == 0.0 ){ break; }
	    // first 3 sweeps
	    tresh = (i < 3)? 0.2*sm*one_ninth : 0.0;
	    for (ip=0; ip<2; ip++) {
	      for (iq=ip+1; iq<3; iq++) {
	        g = 100.0*fabs(a[ip][iq]);

	        // after 4 sweeps
	        if ( i > 3 && (fabs(w[ip])+g) == fabs(w[ip])
	             && (fabs(w[iq])+g) == fabs(w[iq])) {
	          a[ip][iq] = 0.0;
	  }
	        else if ( fabs(a[ip][iq]) > tresh) {
	          h = w[iq] - w[ip];
	          if ( (fabs(h)+g) == fabs(h)){ t = (a[ip][iq]) / h; }
	          else {
	            theta = 0.5*h / (a[ip][iq]);
	            t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
	            if (theta < 0.0) { t = -t;}
	          }
	          c = 1.0 / sqrt(1+t*t);
	          s = t*c;
	          tau = s/(1.0+c);
	          h = t*a[ip][iq];
	          z[ip] -= h;
	          z[iq] += h;
	          w[ip] -= h;
	          w[iq] += h;
	          a[ip][iq]=0.0;
	          // ip already shifted left by 1 unit
	          for (j = 0;j <= ip-1;j++) { VTK_ROTATE(a,j,ip,j,iq); }
	          // ip and iq already shifted left by 1 unit
	          for (j = ip+1;j <= iq-1;j++) { VTK_ROTATE(a,ip,j,j,iq); }
	          // iq already shifted left by 1 unit
	          for (j=iq+1; j<3; j++) { VTK_ROTATE(a,ip,j,iq,j); }
	          for (j=0; j<3; j++) { VTK_ROTATE(v,j,ip,j,iq); }
	          }
	        }
	      }

	    for (ip=0; ip<3; ip++) {
	      b[ip] += z[ip];
	      w[ip] = b[ip];
	      z[ip] = 0.0;
	    }
	  }

	  //// this is NEVER called
	  if ( i >= MAX_ROTATIONS ) {
	      std::cerr << "Matrix3D: Error extracting eigenfunctions" << std::endl;
	      return MB_FAILURE;
	  }

	  // sort eigenfunctions                 these changes do not affect accuracy
	  for (j=0; j<2; j++){                  // boundary incorrect
	    k = j;
	    tmp = w[k];
	    for (i=j+1; i<3; i++){                // boundary incorrect, shifted already
	      if (w[i] >= tmp){                  // why exchage if same?
	        k = i;
	        tmp = w[k];
	        }
	    }
	    if (k != j){
	      w[k] = w[j];
	      w[j] = tmp;
	      for (i=0; i<3; i++){
	        tmp = v[i][j];
	        v[i][j] = v[i][k];
	        v[i][k] = tmp;
	        }
	    }
	  }
	  // insure eigenvector consistency (i.e., Jacobi can compute vectors that
	  // are negative of one another (.707,.707,0) and (-.707,-.707,0). This can
	  // reek havoc in hyperstreamline/other stuff. We will select the most
	  // positive eigenvector.
	  int ceil_half_n = (3 >> 1) + (3 & 1);
	  for (j=0; j<3; j++) {
	    for (numPos=0, i=0; i<3; i++) {
	      if ( v[i][j] >= 0.0 ) { numPos++; }
	    }
	//    if ( numPos < ceil(double(n)/double(2.0)) )
	    if ( numPos < ceil_half_n) {
	      for(i=0; i<3; i++) { v[i][j] *= -1.0; }
	    }
	  }

	  //transpose the vector array for output
	  for(i=0; i<3; i++)
	    {
	    for(j=0; j<3; j++)
	      { o[i][j]=v[j][i]; }
	    }

	  return MB_SUCCESS;
	}
} //namespace Matrix

class Matrix3  {
  //TODO: std::array when we can use C++11
  double d[9];

public:
  //Default Constructor
  inline Matrix3(){
	for(int i = 0; i < 9; ++i){ d[ i] = 0; }
  }
  //TODO: Deprecate this.
  //Then we can go from three Constructors to one. 
  inline Matrix3( double diagonal ){ 
      d[0] = d[4] = d[8] = diagonal;
      d[1] = d[2] = d[3] = 0.0;
      d[5] = d[6] = d[7] = 0.0;
  }
  inline Matrix3( const CartVect & diagonal ){ 
      d[0] = diagonal[0];
      d[4] = diagonal[1],
      d[8] = diagonal[2];
      d[1] = d[2] = d[3] = 0.0;
      d[5] = d[6] = d[7] = 0.0;
  }
  //TODO: not strictly correct as the Matrix3 object
  //is a double d[ 9] so the only valid model of T is
  //double, or any refinement (int, float) 
  //*but* it doesn't really matter anything else
  //will fail to compile.
  template< typename T> 
  inline Matrix3( const std::vector< T> & diagonal ){ 
      d[0] = diagonal[0];
      d[4] = diagonal[1],
      d[8] = diagonal[2];
      d[1] = d[2] = d[3] = 0.0;
      d[5] = d[6] = d[7] = 0.0;
  }

inline Matrix3( double v00, double v01, double v02,
                double v10, double v11, double v12,
                double v20, double v21, double v22 ){
    d[0] = v00; d[1] = v01; d[2] = v02;
    d[3] = v10; d[4] = v11; d[5] = v12;
    d[6] = v20; d[7] = v21; d[8] = v22;
}

  //Copy constructor 
  Matrix3 ( const Matrix3 & f){
	for(int i = 0; i < 9; ++i) { d[ i] = f.d[ i]; }
  }
  //Weird constructors 
  template< typename Vector> 
  inline Matrix3(   const Vector & row0,
                    const Vector & row1,
                    const Vector & row2 ) {
      for(std::size_t i = 0; i < 3; ++i){
	d[ i] = row0[ i];
	d[ i+3]= row1[ i];
	d[ i+6] = row2[ i];
      }
  }
  
  inline Matrix3( const double* v ){ 
      d[0] = v[0]; d[1] = v[1]; d[2] = v[2];
      d[3] = v[3]; d[4] = v[4]; d[5] = v[5]; 
      d[6] = v[6]; d[7] = v[7]; d[8] = v[8];
  }
  
  inline Matrix3& operator=( const Matrix3& m ){
      d[0] = m.d[0]; d[1] = m.d[1]; d[2] = m.d[2];
      d[3] = m.d[3]; d[4] = m.d[4]; d[5] = m.d[5];
      d[6] = m.d[6]; d[7] = m.d[7]; d[8] = m.d[8];
      return *this;
  }
  
  inline Matrix3& operator=( const double* v ){ 
      d[0] = v[0]; d[1] = v[1]; d[2] = v[2];
      d[3] = v[3]; d[4] = v[4]; d[5] = v[5]; 
      d[6] = v[6]; d[7] = v[7]; d[8] = v[8];
      return *this;
 }

  inline double* operator[]( unsigned i ){ return d + 3*i; }
  inline const double* operator[]( unsigned i ) const{ return d + 3*i; }
  inline double& operator()(unsigned r, unsigned c) { return d[3*r+c]; }
  inline double operator()(unsigned r, unsigned c) const { return d[3*r+c]; }
  inline double& operator()(unsigned i) { return d[i]; }
  inline double operator()(unsigned i) const { return d[i]; }
  
    // get pointer to array of nine doubles
  inline double* array()
      { return d; }
  inline const double* array() const
      { return d; }

  inline Matrix3& operator+=( const Matrix3& m ){
      d[0] += m.d[0]; d[1] += m.d[1]; d[2] += m.d[2];
      d[3] += m.d[3]; d[4] += m.d[4]; d[5] += m.d[5];
      d[6] += m.d[6]; d[7] += m.d[7]; d[8] += m.d[8];
      return *this;
  }
  
  inline Matrix3& operator-=( const Matrix3& m ){
      d[0] -= m.d[0]; d[1] -= m.d[1]; d[2] -= m.d[2];
      d[3] -= m.d[3]; d[4] -= m.d[4]; d[5] -= m.d[5];
      d[6] -= m.d[6]; d[7] -= m.d[7]; d[8] -= m.d[8];
      return *this;
  }
  
  inline Matrix3& operator*=( double s ){
      d[0] *= s; d[1] *= s; d[2] *= s;
      d[3] *= s; d[4] *= s; d[5] *= s;
      d[6] *= s; d[7] *= s; d[8] *= s;
      return *this;
 }
  
  inline Matrix3& operator/=( double s ){
      d[0] /= s; d[1] /= s; d[2] /= s;
      d[3] /= s; d[4] /= s; d[5] /= s;
      d[6] /= s; d[7] /= s; d[8] /= s;
      return *this;
  }
 
  inline Matrix3& operator*=( const Matrix3& m ){
	(*this) = moab::Matrix::mmult3((*this),m); 
	return *this;
  }
  
  inline double determinant() const{
  	return moab::Matrix::determinant3( *this);
  }
 
  inline Matrix3 inverse() const { 
	const double i = 1.0/determinant();
	return moab::Matrix::inverse( *this, i); 
  }
  inline Matrix3 inverse( double i ) const {
  	return moab::Matrix::inverse( *this, i); 
  }
  
  inline bool positive_definite() const{
  	double tmp;
  	return positive_definite( tmp);
  }
  
  inline bool positive_definite( double& det ) const{
	  return moab::Matrix::positive_definite( *this, det);
  }
  
  inline Matrix3 transpose() const{ return moab::Matrix::transpose( *this); }
  
  inline bool invert() {
    double i = 1.0 / determinant();
    if (!Util::is_finite(i) || fabs(i) < std::numeric_limits<double>::epsilon())
      return false;
    *this = inverse( i );
    return true;
  }
    // Calculate determinant of 2x2 submatrix composed of the
    // elements not in the passed row or column.
  inline double subdet( int r, int c ) const{
	const int r1 = (r+1)%3, r2 = (r+2)%3;
	const int c1 = (c+1)%3, c2 = (c+2)%3;
	assert(r >= 0 && c >= 0);
	if (r < 0 || c < 0) return DBL_MAX;
	return d[3*r1+c1]*d[3*r2+c2] - d[3*r1+c2]*d[3*r2+c1];
  }
}; //class Matrix3

inline Matrix3 operator+( const Matrix3& a, const Matrix3& b ){ 
	return Matrix3(a) += b; 
}
inline Matrix3 operator-( const Matrix3& a, const Matrix3& b ){ 
	return Matrix3(a) -= b; 
}

inline Matrix3 operator*( const Matrix3& a, const Matrix3& b ) {
	return moab::Matrix::mmult3( a, b);
}

template< typename Vector>
inline Matrix3 outer_product( const Vector & u,
                              const Vector & v ) {
  return Matrix3( u[0] * v[0], u[0] * v[1], u[0] * v[2],
                  u[1] * v[0], u[1] * v[1], u[1] * v[2],
                  u[2] * v[0], u[2] * v[1], u[2] * v[2] );
}

template< typename T>
inline std::vector< T> operator*( const Matrix3&m, const std::vector< T> & v){
		return moab::Matrix::matrix_vector( m, v);
}

template< typename T>
inline std::vector< T> operator*( const std::vector< T>& v, const Matrix3&m){
		return moab::Matrix::vector_matrix( v, m);
}

inline CartVect operator*( const Matrix3&m,  const CartVect& v){
		return moab::Matrix::matrix_vector( m, v);
}

inline CartVect operator*( const CartVect& v, const Matrix3& m){
		return moab::Matrix::vector_matrix( v, m);
}

} // namespace moab

#ifndef MOAB_MATRIX3_OPERATORLESS
#define MOAB_MATRIX3_OPERATORLESS
inline std::ostream& operator<<( std::ostream& s, const moab::Matrix3& m ){
  return s <<  "| " << m(0,0) << " " << m(0,1) << " " << m(0,2) 
           << " | " << m(1,0) << " " << m(1,1) << " " << m(1,2) 
           << " | " << m(2,0) << " " << m(2,1) << " " << m(2,2) 
           << " |" ;
}
#endif//MOAB_MATRIX3_OPERATORLESS
#endif //MOAB_MATRIX3_HPP
