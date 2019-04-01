/**
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
 * 
 */

#ifndef MOAB_HOMXFORM
#define MOAB_HOMXFORM

/*
 * \file HomXform.hpp
 *
 * \brief Representation and functions for homogeneous transforms
 *
 * A homogeneous transform is a 4x4 matrix for representing and manipulating
 * homogeneous coordinates, which are x' = (x/h, y/h, z/h, h).  See Mortenson,
 * Geometric Modeling, or other texts on geometric modeling for details of homogeneous
 * transforms.
 */

#define XFORM(a,b) xForm[4*a+b]
#define XFORM_INDEX(a,b) 4*a+b

#include <math.h>
#include <ostream>

namespace moab {

class HomXform;

/** \class HomCoord
 * \brief Homogeneous coordinate vector
 */
class HomCoord 
{
private:

    //! coordinate data
#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER < 1310)
  // Hack Intel compiler 12 issues with -O2 optimization
  int homCoord[5];
#else
  int homCoord[4];
#endif

public:
  friend class HomXform;

  static HomCoord unitv[3];
  static HomCoord IDENTITY;

    //! constructors
  HomCoord();
  HomCoord(const int coords[], const int num_coords = 4);
  HomCoord(const int coord0, const int coord1, 
           const int coord2, const int coord3);
  HomCoord(const int coord0, const int coord1, 
           const int coord2);
  HomCoord(const HomCoord &coord);

    //! set function
  void set(const int coords[]);
  void set(const int i, const int j, const int k, 
           const int h = 1);

    //! get function
  const int *hom_coord() const {return homCoord;}

    //! parameter-based access functions
  int i() const {return homCoord[0];}
  int j() const {return homCoord[1];}
  int k() const {return homCoord[2];}
  int h() const {return homCoord[3];}

    //! squared length
  int length_squared() const;

    //! length
  int length() const;

    //! normalize
  void normalize();

    //! operators
  HomCoord &operator*=(const HomXform &rhs2);
  HomCoord operator*(const HomXform &rhs2) const;
  HomCoord &operator*=(const int mult);
  HomCoord operator*(const int mult) const;
  inline HomCoord &operator/=(const HomXform &rhs2);
  HomCoord operator/(const HomXform &rhs2) const;
  inline HomCoord &operator/=(const int mult);
  HomCoord operator/(const int mult) const;
  inline HomCoord &operator+=(const HomCoord &rhs1);
  HomCoord operator+(const HomCoord &rhs1) const;
  inline HomCoord &operator-=(const HomCoord &rhs1);
  HomCoord operator-(const HomCoord &rhs1) const;
  HomCoord &operator=(const HomCoord &rhs);

    // dot product
  int operator%(const HomCoord &rhs) const;
  
    // cross product
  HomCoord operator*(const HomCoord &rhs) const;
  HomCoord &operator*=(const HomCoord &rhs);
  
  bool operator==(const HomCoord &rhs1) const;
  bool operator!=(const HomCoord &rhs1) const;
  bool operator>=(const HomCoord &rhs1) const;
  bool operator<=(const HomCoord &rhs1) const;
  bool operator>(const HomCoord &rhs1) const;
  bool operator<(const HomCoord &rhs1) const;
  inline int operator[](const int &param) const;
  int &operator[](const int &param);

};
  
/** \class HomXform
 * \brief Homogeneous coordinate transformation matrix
 */
class HomXform
{

private:

    //! the matrix; don't bother storing the last column, since we assume for now it's
    //! always unused
  int xForm[16];

public:

  friend class HomCoord;

  static HomXform IDENTITY;

    //! constructor from matrix
  HomXform(const int matrix[16]);

    //! bare constructor
  HomXform();
  
    //! constructor from rotation, scaling, translation
  HomXform(const int rotate[9], const int scale[3], const int translate[3]);

    //! constructor taking 16 ints, useful for efficient operators
  HomXform(int i1, int i2, int i3, int i4,
           int i5, int i6, int i7, int i8,
           int i9, int i10, int i11, int i12,
           int i13, int i14, int i15, int i16);

    //! return this.inverse
  inline HomXform inverse() const;

    //! compute a transform from three points
  void three_pt_xform(const HomCoord &p1, const HomCoord &q1,
                      const HomCoord &p2, const HomCoord &q2,
                      const HomCoord &p3, const HomCoord &q3);

    //! operators
  int operator[](const int &count) const;
  int &operator[](const int &count);
  bool operator==(const HomXform &rhs) const;
  bool operator!=(const HomXform &rhs) const;

  HomXform &operator=(const HomXform &rhs);
  HomXform &operator*=(const HomXform &rhs);
  HomXform operator*(const HomXform &rhs2) const;
};

inline HomCoord::HomCoord() 
{
  homCoord[0] = 0;
  homCoord[1] = 0;
  homCoord[2] = 0;
  homCoord[3] = 0;
}
  
inline HomCoord::HomCoord(const int coords[], const int num_coords) 
{
  for (int tmpj = 0; tmpj < num_coords; tmpj++) homCoord[tmpj] = coords[tmpj];
  if (num_coords != 4) homCoord[3] = 1;
}

inline HomCoord::HomCoord(const int coord0, const int coord1, 
                          const int coord2, const int coord3) 
{
  homCoord[0] = coord0;
  homCoord[1] = coord1;
  homCoord[2] = coord2;
  homCoord[3] = coord3;
}
  
inline HomCoord::HomCoord(const int coord0, const int coord1, 
                          const int coord2) 
{
  homCoord[0] = coord0;
  homCoord[1] = coord1;
  homCoord[2] = coord2;
  homCoord[3] = 1;
}

inline HomCoord::HomCoord(const HomCoord &coords) 
{
#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER < 1310)
  // Hack Intel compiler 12 issues with -O2 optimization
  int coord0 = coords[0];
  int coord1 = coords[1];
  int coord2 = coords[2];
  int coord3 = coords[3];
  homCoord[0] = coord0;
  homCoord[1] = coord1;
  homCoord[2] = coord2;
  homCoord[3] = coord3;
#else
  homCoord[0] = coords[0];
  homCoord[1] = coords[1];
  homCoord[2] = coords[2];
  homCoord[3] = coords[3];
#endif
}
  
inline void HomCoord::set(const int coords[]) 
{
    homCoord[0] = coords[0];
    homCoord[1] = coords[1];
    homCoord[2] = coords[2];
    homCoord[3] = coords[3];
}
  
inline void HomCoord::set(const int ip, const int jp, const int kp,
                          const int hp)
{
    homCoord[0] = ip;
    homCoord[1] = jp;
    homCoord[2] = kp;
    homCoord[3] = hp;
}
  
inline HomCoord &HomCoord::operator=(const HomCoord &rhs1)
{
  homCoord[0] = rhs1.homCoord[0];
  homCoord[1] = rhs1.homCoord[1];
  homCoord[2] = rhs1.homCoord[2];
  homCoord[3] = rhs1.homCoord[3];
  return *this;
}

    //! squared length
inline int HomCoord::length_squared() const 
{
  return homCoord[0]*homCoord[0] +
    homCoord[1]*homCoord[1] +
    homCoord[2]*homCoord[2];
}


    //! length
inline int HomCoord::length() const 
{
  return (int) sqrt((float)length_squared());
}

    //! normalize
inline void HomCoord::normalize()
{
  *this /= length();
}

    // dot product
inline int HomCoord::operator%(const HomCoord &rhs) const 
{
  return homCoord[0]*rhs.homCoord[0] +
    homCoord[1]*rhs.homCoord[1] +
    homCoord[2]*rhs.homCoord[2];
}

    // cross product
inline HomCoord HomCoord::operator*(const HomCoord &rhs) const 
{
  return HomCoord(
    homCoord[1]*rhs.homCoord[2] - homCoord[2]*rhs.homCoord[1],
    homCoord[2]*rhs.homCoord[0] - homCoord[0]*rhs.homCoord[2],
    homCoord[0]*rhs.homCoord[1] - homCoord[1]*rhs.homCoord[0]);
}

inline HomCoord &HomCoord::operator*=(const HomCoord &rhs) 
{
  *this = HomCoord(
    homCoord[1]*rhs.homCoord[2] - homCoord[2]*rhs.homCoord[1],
    homCoord[2]*rhs.homCoord[0] - homCoord[0]*rhs.homCoord[2],
    homCoord[0]*rhs.homCoord[1] - homCoord[1]*rhs.homCoord[0]);
  
  return *this;
}
  
inline bool HomCoord::operator==(const HomCoord &rhs1) const 
{
  return (homCoord[0] == rhs1.homCoord[0] &&
          homCoord[1] == rhs1.homCoord[1] &&
          homCoord[2] == rhs1.homCoord[2] &&
          homCoord[3] == rhs1.homCoord[3]);
}

inline bool HomCoord::operator!=(const HomCoord &rhs1) const 
{
  return (homCoord[0] != rhs1.homCoord[0] ||
          homCoord[1] != rhs1.homCoord[1] ||
          homCoord[2] != rhs1.homCoord[2] ||
          homCoord[3] != rhs1.homCoord[3]);
}

inline bool HomCoord::operator>=(const HomCoord &rhs1) const 
{
  return (homCoord[0] >= rhs1.homCoord[0] &&
          homCoord[1] >= rhs1.homCoord[1] &&
          homCoord[2] >= rhs1.homCoord[2] &&
          homCoord[3] == rhs1.homCoord[3]);
}

inline bool HomCoord::operator<=(const HomCoord &rhs1) const 
{
  return (homCoord[0] <= rhs1.homCoord[0] &&
          homCoord[1] <= rhs1.homCoord[1] &&
          homCoord[2] <= rhs1.homCoord[2] &&
          homCoord[3] == rhs1.homCoord[3]);
}

inline bool HomCoord::operator<(const HomCoord &rhs1) const 
{
  return (homCoord[0] < rhs1.homCoord[0] &&
          homCoord[1] < rhs1.homCoord[1] &&
          homCoord[2] < rhs1.homCoord[2] &&
          homCoord[3] == rhs1.homCoord[3]);
}

inline bool HomCoord::operator>(const HomCoord &rhs1) const 
{
  return (homCoord[0] > rhs1.homCoord[0] &&
          homCoord[1] > rhs1.homCoord[1] &&
          homCoord[2] > rhs1.homCoord[2] &&
          homCoord[3] == rhs1.homCoord[3]);
}

inline HomCoord HomCoord::operator*(const HomXform &rhs2) const 
{
  return HomCoord(
//    homCoord[0]*rhs2[4*0+0] + homCoord[1]*rhs2[4*1+0] + 
//    homCoord[2]*rhs2[4*2+0] + homCoord[3]*rhs2[4*3+0],
    homCoord[0]*rhs2.xForm[0] + homCoord[1]*rhs2.xForm[4] + 
    homCoord[2]*rhs2.xForm[8] + homCoord[3]*rhs2.xForm[12],

//    homCoord[0]*rhs2.xForm[4*0+1] + homCoord[1]*rhs2.xForm[4*1+1] + 
//    homCoord[2]*rhs2.xForm[4*2+1] + homCoord[3]*rhs2.xForm[4*3+1],
    homCoord[0]*rhs2.xForm[1] + homCoord[1]*rhs2.xForm[5] + 
    homCoord[2]*rhs2.xForm[9] + homCoord[3]*rhs2.xForm[13],

//    homCoord[0]*rhs2.xForm[4*0+2] + homCoord[1]*rhs2.xForm[4*1+2] + 
//    homCoord[2]*rhs2.xForm[4*2+2] + homCoord[3]*rhs2.xForm[4*3+2],
    homCoord[0]*rhs2.xForm[2] + homCoord[1]*rhs2.xForm[6] + 
    homCoord[2]*rhs2.xForm[10] + homCoord[3]*rhs2.xForm[14],

//    homCoord[0]*rhs2.xForm[4*0+3] + homCoord[1]*rhs2.xForm[4*1+3] + 
//    homCoord[2]*rhs2.xForm[4*2+3] + homCoord[3]*rhs2.xForm[4*3+3]
    homCoord[0]*rhs2.xForm[3] + homCoord[1]*rhs2.xForm[7] + 
    homCoord[2]*rhs2.xForm[11] + homCoord[3]*rhs2.xForm[15]
    );
}

inline HomCoord &HomCoord::operator*=(const HomXform &rhs2)
{
  *this = HomCoord(
//    homCoord[0]*rhs2.xForm[4*0+0] + homCoord[1]*rhs2.xForm[4*1+0] + 
//    homCoord[2]*rhs2.xForm[4*2+0] + homCoord[3]*rhs2.xForm[4*3+0],
    homCoord[0]*rhs2.xForm[0] + homCoord[1]*rhs2.xForm[4] + 
    homCoord[2]*rhs2.xForm[8] + homCoord[3]*rhs2.xForm[12],

//    homCoord[0]*rhs2.xForm[4*0+1] + homCoord[1]*rhs2.xForm[4*1+1] + 
//    homCoord[2]*rhs2.xForm[4*2+1] + homCoord[3]*rhs2.xForm[4*3+1],
    homCoord[0]*rhs2.xForm[1] + homCoord[1]*rhs2.xForm[5] + 
    homCoord[2]*rhs2.xForm[9] + homCoord[3]*rhs2.xForm[13],

//    homCoord[0]*rhs2.xForm[4*0+2] + homCoord[1]*rhs2.xForm[4*1+2] + 
//    homCoord[2]*rhs2.xForm[4*2+2] + homCoord[3]*rhs2.xForm[4*3+2],
    homCoord[0]*rhs2.xForm[2] + homCoord[1]*rhs2.xForm[6] + 
    homCoord[2]*rhs2.xForm[10] + homCoord[3]*rhs2.xForm[14],

//    homCoord[0]*rhs2.xForm[4*0+3] + homCoord[1]*rhs2.xForm[4*1+3] + 
//    homCoord[2]*rhs2.xForm[4*2+3] + homCoord[3]*rhs2.xForm[4*3+3]
    homCoord[0]*rhs2.xForm[3] + homCoord[1]*rhs2.xForm[7] + 
    homCoord[2]*rhs2.xForm[11] + homCoord[3]*rhs2.xForm[15]
    );
  return *this;
}

inline HomCoord HomCoord::operator*(const int mult) const 
{
  return HomCoord(mult*homCoord[0], mult*homCoord[1], mult*homCoord[2]);
}

inline HomCoord &HomCoord::operator*=(const int mult)
{
  homCoord[0] *= mult;
  homCoord[1] *= mult;
  homCoord[2] *= mult;
  return *this;
}

inline HomCoord HomCoord::operator/(const int div) const 
{
  return HomCoord(homCoord[0]/div, homCoord[1]/div, homCoord[2]/div);
}

inline HomCoord &HomCoord::operator/=(const int div)
{
  homCoord[0] /= div;
  homCoord[1] /= div;
  homCoord[2] /= div;
  return *this;
}

inline HomCoord HomCoord::operator-(const HomCoord &rhs2) const 
{
  return HomCoord(*this) -= rhs2;
}

inline HomCoord &HomCoord::operator-=(const HomCoord &rhs2)
{
  homCoord[0] -= rhs2[0];
  homCoord[1] -= rhs2[1];
  homCoord[2] -= rhs2[2];
  return *this;
}

inline HomCoord HomCoord::operator+(const HomCoord &rhs2) const 
{
  return HomCoord(*this) += rhs2;
}

inline HomCoord &HomCoord::operator+=(const HomCoord &rhs2)
{
  homCoord[0] += rhs2[0];
  homCoord[1] += rhs2[1];
  homCoord[2] += rhs2[2];
  return *this;
}

inline HomCoord HomCoord::operator/(const HomXform &rhs2) const 
{
  return HomCoord(*this) /= rhs2;
}

inline HomCoord &HomCoord::operator/=(const HomXform &rhs2)
{
  HomXform inv = rhs2.inverse();
  *this *= inv;
  return *this;
}

inline int HomCoord::operator[](const int &param) const
{
  return homCoord[param];
}

inline int &HomCoord::operator[](const int &param)
{
  return homCoord[param];
}

inline std::ostream &operator<<(std::ostream &str, const HomCoord &hc)
{
  str << "(" << hc.i() << "," << hc.j() << "," << hc.k() << ")";
  return str;
}

inline HomXform::HomXform(const int matrix[16]) 
{
  for (int i = 0; i < 16; i++)
    xForm[i] = matrix[i];
}

inline HomXform::HomXform() 
{
}

inline HomXform::HomXform(const int rotate[9], const int scale[3], 
                          const int translate[3]) 
{
  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++)
      xForm[i*4+j] = rotate[i*3+j]*scale[j];

    xForm[12+i] = translate[i];
  }
  xForm[3] = 0;
  xForm[7] = 0;
  xForm[11] = 0;
  xForm[15] = 1;
  
}

inline HomXform::HomXform(int i1, int i2, int i3, int i4,
                          int i5, int i6, int i7, int i8,
                          int i9, int i10, int i11, int i12,
                          int i13, int i14, int i15, int i16) 
{
  xForm[0] = i1; xForm[1] = i2; xForm[2] = i3; xForm[3] = i4;
  xForm[4] = i5; xForm[5] = i6; xForm[6] = i7; xForm[7] = i8;
  xForm[8] = i9; xForm[9] = i10; xForm[10] = i11; xForm[11] = i12;
  xForm[12] = i13; xForm[13] = i14; xForm[14] = i15; xForm[15] = i16;
}

inline HomXform &HomXform::operator=(const HomXform &rhs) 
{
  for (int i = 0; i < 16; i++)
    xForm[i] = rhs.xForm[i];

  return *this;
}

inline HomXform HomXform::operator*(const HomXform &rhs2) const
{
  return HomXform(
//  temp.XFORM(0,0)
    XFORM(0,0)*rhs2.XFORM(0,0) + XFORM(0,1)*rhs2.XFORM(1,0) + 
    XFORM(0,2)*rhs2.XFORM(2,0) + XFORM(0,3)*rhs2.XFORM(3,0),
//  temp.XFORM(0,1)  
    XFORM(0,0)*rhs2.XFORM(0,1) + XFORM(0,1)*rhs2.XFORM(1,1) + 
    XFORM(0,2)*rhs2.XFORM(2,1) + XFORM(0,3)*rhs2.XFORM(3,1),
//  temp.XFORM(0,2)  
    XFORM(0,0)*rhs2.XFORM(0,2) + XFORM(0,1)*rhs2.XFORM(1,2) + 
    XFORM(0,2)*rhs2.XFORM(2,2) + XFORM(0,3)*rhs2.XFORM(3,2),
//  temp.XFORM(0,3)  
    XFORM(0,0)*rhs2.XFORM(0,3) + XFORM(0,1)*rhs2.XFORM(1,3) + 
    XFORM(0,2)*rhs2.XFORM(2,3) + XFORM(0,3)*rhs2.XFORM(3,3),

//  temp.XFORM(1,0)  
    XFORM(1,0)*rhs2.XFORM(0,0) + XFORM(1,1)*rhs2.XFORM(1,0) + 
    XFORM(1,2)*rhs2.XFORM(2,0) + XFORM(1,3)*rhs2.XFORM(3,0),
//  temp.XFORM(1,1)  
    XFORM(1,0)*rhs2.XFORM(0,1) + XFORM(1,1)*rhs2.XFORM(1,1) + 
    XFORM(1,2)*rhs2.XFORM(2,1) + XFORM(1,3)*rhs2.XFORM(3,1),
//  temp.XFORM(1,2)  
    XFORM(1,0)*rhs2.XFORM(0,2) + XFORM(1,1)*rhs2.XFORM(1,2) + 
    XFORM(1,2)*rhs2.XFORM(2,2) + XFORM(1,3)*rhs2.XFORM(3,2),
//  temp.XFORM(1,3)  
    XFORM(1,0)*rhs2.XFORM(0,3) + XFORM(1,1)*rhs2.XFORM(1,3) + 
    XFORM(1,2)*rhs2.XFORM(2,3) + XFORM(1,3)*rhs2.XFORM(3,3),

//  temp.XFORM(2,0)  
    XFORM(2,0)*rhs2.XFORM(0,0) + XFORM(2,1)*rhs2.XFORM(1,0) + 
    XFORM(2,2)*rhs2.XFORM(2,0) + XFORM(2,3)*rhs2.XFORM(3,0),
//  temp.XFORM(2,1)  
    XFORM(2,0)*rhs2.XFORM(0,1) + XFORM(2,1)*rhs2.XFORM(1,1) + 
    XFORM(2,2)*rhs2.XFORM(2,1) + XFORM(2,3)*rhs2.XFORM(3,1),
//  temp.XFORM(2,2)  
    XFORM(2,0)*rhs2.XFORM(0,2) + XFORM(2,1)*rhs2.XFORM(1,2) + 
    XFORM(2,2)*rhs2.XFORM(2,2) + XFORM(2,3)*rhs2.XFORM(3,2),
//  temp.XFORM(2,3)  
    XFORM(2,0)*rhs2.XFORM(0,3) + XFORM(2,1)*rhs2.XFORM(1,3) + 
    XFORM(2,2)*rhs2.XFORM(2,3) + XFORM(2,3)*rhs2.XFORM(3,3),

//  temp.XFORM(3,0)  
//  xForm[12]*rhs2.xForm[0] + xForm[13]*rhs2.xForm[4] + xForm[14]*rhs2.xForm[8] + xForm[15]*rhs2.xForm[12]
    XFORM(3,0)*rhs2.XFORM(0,0) + XFORM(3,1)*rhs2.XFORM(1,0) + 
    XFORM(3,2)*rhs2.XFORM(2,0) + XFORM(3,3)*rhs2.XFORM(3,0),
//  temp.XFORM(3,1)  
//  xForm[12]*rhs2.xForm[1] + xForm[13]*rhs2.xForm[5] + xForm[14]*rhs2.xForm[9] + xForm[15]*rhs2.xForm[13]
    XFORM(3,0)*rhs2.XFORM(0,1) + XFORM(3,1)*rhs2.XFORM(1,1) + 
    XFORM(3,2)*rhs2.XFORM(2,1) + XFORM(3,3)*rhs2.XFORM(3,1),
//  temp.XFORM(3,2)  
//  xForm[12]*rhs2.xForm[2] + xForm[13]*rhs2.xForm[6] + xForm[14]*rhs2.xForm[10] + xForm[15]*rhs2.xForm[14]
    XFORM(3,0)*rhs2.XFORM(0,2) + XFORM(3,1)*rhs2.XFORM(1,2) + 
    XFORM(3,2)*rhs2.XFORM(2,2) + XFORM(3,3)*rhs2.XFORM(3,2),
//  temp.XFORM(3,3)  
//  xForm[12]*rhs2.xForm[3] + xForm[13]*rhs2.xForm[7] + xForm[14]*rhs2.xForm[11] + xForm[15]*rhs2.xForm[15]
    XFORM(3,0)*rhs2.XFORM(0,3) + XFORM(3,1)*rhs2.XFORM(1,3) + 
    XFORM(3,2)*rhs2.XFORM(2,3) + XFORM(3,3)*rhs2.XFORM(3,3));
}

inline HomXform &HomXform::operator*=(const HomXform &rhs2)
{
  *this =  HomXform(
//  temp.XFORM(0,0)
    XFORM(0,0)*rhs2.XFORM(0,0) + XFORM(0,1)*rhs2.XFORM(1,0) + 
    XFORM(0,2)*rhs2.XFORM(2,0) + XFORM(0,3)*rhs2.XFORM(3,0),
//  temp.XFORM(0,1)  
    XFORM(0,0)*rhs2.XFORM(0,1) + XFORM(0,1)*rhs2.XFORM(1,1) + 
    XFORM(0,2)*rhs2.XFORM(2,1) + XFORM(0,3)*rhs2.XFORM(3,1),
//  temp.XFORM(0,2)  
    XFORM(0,0)*rhs2.XFORM(0,2) + XFORM(0,1)*rhs2.XFORM(1,2) + 
    XFORM(0,2)*rhs2.XFORM(2,2) + XFORM(0,3)*rhs2.XFORM(3,2),
//  temp.XFORM(0,3)  
    XFORM(0,0)*rhs2.XFORM(0,3) + XFORM(0,1)*rhs2.XFORM(1,3) + 
    XFORM(0,2)*rhs2.XFORM(2,3) + XFORM(0,3)*rhs2.XFORM(3,3),

//  temp.XFORM(1,0)  
    XFORM(1,0)*rhs2.XFORM(0,0) + XFORM(1,1)*rhs2.XFORM(1,0) + 
    XFORM(1,2)*rhs2.XFORM(2,0) + XFORM(1,3)*rhs2.XFORM(3,0),
//  temp.XFORM(1,1)  
    XFORM(1,0)*rhs2.XFORM(0,1) + XFORM(1,1)*rhs2.XFORM(1,1) + 
    XFORM(1,2)*rhs2.XFORM(2,1) + XFORM(1,3)*rhs2.XFORM(3,1),
//  temp.XFORM(1,2)  
    XFORM(1,0)*rhs2.XFORM(0,2) + XFORM(1,1)*rhs2.XFORM(1,2) + 
    XFORM(1,2)*rhs2.XFORM(2,2) + XFORM(1,3)*rhs2.XFORM(3,2),
//  temp.XFORM(1,3)  
    XFORM(1,0)*rhs2.XFORM(0,3) + XFORM(1,1)*rhs2.XFORM(1,3) + 
    XFORM(1,2)*rhs2.XFORM(2,3) + XFORM(1,3)*rhs2.XFORM(3,3),

//  temp.XFORM(2,0)  
    XFORM(2,0)*rhs2.XFORM(0,0) + XFORM(2,1)*rhs2.XFORM(1,0) + 
    XFORM(2,2)*rhs2.XFORM(2,0) + XFORM(2,3)*rhs2.XFORM(3,0),
//  temp.XFORM(2,1)  
    XFORM(2,0)*rhs2.XFORM(0,1) + XFORM(2,1)*rhs2.XFORM(1,1) + 
    XFORM(2,2)*rhs2.XFORM(2,1) + XFORM(2,3)*rhs2.XFORM(3,1),
//  temp.XFORM(2,2)  
    XFORM(2,0)*rhs2.XFORM(0,2) + XFORM(2,1)*rhs2.XFORM(1,2) + 
    XFORM(2,2)*rhs2.XFORM(2,2) + XFORM(2,3)*rhs2.XFORM(3,2),
//  temp.XFORM(2,3)  
    XFORM(2,0)*rhs2.XFORM(0,3) + XFORM(2,1)*rhs2.XFORM(1,3) + 
    XFORM(2,2)*rhs2.XFORM(2,3) + XFORM(2,3)*rhs2.XFORM(3,3),

//  temp.XFORM(3,0)  
    XFORM(3,0)*rhs2.XFORM(0,0) + XFORM(3,1)*rhs2.XFORM(1,0) + 
    XFORM(3,2)*rhs2.XFORM(2,0) + XFORM(3,3)*rhs2.XFORM(3,0),
//  temp.XFORM(3,1)  
    XFORM(3,0)*rhs2.XFORM(0,1) + XFORM(3,1)*rhs2.XFORM(1,1) + 
    XFORM(3,2)*rhs2.XFORM(2,1) + XFORM(3,3)*rhs2.XFORM(3,1),
//  temp.XFORM(3,2)  
    XFORM(3,0)*rhs2.XFORM(0,2) + XFORM(3,1)*rhs2.XFORM(1,2) + 
    XFORM(3,2)*rhs2.XFORM(2,2) + XFORM(3,3)*rhs2.XFORM(3,2),
//  temp.XFORM(3,3)  
    XFORM(3,0)*rhs2.XFORM(0,3) + XFORM(3,1)*rhs2.XFORM(1,3) + 
    XFORM(3,2)*rhs2.XFORM(2,3) + XFORM(3,3)*rhs2.XFORM(3,3));

  return *this;
}

inline int HomXform::operator[](const int &count) const
{
  return xForm[count];
}

inline int &HomXform::operator[](const int &count)
{
  return xForm[count];
}

inline bool HomXform::operator==(const HomXform &rhs) const 
{
  return (
    xForm[0] == rhs.xForm[0] && xForm[1] == rhs.xForm[1] &&
    xForm[2] == rhs.xForm[2] && xForm[3] == rhs.xForm[3] &&
    xForm[4] == rhs.xForm[4] && xForm[5] == rhs.xForm[5] &&
    xForm[6] == rhs.xForm[6] && xForm[7] == rhs.xForm[7] &&
    xForm[8] == rhs.xForm[8] && xForm[9] == rhs.xForm[9] &&
    xForm[10] == rhs.xForm[10] && xForm[11] == rhs.xForm[11] &&
    xForm[12] == rhs.xForm[12] && xForm[13] == rhs.xForm[13] &&
    xForm[14] == rhs.xForm[14] && xForm[15] == rhs.xForm[15]);
}

inline bool HomXform::operator!=(const HomXform &rhs) const 
{
  return (
    xForm[0] != rhs.xForm[0] || xForm[1] != rhs.xForm[1] ||
    xForm[2] != rhs.xForm[2] || xForm[3] != rhs.xForm[3] ||
    xForm[4] != rhs.xForm[4] || xForm[5] != rhs.xForm[5] ||
    xForm[6] != rhs.xForm[6] || xForm[7] != rhs.xForm[7] ||
    xForm[8] != rhs.xForm[8] || xForm[9] != rhs.xForm[9] ||
    xForm[10] != rhs.xForm[10] || xForm[11] != rhs.xForm[11] ||
    xForm[12] != rhs.xForm[12] || xForm[13] != rhs.xForm[13] ||
    xForm[14] != rhs.xForm[14] || xForm[15] != rhs.xForm[15]);
}

inline HomXform HomXform::inverse() const 
{

/*
// original code:

  HomXform tmp;
  
    // assign the diagonal
  tmp[0] = xForm[0];
  tmp[5] = xForm[5];
  tmp[10] = xForm[10];
  tmp[15] = xForm[15];
  
    // invert the rotation matrix
  tmp[XFORM_INDEX(0,1)] = XFORM(1,0);
  tmp[XFORM_INDEX(0,2)] = XFORM(2,0);
  tmp[XFORM_INDEX(1,0)] = XFORM(0,1);
  tmp[XFORM_INDEX(1,2)] = XFORM(2,1);
  tmp[XFORM_INDEX(2,0)] = XFORM(0,2);
  tmp[XFORM_INDEX(2,1)] = XFORM(1,2);

    // negative translate * Rinv
  tmp[XFORM_INDEX(3,0)] = -(XFORM(3,0)*tmp.XFORM(0,0) + XFORM(3,1)*tmp.XFORM(1,0) + XFORM(3,2)*tmp.XFORM(2,0));
  tmp[XFORM_INDEX(3,1)] = -(XFORM(3,0)*tmp.XFORM(0,1) + XFORM(3,1)*tmp.XFORM(1,1) + XFORM(3,2)*tmp.XFORM(2,1));
  tmp[XFORM_INDEX(3,2)] = -(XFORM(3,0)*tmp.XFORM(0,2) + XFORM(3,1)*tmp.XFORM(1,2) + XFORM(3,2)*tmp.XFORM(2,2));
  
    // zero last column
  tmp[XFORM_INDEX(0,3)] = 0;
  tmp[XFORM_INDEX(1,3)] = 0;
  tmp[XFORM_INDEX(2,3)] = 0;

    // h factor
  tmp[XFORM_INDEX(3,3)] = 1;

  return tmp;
*/

// more efficient, but somewhat confusing (remember, column-major):

  return HomXform(
      // row 0
    xForm[0], XFORM(1,0), XFORM(2,0), 0,
      // row 1
    XFORM(0,1), xForm[5], XFORM(2,1), 0,
      // row 2
    XFORM(0,2), XFORM(1,2), xForm[10], 0,
      // row 3
    -(XFORM(3,0)*xForm[0] + XFORM(3,1)*XFORM(0,1) + XFORM(3,2)*XFORM(0,2)),
    -(XFORM(3,0)*XFORM(1,0) + XFORM(3,1)*xForm[5] + XFORM(3,2)*XFORM(1,2)),
    -(XFORM(3,0)*XFORM(2,0) + XFORM(3,1)*XFORM(2,1) + XFORM(3,2)*xForm[10]),
    1);
}

} // namespace moab 

#endif
