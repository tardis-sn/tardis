#ifndef BOUND_BOX_HPP
#define BOUND_BOX_HPP

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"

#include <cfloat>

namespace moab {

    class BoundBox {
  public:
      BoundBox() : bMin(DBL_MAX), bMax(-DBL_MAX) {}
      BoundBox(const CartVect &min, const CartVect &max) : 
              bMin(min), bMax(max) {}
      BoundBox(const double *corners);
      // constructor used in element maps
      BoundBox(std::vector<CartVect> points): bMin(DBL_MAX), bMax(-DBL_MAX)
      {
        for (size_t i=0; i<points.size(); i++)
        {
          update_min( points[i].array() );
          update_max( points[i].array() );
        }
      }
      ~BoundBox() {}

      bool contains_point(const double *point, const double tol = 0.0) const;
      bool intersects_box(const BoundBox &b, const double tol = 0.0) const;
      void compute_center(CartVect &center);
      void update(const BoundBox &other_box);
      void update(const double *coords);
      ErrorCode update(Interface &iface, const Range& elems);
      ErrorCode update(Interface &iface, const EntityHandle ent);
      void update_min(const BoundBox &other_box);
      void update_min(const double *coords);
      void update_max(const BoundBox &other_box);
      void update_max(const double *coords);
      ErrorCode get(double *coords);

        /** \brief Return the diagonal length of this box
         */
      double diagonal_length() const;
      
        /** \brief Return the square of the diagonal length of this box
         */
      double diagonal_squared() const;
      
        /** \brief Return square of distance from box, or zero if inside 
         * \param from_point Point from which you want distance_sq
         */
      double distance_squared(const double *from_point) const;

        /** \brief Return distance from box, or zero if inside 
         * \param from_point Point from which you want distance
         */
      double distance(const double *from_point) const;
      
      BoundBox &operator=(const BoundBox &from) {
        bMin = from.bMin;
        bMax = from.bMax;
        return *this;
      }
      inline bool operator==(const BoundBox &box) const {
        return (bMin == box.bMin && bMax == box.bMax);
      }

      CartVect bMin, bMax;
    };
    
    inline BoundBox::BoundBox(const double *corners) 
    {
        // relies on CartVect being Plain Old Data, no virtual table
      double *arr = bMin.array();
      for (int i = 0; i < 6; i++)
        arr[i] = corners[i];
    }

    inline bool BoundBox::contains_point(const double *point, const double tol) const {
      if (point[0] < bMin[0]-tol || point[0] > bMax[0]+tol ||
          point[1] < bMin[1]-tol || point[1] > bMax[1]+tol ||
          point[2] < bMin[2]-tol || point[2] > bMax[2]+tol)
        return false;
      else return true;
    }

    inline bool BoundBox::intersects_box(const BoundBox &b, const double tol) const {
      if (b.bMax[0] < bMin[0]-tol || b.bMin[0] > bMax[0]+tol ||
          b.bMax[1] < bMin[1]-tol || b.bMin[1] > bMax[1]+tol ||
          b.bMax[2] < bMin[2]-tol || b.bMin[2] > bMax[2]+tol) 
        return false;

      else return true;
    }

    inline void BoundBox::update(const BoundBox &other_box) 
    {
      update_min(other_box);
      update_max(other_box);
    }

    inline void BoundBox::update(const double *coords) 
    {
      update_min(coords);
      update_max(coords+3);
    }

    inline void BoundBox::update_min(const BoundBox &other_box) 
    {
      bMin[0] = std::min(bMin[0], other_box.bMin[0]);
      bMin[1] = std::min(bMin[1], other_box.bMin[1]);
      bMin[2] = std::min(bMin[2], other_box.bMin[2]);
    }
      
    inline void BoundBox::update_min(const double *coords) 
    {
      bMin[0] = std::min(bMin[0], coords[0]);
      bMin[1] = std::min(bMin[1], coords[1]);
      bMin[2] = std::min(bMin[2], coords[2]);
    }
      
    inline void BoundBox::update_max(const BoundBox &other_box)
    {
      bMax[0] = std::max(bMax[0], other_box.bMax[0]);
      bMax[1] = std::max(bMax[1], other_box.bMax[1]);
      bMax[2] = std::max(bMax[2], other_box.bMax[2]);
    }

    inline void BoundBox::update_max(const double *coords)
    {
      bMax[0] = std::max(bMax[0], coords[0]);
      bMax[1] = std::max(bMax[1], coords[1]);
      bMax[2] = std::max(bMax[2], coords[2]);
    }

    inline ErrorCode BoundBox::get(double *coords)
    {
      bMin.get(coords);
      bMax.get(coords+3);
      return MB_SUCCESS;
    }

    inline void BoundBox::compute_center(CartVect &center){
      center = 0.5 * (bMin + bMax);
    }

    inline std::ostream &operator<<(std::ostream& out, const BoundBox &box) {
      out << (std::string) "Min: ";
      out << box.bMin;
      out << (std::string) ", Max: ";
      out << box.bMax;
      return out;
    }

    inline ErrorCode BoundBox::update(Interface &iface, const EntityHandle ent) 
    {
      Range tmp_range(ent, ent);
      return update(iface, tmp_range);
    }

    inline double BoundBox::distance_squared(const double *from_point) const
    {
      double dist_sq = 0.0;
      for (int i = 0; i < 3; ++i) {
        if (from_point[i] < bMin[i])
          dist_sq += (bMin[i] - from_point[i]) * (bMin[i] - from_point[i]);
        else if (from_point[i] > bMax[i]) 
          dist_sq += (bMax[i] - from_point[i]) * (bMax[i] - from_point[i]);
      }
      return dist_sq;
    }

    inline double BoundBox::distance(const double *from_point) const
    {
      double dist_sq = distance_squared(from_point);
      return sqrt(dist_sq);
    }

    inline double BoundBox::diagonal_length() const 
    {
      if (DBL_MAX == bMax[0] || DBL_MAX == bMax[1] || DBL_MAX == bMax[2] ||
          DBL_MAX == bMin[0] || DBL_MAX == bMin[1] || DBL_MAX == bMin[2]) return DBL_MAX;
      return (bMax - bMin).length();
    }

    inline double BoundBox::diagonal_squared() const 
    {
      if (DBL_MAX == bMax[0] || DBL_MAX == bMax[1] || DBL_MAX == bMax[2] ||
          DBL_MAX == bMin[0] || DBL_MAX == bMin[1] || DBL_MAX == bMin[2]) return DBL_MAX;
      return (bMax - bMin).length_squared();
    }
}

#endif
