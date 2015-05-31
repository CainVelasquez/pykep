#ifndef KEP_TOOLBOX_VECTORS_H
#define KEP_TOOLBOX_VECTORS_H

#include <vector>
#include "../core_functions/array3D_operations.h"

namespace kep_toolbox
{
  void scale(array3D &vect, array3D &out, double scalar)
  {
    for (int i = 0; i < 3; i++)
    {
      out[i] = scalar * vect[i];
    }
  };

  /* Project 2D cartesian coordinates into 3D space given x and z 3D axes and x', y' */
  void rel_coord(array3D &out, array3D z_axis, array3D x_axis, double x_prime, double y_prime)
  {
    /* Unit radial and tangential components */
    array3D ur;
    array3D ut;

    /* Construct unit components */
    vers(ur, x_axis);
    cross(ut, z_axis, x_axis);
    vers(ut, ut);

    /* Construct final vector, r = x' * <X> + y' * <Y> */
    scale(ur, ur, x_prime);
    scale(ut, ut, y_prime);
    sum(out, ur, ut);
  };
}

#endif
