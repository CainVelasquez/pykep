#ifndef KEP_TOOLBOX_LAMBERT_EXPOSIN_H
#define KEP_TOOLBOX_LAMBERT_EXPOSIN_H

#include <vector>
#include <cmath>
#include "../core_functions/array3D_operations.h"
#include "class_exposin.h"
#include "exposin.h"
#include "vectors.h"

namespace kep_toolbox
{
  /// Lambert Solver for Exponential Sinusoid Trajectories
  /**
  * Represents lambert's boundary value problem for a class of trajectories
  * known as exponential sinusoids. These provide a sub-optimal but analytical
  * solution to the optimal control problem.
  *
  * @author Chris Andre (chrisandre01 _AT_ gmail.com)
  */
  class __KEP_TOOL_VISIBLE lambert_exposin
  {
  private:
    class_exposin k2_class;
    array3D r1, r2;
    int lw, multi_revs, error = 0;
    double tof, mu, tany1;
  public:
    lambert_exposin(const array3D &r1, const array3D &r2, const double tof, const double mu, const int lw, const int multi_revs, double k2 = 0.1)
    {
      this->r1 = r1;
      this->r2 = r2;
      this->tof = tof;
      this->mu = mu;
      this->lw = lw;
      this->multi_revs = multi_revs;

      double r1_m = norm(r1);
      double r2_m = norm(r2);
      double angle = acos(dot(r1, r2) / r1_m / r2_m);
      if (lw)
      {
        angle = 2.0 * M_PI - angle;
      }
      k2_class = class_exposin(k2, r1_m, r2_m, angle, multi_revs);
      error = k2_class.search_tany1(tof, mu, tany1);
    };
    int get_v1(array3D &v1) const
    {
      exposin expsn{};
      k2_class.create_exposin(tany1, expsn);
      double v1_m = expsn.v_m(0.0, mu);
      double v1_t = 1.0 / sqrt(1.0 + tany1*tany1) * v1_m; // tangential component
      double v1_r = tany1 / sqrt(1.0 + tany1*tany1) * v1_m; // radial component
      array3D normal;
      cross(normal, r1, r2);
      rel_coord(v1, normal, r1, v1_r, v1_t);
      if (lw)
      {
        scale(v1, v1, -1.0);
      }
      return error;
    };
    int get_v2(array3D &v2) const
    {
      exposin expsn{};
      k2_class.create_exposin(tany1, expsn);

      double tany1range[3]{};
      k2_class.tany1_range(tany1range);
      double tany2 = tany1range[0] + tany1range[1] - tany1;

      double v2_m = expsn.v_m(k2_class.get_psi(), mu);

      double v2_t = 1.0 / sqrt(1.0 + tany2*tany2) * v2_m;
      double v2_r = tany2 / sqrt(1.0 + tany2*tany2) * v2_m;

      array3D normal;
      cross(normal, r1, r2);

      rel_coord(v2, normal, r2, v2_r, v2_t);

      if (lw)
      {
        scale(v2, v2, -1.0);
      }
      return error;
    };
    double get_traversal_final_mass(const double isp, const double mempty, const double m) const;
    double get_max_thrust(const double isp, const double mempty, const double m) const;
  };
}

#endif
