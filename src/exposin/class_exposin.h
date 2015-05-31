#ifndef KEP_TOOLBOX_CLASS_EXPOSIN_H
#define KEP_TOOLBOX_CLASS_EXPOSIN_H

#include <cmath>
#include "exposin.h"
#include "../astro_constants.h"

namespace kep_toolbox
{
  class class_exposin
  {
  private:
    double k2, r1_m, r2_m, psi;
    int multi_revs;
  public:
    class_exposin(){};
    class_exposin(double k2, double r1_m, double r2_m, double angle, int multi_revs)
    {
      this->k2 = k2;
      this->r1_m = r1_m;
      this->r2_m = r2_m;
      this->psi = 2.0 * M_PI * multi_revs + angle;
      this->multi_revs = multi_revs;
    };
    double get_psi() const
    {
      return psi;
    };
    int tany1_range(double (&min_max_mid)[3]) const
    {
      double logr1r2 = log(r1_m / r2_m);
      double cosk2O = cos(k2 * psi);
      double delta = 2.0 * (1.0 - cosk2O) / pow(k2, 4.0) - logr1r2*logr1r2;
      if (delta < 0.0)
      {
        return 1;
      }
      double tany1min = k2/2.0 * (-logr1r2 / tan(k2*psi/2) - sqrt(delta));
      double tany1max = k2/2.0 * (-logr1r2 / tan(k2*psi/2) + sqrt(delta));
      min_max_mid[0] = tany1min;
      min_max_mid[1] = tany1max;
      min_max_mid[2] = -k2 / 2.0 * logr1r2 / tan(k2*psi/2);

      return 0;
    };
    int create_exposin(double tany1, exposin &expsn) const
    {
      double range[3]{};
      tany1_range(range);
      if (tany1 < range[0] || tany1 > range[1])
      {
        return 1;
      }
      double logr1r2 = log(r1_m / r2_m);
      double sink2O = sin(k2 * psi);
      double cosk2O = cos(k2 * psi);

      double k1_sqr = pow((logr1r2 + tany1 / k2 * sink2O)/(1.0 - cosk2O), 2.0) + pow(tany1 / k2, 2.0);
      double k1_sign = (logr1r2 + tany1 / k2 * sink2O)/(1.0 - cosk2O);

      double k1;
      if (k1_sign < 0)
      {
        k1 = -sqrt(k1_sqr);
      }
      else
      {
        k1 = sqrt(k1_sqr);
      }
      double phi = acos(tany1/k1/k2);
      double k0 = r1_m/exp(k1*sin(phi));
      expsn.set(k0, k1, k2, phi);
      return 0;
    };
    int search_tany1(double dT, double mu, double &tany1) const
    {
      const double STOP_CRITERION = 10000.0; // seconds
      double MAX_ITERS = 10000;
      const double MAX_DT = 1000 * 365.0 * ASTRO_DAY2SEC;

      double range[3];
      if (tany1_range(range))
      {
        return 2;
      }
      exposin exps;

      double tany1_a = range[0];
      double tany1_b = range[1];

      double tany1_c, tof_c;

      do
      {
        if(create_exposin(tany1_a, exps))
        {
          return 3;
        }
        double tof_a = exps.dT(psi, mu);

        if(create_exposin(tany1_b, exps))
        {
          return 4;
        }
        double tof_b = exps.dT(psi, mu);
        tof_b = fmin(tof_b, MAX_DT);

        tany1_c = tany1_b - (tof_b - dT) * (tany1_b - tany1_a) / (tof_b - tof_a);
        create_exposin(tany1_c, exps);
        tof_c = exps.dT(psi, mu);

        if ((tof_a > 0.0 && tof_c > 0.0) || (tof_a < 0.0 && tof_c < 0.0))
        {
          tany1_a = tany1_c;
        }
        else
        {
          tany1_b = tany1_c;
        }
      }
      while(fabs(tof_c - dT) > STOP_CRITERION && (--MAX_ITERS));
      if (MAX_ITERS <= 0)
      {
        return 1;
      }
      tany1 = tany1_c;
      return 0;
    };
  };
}

#endif
