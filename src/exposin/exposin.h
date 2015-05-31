#ifndef KEP_TOOLBOX_EXPOSIN_H
#define KEP_TOOLBOX_EXPOSIN_H

#include <cmath>

namespace kep_toolbox
{
  /// Exposin
  /*
  * Represents an exponential sinusoid.
  */
  class exposin
  {
  private:
    double k0, k1, k2, phi;
  public:
    exposin()
    {
      this->k0 = 0.0;
      this->k1 = 0.0;
      this->k2 = 0.0;
      this->phi = 0.0;
    };
    exposin(double k0, double k1, double k2, double phi)
    {
      this->k0 = k0;
      this->k1 = k1;
      this->k2 = k2;
      this->phi = phi;
    };
    void set(double k0, double k1, double k2, double phi)
    {
      this->k0 = k0;
      this->k1 = k1;
      this->k2 = k2;
      this->phi = phi;
    };
    double tany(double theta) const
    {
      return k1 * k2 * cos(k2 * theta + phi);
    };
    double radius(double theta) const
    {
      return k0 * exp(k1 * sin(k2 * theta + phi));
    };
    double radius_dot(double theta, double mu) const
    {
      return radius(theta) * k1 * cos(k2 * theta + phi) * k2 * theta_dot(theta, mu);
    };
    double theta_dot(double theta, double mu) const
    {
      return sqrt(mu / pow(radius(theta),3.0) / (pow(tany(theta), 2.0) + k1 * k2 * k2 * sin(k2 * theta + phi) + 1));
    };
    double dT(double psi, double mu) const
    {
      const int ABSCISSAS = 100;
      double d_theta = psi / ABSCISSAS;
      double integral = 0.0;
      for (int i = 0; i < ABSCISSAS; i++)
      {
        integral += 1.0 / theta_dot(d_theta * i + d_theta / 2.0, mu) * d_theta;
      }
      return integral;
    };
    double v_m(double theta, double mu) const
    {
      double tangential = theta_dot(theta, mu) * radius(theta);
      double radial = radius_dot(theta, mu);
      return sqrt(tangential*tangential + radial*radial);
    };
    double local_accel_m(double theta) const
    {
      double tan_y = tany(theta);
      double s = sin(k2 * theta + phi);
      double cosy = 1.0 / sqrt(1.0 + tan_y*tan_y);
      double tan2yk1k22s1 = tan_y*tan_y + k1 * k2*k2 * s + 1;
      return tan_y / 2.0 / cosy * (1 / tan2yk1k22s1 - (k2*k2 * (1 - 2 * k1 * s))/tan2yk1k22s1*tan2yk1k22s1);
    };
    double accel_m(double theta, double mu) const
    {
      return mu / pow(radius(theta), 2.0) * local_accel_m(theta);
    };
  };
}

#endif
