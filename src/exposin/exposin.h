/*****************************************************************************
 *   Copyright (C) 2004-2015 The PyKEP development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef KEP_TOOLBOX_EXPOSIN_H
#define KEP_TOOLBOX_EXPOSIN_H

#include <cmath>
#include "vectors.h"
#include "../serialization.h"
#include "../config.h"

#define PROPAGATION_STEPS 100

typedef boost::array<double, 9> array9D;

namespace kep_toolbox {
    /// Exponential Sinusoid Trajectory Model
    /**
     * The exposin class represents a single trajectory in the form of an exponential sinusoid curve in 3D space.
     *
     * @author Chris Andre (chrisandre01 _AT_ gmail.com)
     */
    class __KEP_TOOL_VISIBLE exposin {
    private:
        /// Shape parameters of exponential sinusoid
        double k0, k1, k2, phi;

        /// Geometry for 3D interpretation as a trajectory
        array3D r1, r2;
        double angle, psi;
        bool lw;
        int revs;
    public:

        exposin(const double &k0 = 0, const double &k1 = 0, const double &k2 = 0, const double &phi = 0) {
            set(k0, k1, k2, phi);
        }

        /// Set shape parameters
        void set(const double &k0, const double &k1, const double &k2, const double &phi) {
            this->k0 = k0;
            this->k1 = k1;
            this->k2 = k2;
            this->phi = phi;
        }

        /// Extra data for 3D projection
        void projection(const array3D &r1, const array3D &r2, const double &angle, const bool &lw, const double &psi,
                        const int &revs) {
            this->r1 = r1;
            this->r2 = r2;
            this->angle = angle;
            this->lw = lw;
            this->psi = psi;
            this->revs = revs;
        }

        /// Calculate tangent of flight path angle
        double tany(const double &theta) const {
            return k1 * k2 * cos(k2 * theta + phi);
        }

        /// Calculate radial distance magnitude
        double r_m(const double &theta) const {
            return k0 * exp(k1 * sin(k2 * theta + phi));
        }

        /// Calculate rate of change of radial distance magnitude w.r.t. time
        double r_m_dot(const double &theta, const double &mu) const {
            return r_m(theta) * k1 * cos(k2 * theta + phi) * k2 * theta_dot(theta, mu);
        }

        /// Calculate rate of change of theta w.r.t. time
        double theta_dot(const double &theta, const double &mu) const {
            return sqrt(
                    mu / pow(r_m(theta), 3.0) / (pow(tany(theta), 2.0) + k1 * k2 * k2 * sin(k2 * theta + phi) + 1));
        }

        /// Calculate time of flight from theta=0 to theta=psi
        double tof(const double &psi, const double &mu, const int &abscissas = 10) const {
            double d_theta = psi / abscissas;
            double tof_quadrature = 0.0;
            for (int i = 0; i < abscissas; i++) {
                tof_quadrature += d_theta / theta_dot(d_theta * i + d_theta / 2.0, mu);
            }
            return tof_quadrature;
        }

        /// Calculate magnitude of velocity
        double v_m(const double &theta, const double &mu) const {
            double tangential = theta_dot(theta, mu) * r_m(theta);
            double radial = r_m_dot(theta, mu);
            return sqrt(tangential * tangential + radial * radial);
        }

        /// Calculate local (non-dimensional) acceleration magnitude
        double local_a_m(const double &theta) const {
            double tan_y = tany(theta);
            double s = sin(k2 * theta + phi);
            double cosy = 1.0 / sqrt(1.0 + tan_y * tan_y);
            double tan2yk1k22s1 = tan_y * tan_y + k1 * k2 * k2 * s + 1;
            return fabs(tan_y) / 2.0 / cosy *
                   (1 / tan2yk1k22s1 - (k2 * k2 * (1 - 2 * k1 * s)) / tan2yk1k22s1 * tan2yk1k22s1);
        }

        /// Calculate dimensional acceleration magnitude
        double a_m(const double &theta, const double &mu) const {
            return mu / pow(r_m(theta), 2.0) * local_a_m(theta);
        }

        /// Calculate position vector
        void r_vec(array3D &out, const double &theta) const {
            double yp = r_m(theta) * sin(theta);
            double xp = r_m(theta) * cos(theta);
            array3D normal;
            cross(normal, r1, r2);
            if (lw) {
                scale(normal, normal, -1.0);
            }
            rel_coord(out, normal, r1, xp, yp);
        }

        /// Calculate velocity vector
        void v_vec(array3D &out, const double &theta, const double &mu) const {
            double yp = r_m_dot(theta, mu) * sin(theta) + r_m(theta) * cos(theta) * theta_dot(theta, mu);
            double xp = r_m_dot(theta, mu) * cos(theta) - r_m(theta) * sin(theta) * theta_dot(theta, mu);
            array3D normal;
            cross(normal, r1, r2);
            if (lw) {
                scale(normal, normal, -1.0);
            }
            rel_coord(out, normal, r1, xp, yp);
        }

        /// Calculate required acceleration vector (no central body attraction)
        void a_vec(array3D &out, const double &theta, const double &mu) const {
            v_vec(out, theta, mu);
            vers(out, out);
            scale(out, out, a_m(theta, mu));
            if (tany(theta) < 0.0) {
                scale(out, out, -1.0);
            }
        }

        /// Get traversed angle
        const double &get_psi() const {
            return psi;
        }

        /// Get revolutions
        const int &get_revs() const {
            return revs;
        }

        /// Calculate final mass
        double get_final_mass(const double &mu, const double &isp, const double &m) {
            double d_theta = psi / PROPAGATION_STEPS;
            double mass_quadrature = m;
            double max_thrust = 0.0;
            for (int i = 0; i < PROPAGATION_STEPS; i++) {
                double theta = d_theta * i + d_theta / 2.0;
                double dt = d_theta / theta_dot(theta, mu);
                double accel = a_m(theta, mu);
                double thrust = accel * mass_quadrature;
                double dm = -thrust / isp / ASTRO_G0;
                double delm = dt * dm;
                mass_quadrature += delm;
                max_thrust = fmax(max_thrust, thrust);
            }
            return mass_quadrature;
        }

        /// Calculate maximum thrust required along the arc
        double get_maximum_thrust(const double &mu, const double &isp, const double &m) {
            double d_theta = psi / PROPAGATION_STEPS;
            double mass_quadrature = m;
            double max_thrust = 0.0;
            for (int i = 0; i < PROPAGATION_STEPS; i++) {
                double theta = d_theta * i + d_theta / 2.0;
                double dt = d_theta / theta_dot(theta, mu);
                double accel = a_m(theta, mu);
                double thrust = accel * mass_quadrature;
                double dm = -thrust / isp / ASTRO_G0;
                double delm = dt * dm;
                mass_quadrature += delm;
                max_thrust = fmax(max_thrust, thrust);
            }
            return max_thrust;
        }

        /*
         * Helper methods and misc.
         */

    private:
        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive &ar, const unsigned int) {
            ar &k0;
            ar &k1;
            ar &k2;
            ar &phi;
            ar &r1;
            ar &r2;
            ar &angle;
            ar &lw;
            ar &psi;
            ar &revs;
        }

        friend std::ostream &operator<<(std::ostream &s, const exposin &exps) {
            s << std::setprecision(14) << "Exponential sinusoid object:" << std::endl;
            s << "k0 = " << exps.k0 << std::endl;
            s << "k1 = " << exps.k1 << std::endl;
            s << "k2 = " << exps.k2 << std::endl;
            s << "phi = " << exps.phi << std::endl;
            return s;
        }
    };
}

#endif
