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

#include "exposin.h"

#include <cmath>
#include "../core_functions/array3D_operations.h"

#include "../astro_constants.h"
#include "../serialization.h"
#include "../config.h"

namespace kep_toolbox {
    /// Constructor
    exposin::exposin(const double &k0, const double &k1, const double &k2, const double &phi) : m_k0(
            k0), m_k1(k1), m_k2(k2), m_phi(phi) { }

    /// Set shape parameters
    void exposin::set(const double &k0, const double &k1, const double &k2, const double &phi) {
        m_k0 = k0;
        m_k1 = k1;
        m_k2 = k2;
        m_phi = phi;
    }

    /// Extra data for 3D projection
    /** Adds extra data for 3D interpretation
     * \param[in] r1 first cartesian position
     * \param[in] r2 second cartesian position
     * \param[in] lw should long way between r1, r2 be taken in solutions
     * \param[in] revs revolutions between r1, r2
     */
    void exposin::projection(const array3D &r1, const array3D &r2, const bool &lw, const int &revs) {
        m_r1 = r1;
        m_r2 = r2;
        m_transfer_angle = angle(r1, r2, lw);
        m_lw = lw;
        m_psi = m_transfer_angle + revs * 2 * M_PI;
        m_revs = revs;
    }

    /// Calculate tangent of flight path angle
    double exposin::tany(const double &theta) const {
        return m_k1 * m_k2 * cos(m_k2 * theta + m_phi);
    }

    /// Calculate radial distance magnitude
    double exposin::r(const double &theta) const {
        return m_k0 * exp(m_k1 * sin(m_k2 * theta + m_phi));
    }

    /// Calculate rate of change of radial distance magnitude w.r.t. time
    double exposin::r_dot(const double &theta, const double &mu) const {
        return r(theta) * m_k1 * cos(m_k2 * theta + m_phi) * m_k2 * theta_dot(theta, mu);
    }

    /// Calculate rate of change of theta w.r.t. time
    double exposin::theta_dot(const double &theta, const double &mu) const {
        return sqrt(
                mu / pow(r(theta), 3.0) /
                (pow(tany(theta), 2.0) + m_k1 * m_k2 * m_k2 * sin(m_k2 * theta + m_phi) + 1.0));
    }

    /// Calculate time of flight from theta=0 to theta=psi
    double exposin::tof(const double &psi, const double &mu, const int &abscissas) const {
        // abscissas = 10 results in <1% error
        double d_theta = psi / abscissas;
        double tof_quadrature = 0.0;
        for (int i = 0; i < abscissas; i++) {
            tof_quadrature += d_theta / theta_dot(d_theta * i + d_theta / 2.0, mu);
        }
        return tof_quadrature;
    }

    /// Calculate magnitude of velocity
    double exposin::v(const double &theta, const double &mu) const {
        double tangential = theta_dot(theta, mu) * r(theta);
        double radial = r_dot(theta, mu);
        return sqrt(tangential * tangential + radial * radial);
    }

    /// Calculate local (non-dimensional) acceleration magnitude
    double exposin::local_a(const double &theta) const {
        double tan_y = tany(theta);
        double s = sin(m_k2 * theta + m_phi);
        double cosy = 1.0 / sqrt(1.0 + tan_y * tan_y);
        double tan2yk1k22s1 = tan_y * tan_y + m_k1 * m_k2 * m_k2 * s + 1.0;
        return fabs(tan_y) / 2.0 / cosy *
               (1 / tan2yk1k22s1 - (m_k2 * m_k2 * (1 - 2 * m_k1 * s)) / tan2yk1k22s1 * tan2yk1k22s1);
    }

    /// Calculate dimensional acceleration magnitude
    double exposin::a(const double &theta, const double &mu) const {
        return mu / pow(r(theta), 2.0) * local_a(theta);
    }

    /// Calculate position vector
    void exposin::r_vec(array3D &out, const double &theta) const {
        double yp = r(theta) * sin(theta);
        double xp = r(theta) * cos(theta);
        array3D normal;
        cross(normal, m_r1, m_r2);
        if (m_lw) {
            scale(normal, normal, -1.0);
        }
        rel_coord(out, normal, m_r1, xp, yp);
    }

    /// Calculate velocity vector
    void exposin::v_vec(array3D &out, const double &theta, const double &mu) const {
        double yp = r_dot(theta, mu) * sin(theta) + r(theta) * cos(theta) * theta_dot(theta, mu);
        double xp = r_dot(theta, mu) * cos(theta) - r(theta) * sin(theta) * theta_dot(theta, mu);
        array3D normal;
        cross(normal, m_r1, m_r2);
        if (m_lw) {
            scale(normal, normal, -1.0);
        }
        rel_coord(out, normal, m_r1, xp, yp);
    }

    /// Calculate required acceleration vector (no central body attraction)
    void exposin::a_vec(array3D &out, const double &theta, const double &mu) const {
        v_vec(out, theta, mu);
        vers(out, out);
        scale(out, out, a(theta, mu));
        if (tany(theta) < 0.0) {
            scale(out, out, -1.0);
        }
    }

    /// Get traversed angle
    const double &exposin::get_psi() const {
        return m_psi;
    }

    /// Get revolutions
    const int &exposin::get_revs() const {
        return m_revs;
    }

    const double &exposin::get_transfer_angle() const {
        return m_transfer_angle;
    }

    /// Calculate final mass
    double exposin::get_final_mass(const double &mu, const double &isp, const double &m) const {
        // Midpoint rule quadrature
        double d_theta = m_psi / PROPAGATION_STEPS;
        double mass_quadrature = m;
        for (int i = 0; i < PROPAGATION_STEPS; i++) {
            double theta = d_theta * i + d_theta / 2.0;
            double dt = d_theta / theta_dot(theta, mu);
            double accel = a(theta, mu);
            double thrust = accel * mass_quadrature;
            double dm = -thrust / isp / ASTRO_G0;
            double delm = dt * dm;
            mass_quadrature += delm;
        }
        return mass_quadrature;
    }

    /// Calculate maximum thrust required along the arc
    double exposin::get_maximum_thrust(const double &mu, const double &isp, const double &m) const {
        // Midpoint rule quadrature
        double d_theta = m_psi / PROPAGATION_STEPS;
        double mass_quadrature = m;
        double max_thrust = 0.0;
        for (int i = 0; i < PROPAGATION_STEPS; i++) {
            double theta = d_theta * i + d_theta / 2.0;
            double dt = d_theta / theta_dot(theta, mu);
            double accel = a(theta, mu);
            double thrust = accel * mass_quadrature;
            double dm = -thrust / isp / ASTRO_G0;
            double delm = dt * dm;
            mass_quadrature += delm;
            max_thrust = fmax(max_thrust, thrust);
        }
        return max_thrust;
    }

    /// Calculate required low-thrust delta-v
    double exposin::get_delta_v(const double &mu) const {
        // Midpoint rule quadrature
        double d_theta = m_psi / PROPAGATION_STEPS;
        double dv_quadrature = 0.0;
        for (int i = 0; i < PROPAGATION_STEPS; i++) {
            double theta = d_theta * i + d_theta / 2.0;
            double dt = d_theta / theta_dot(theta, mu);
            double accel = a(theta, mu);
            double dv = dt * accel;
            dv_quadrature += dv;
        }
        return dv_quadrature;
    }

    /// Project 2D cartesian coordinates into 3D space given x and z 3D axes and x', y'
    void exposin::rel_coord(array3D &out, const array3D &z_axis, const array3D &x_axis, const double &x_prime,
                            const double &y_prime) const {
        // Unit radial and tangential components
        array3D ur;
        array3D ut;

        // Construct unit components
        vers(ur, x_axis);
        cross(ut, z_axis, x_axis);
        vers(ut, ut);

        // Construct final vector, r = x' * <X> + y' * <Y>
        scale(ur, ur, x_prime);
        scale(ut, ut, y_prime);
        sum(out, ur, ut);
    }

    /// Printing
    std::ostream &operator<<(std::ostream &s, const exposin &exps) {
        //s << std::setprecision(14) << "Exponential sinusoid object:" << std::endl;
        s << std::setprecision(16) << " k0: " << exps.m_k0;
        s << std::setprecision(16) << " k1: " << exps.m_k1;
        s << std::setprecision(16) << " k2: " << exps.m_k2;
        s << std::setprecision(16) << " phi: " << exps.m_phi;
        s << std::endl;
        return s;
    }
}