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

#include "class_exposin.h"

#include <cmath>
#include "exposin.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "../astro_constants.h"
#include "../serialization.h"
#include "../config.h"

namespace kep_toolbox {
    class_exposin::class_exposin(const double &k2, const double &R1, const double &R2, const double &transfer_angle,
                                 const int &revs) : m_k2(k2), m_R1(R1), m_R2(R2), m_transfer_angle(transfer_angle),
                                                    m_revs(revs), m_psi(revs * 2 * M_PI) { }

    /// Get k2
    const double &class_exposin::get_k2() const {
        return m_k2;
    }

    /// Get R1
    const double &class_exposin::get_R1() const {
        return m_R1;
    }

    /// Get R2
    const double &class_exposin::get_R2() const {
        return m_R2;
    }

    /// Get transfer angle
    const double &class_exposin::get_transfer_angle() const {
        return m_transfer_angle;
    }

    /// Get psi
    const double &class_exposin::get_psi() const {
        return m_psi;
    }

    /// Set required revolutions about central body between initial and final states
    void class_exposin::set_revs(const int &revs) {
        m_psi = 2.0 * M_PI * revs + m_transfer_angle;
        m_revs = revs;
    }

    /// Calculate minimum and maximum permissible tan y1, returns true if range is valid
    bool class_exposin::tany1_range(double &tany1_lb, double &tany1_ub) {
        double logr1r2 = log(m_R1 / m_R2);
        double cos_k2_theta = cos(m_k2 * m_psi);
        double delta = 2.0 * (1.0 - cos_k2_theta) / pow(m_k2, 4.0) - logr1r2 * logr1r2;
        tany1_lb = m_k2 / 2.0 * (-logr1r2 / tan(m_k2 * m_psi / 2.0) - sqrt(delta));
        tany1_ub = m_k2 / 2.0 * (-logr1r2 / tan(m_k2 * m_psi / 2.0) + sqrt(delta));
        return delta >= 0.0;
    }

    /// Build an exposin instance according to a given tany1
    void class_exposin::create_exposin(exposin &expsn, const double &tany1) const {
        double logr1r2 = log(m_R1 / m_R2);
        double sin_k2_theta = sin(m_k2 * m_psi);
        double cos_k2_theta = cos(m_k2 * m_psi);

        double k1_sqr =
                pow((logr1r2 + tany1 / m_k2 * sin_k2_theta) / (1.0 - cos_k2_theta), 2.0) + pow(tany1 / m_k2, 2.0);
        double k1_sign = (logr1r2 + tany1 / m_k2 * sin_k2_theta) / (1.0 - cos_k2_theta);

        double k1;
        if (k1_sign < 0) {
            k1 = -sqrt(k1_sqr);
        }
        else {
            k1 = sqrt(k1_sqr);
        }
        double phi = acos(tany1 / k1 / m_k2);
        double k0 = m_R1 / exp(k1 * sin(phi));
        expsn.set(k0, k1, m_k2, phi);
    }

    /// Find the tany1 that results in a given time of flight. Returns false if none exists in valid range.
    bool class_exposin::search_tany1(double &tany1, const double &dT, const double &mu, const double &stop_tol) {
        // Uses Newton-Raphson with approximated derivatives to find tany1.
        exposin exps;

        // Check that candidate sinusoids exist
        double tany1_lb, tany1_ub;
        if (!tany1_range(tany1_lb, tany1_ub)) return false;

        // Initial x, f(x)
        double tany1_guess = 0.5 * (tany1_lb + tany1_ub);
        create_exposin(exps, tany1_guess);
        double tof_guess = exps.tof(m_psi, mu);

        int iters = NR_MAX_ITERS;
        do {
            // Calculate d_tof / d_tany1
            create_exposin(exps, tany1_guess + NR_EPSILON);
            double dy = exps.tof(m_psi, mu) - tof_guess;
            double dydx = dy / NR_EPSILON;
            // Update guess
            tany1_guess += -(tof_guess - dT) / dydx;
            // Break if NaN comes up
            if (!boost::math::isfinite(tany1_guess)) return false;
            // Get new f(x)
            create_exposin(exps, tany1_guess);
            tof_guess = exps.tof(m_psi, mu);
        } while (fabs(tof_guess - dT) > stop_tol && (--iters));

        if (iters == 0) return false;
        if (tany1_guess < tany1_lb || tany1_guess > tany1_ub) return false;
        tany1 = tany1_guess;
        return true; // success!
    }

    /// Build an exposin instance from a given tof. Returns false at failure, true otherwise.
    bool class_exposin::tof_to_exposin(exposin &exps, const double &tof, const double &mu, const double &stop_tol) {
        double candidate_tany1;
        if (!search_tany1(candidate_tany1, tof, mu, stop_tol)) return false;
        create_exposin(exps, candidate_tany1);
        return true;
    }

    /// Printing
    std::ostream &operator<<(std::ostream &s, const class_exposin &ce) {
        s << "k2 Exposin Class" << std::endl;
        s << "k2: " << std::setprecision(16) << ce.m_k2 << std::endl;
        s << "R1: " << std::setprecision(16) << ce.m_R1 << std::endl;
        s << "R2: " << std::setprecision(16) << ce.m_R2 << std::endl;
        s << "transfer angle: " << std::setprecision(16) << ce.m_transfer_angle << std::endl;
        s << "psi: " << std::setprecision(16) << ce.m_psi << std::endl;
        s << "revs: " << ce.m_revs << std::endl;
        return s;
    }
};