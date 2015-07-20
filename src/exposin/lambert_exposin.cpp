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

#include "lambert_exposin.h"

#include "../core_functions/array3D_operations.h"
#include <cmath>
#include "class_exposin.h"
#include "exposin.h"
#include <vector>

#include "../astro_constants.h"
#include "../serialization.h"
#include "../config.h"

namespace kep_toolbox {

    const array3D lambert_exposin::default_r1 = {{1.0, 0.0, 0.0}};
    const array3D lambert_exposin::default_r2 = {{0.0, 1.0, 0.0}};

    /// Constructor
    /** Constructs and solves a Lambert BVP for exponential sinusoids.
     *
     * \param[in] r1 first cartesian position
     * \param[in] r2 second cartesian position
     * \param[in] tof time of flight
     * \param[in] mu gravity parameter
     * \param[in] lw should long way between r1, r2 be taken in solutions
     * \param[in] multi_revs -1 to autodetect all valid revolutions, N to force a specific revolution N
     * \param[in] k2 winding parameter
     */
    lambert_exposin::lambert_exposin(const array3D &r1, const array3D &r2, const double &tof,
                                     const double &mu, const bool &lw,
                                     const int &revs, const double &k2) : m_r1(r1), m_r2(r2), m_tof(tof),
                                                                                m_mu(mu), m_k2(k2), m_lw(lw),
                                                                                m_transfer_angle(angle(r1, r2, lw)),
                                                                                m_num_solutions(0) {

        double R1 = norm(m_r1);
        double R2 = norm(m_r2);

        class_exposin k2_class = class_exposin(m_k2, R1, R2, m_transfer_angle);

        if (revs == -1) { // Detect maximum revolutions
            // This strategy will iterate through Nmin <= n <= Nmax, where Nmax,min are heuristic limits that should
            // encompass all good solutions. Once the first the solution is found, we continue iterating up until a
            // failing solution or Nmax is hit and stop.
            bool found_one = false;
            int min_revs = get_min_revs_heuristic();
            int max_revs = get_max_revs_heuristic();
            for (int n = min_revs; n <= max_revs; n++) {
                if (build_a_solution(k2_class, n)) { // found a solution at n
                    if (!found_one) {
                        found_one = true;
                        m_min_revs = n;
                    }
                    m_max_revs = n;
                }
                else if (found_one) { // stop search when solutions stop appearing
                    break;
                }
            }
        }
        else { // Force the specified revolution number
            build_a_solution(k2_class, revs);
        }
        m_has_solutions = m_num_solutions != 0;
    }

    /// Get boundary velocities at r1
    const std::vector<array3D> &lambert_exposin::get_v1() const {
        return m_v1;
    }

    /// Get boundary velocities at r2
    const std::vector<array3D> &lambert_exposin::get_v2() const {
        return m_v2;
    }

    /// Get starting position
    const array3D &lambert_exposin::get_r1() const {
        return m_r1;
    }

    /// Get ending position
    const array3D &lambert_exposin::get_r2() const {
        return m_r2;
    }

    /// Get time of flight
    const double &lambert_exposin::get_tof() const {
        return m_tof;
    }

    /// Get gravitational parameter
    const double &lambert_exposin::get_mu() const {
        return m_mu;
    }

    /// Get transfer angle
    const double &lambert_exposin::get_transfer_angle() const {
        return m_transfer_angle;
    }

    /// Get winding parameter k2
    const double &lambert_exposin::get_k2() const {
        return m_k2;
    }

    /// Return whether any solutions have been found for the given problem
    const bool &lambert_exposin::has_solutions() const {
        return m_has_solutions;
    }

    /// Get the smallest revolution number of any solution found
    const int &lambert_exposin::min_revs() const {
        return m_min_revs;
    }

    /// Get the largest revolution number of any solution found
    const int &lambert_exposin::max_revs() const {
        return m_max_revs;
    }

    /// Get the number of solutions found for the given problem
    const int &lambert_exposin::num_solutions() const {
        return m_num_solutions;
    }

    /// Get the exposin curves that solve the given problem
    const std::vector<exposin> &lambert_exposin::get_exposins() const {
        return m_solv_exposins;
    }

    /// Helper function to fully construct and save a solution for a given revolution; returns false if none found.
    bool lambert_exposin::build_a_solution(class_exposin &k2_class, const int revs) {
        k2_class.set_revs(revs);
        // Try to build a valid exposin
        m_solv_exposins.resize((unsigned long)m_num_solutions + 1); // allocate one to test
        if (!k2_class.tof_to_exposin(m_solv_exposins[m_num_solutions], m_tof, m_mu, m_tof * TOF_FRACTION)) {
            m_solv_exposins.resize((unsigned long)m_num_solutions); // deallocate at failure
            return false;
        }
        // Exposin is valid solution...
        m_v1.resize((unsigned long)m_num_solutions + 1);
        m_v2.resize((unsigned long)m_num_solutions + 1);

        // Add 3D interpretation data to exposin
        m_solv_exposins[m_num_solutions].projection(m_r1, m_r2, m_lw, revs);

        // Get boundary velocities
        double psi = m_solv_exposins[m_num_solutions].get_psi();
        m_solv_exposins[m_num_solutions].v_vec(m_v1[m_num_solutions], 0.0, m_mu);
        m_solv_exposins[m_num_solutions].v_vec(m_v2[m_num_solutions], psi, m_mu);

        m_num_solutions++;
        return true;
    }

    /// Helper function to approximate the maximum revolutions for a given transfer problem that are feasible.
    int lambert_exposin::get_max_revs_heuristic() const {
        double period1 = 2.0 * M_PI * sqrt(pow(norm(m_r1), 3.0) / m_mu); // circular orbit period
        double period2 = 2.0 * M_PI * sqrt(pow(norm(m_r2), 3.0) / m_mu);
        return (int) ceil(m_tof / fmin(period1, period2));
    }

    /// Helper function to approximate the minimum revolutions for a given transfer problem that are feasible.
    int lambert_exposin::get_min_revs_heuristic() const {
        double period1 = 2.0 * M_PI * sqrt(pow(norm(m_r1), 3.0) / m_mu); // circular orbit period
        double period2 = 2.0 * M_PI * sqrt(pow(norm(m_r2), 3.0) / m_mu);
        return (int) floor(m_tof / fmax(period1, period2));
    }

    /// Printing
    std::ostream &operator<<(std::ostream &s, const lambert_exposin &lp) {
        s << "Lambert's problem (exponential sinusoid)" << std::endl;
        s << "mu: " << std::setprecision(16) << lp.m_mu << std::endl;
        s << "tof: " << std::setprecision(16) << lp.m_tof << std::endl;
        s << "k2: " << std::setprecision(16) << lp.m_k2 << std::endl;
        s << "r1: " << lp.m_r1 << std::endl;
        s << "r2: " << lp.m_r2 << std::endl;
        s << "lw: " << lp.m_lw << std::endl;
        s << "transfer angle: " << std::setprecision(16) << lp.m_transfer_angle << std::endl;
        s << "has solutions? " << lp.m_has_solutions << std::endl;
        if (lp.has_solutions()) {
            s << "num solutions: " << lp.m_num_solutions << std::endl;
            s << "min revs: " << lp.m_min_revs << std::endl;
            s << "max revs: " << lp.m_max_revs << std::endl;
            s << "solutions: " << lp.m_num_solutions << std::endl;
            for (int i = 0; i < lp.m_num_solutions; i++) {
                s << "Rev: " << lp.m_solv_exposins[i].get_revs() << ", v1: " << lp.get_v1()[i] << ", v2: " <<
                lp.get_v2()[i] << std::endl;
                s << "Exposin: " << lp.m_solv_exposins[i] << std::endl;
            }
        }
        return s;
    }
} // namespace
