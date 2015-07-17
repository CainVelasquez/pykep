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

#ifndef KEP_TOOLBOX_LAMBERT_EXPOSIN_H
#define KEP_TOOLBOX_LAMBERT_EXPOSIN_H

#include <vector>
#include <cmath>
#include "../core_functions/array3D_operations.h"
#include "class_exposin.h"
#include "exposin.h"
#include "../serialization.h"
#include "../config.h"

// Tolerance for TOF search
#define STOP_TOL 1.0e3

namespace kep_toolbox {
    /// Lambert Solver for Exponential Sinusoid Trajectories
    /**
    * Represents lambert's boundary value problem for a class of trajectories
    * known as exponential sinusoids. These provide a sub-optimal but analytical
    * solution to the optimal control problem.
    *
    * See http://www.esa.int/gsp/ACT/doc/MAD/pub/ACT-TNT-MAD-LMSP01.pdf
    * See http://www.esa.int/gsp/ACT/doc/ARI/ARI%20Study%20Report/ACT-RPT-MAD-ARI-05-4106-Spiral%20Trajectories%20in%20Global%20Optimisation.pdf
    *
    * @author Chris Andre (chrisandre01 _AT_ gmail.com)
    */
    class __KEP_TOOL_VISIBLE lambert_exposin {
    private:
        /// Represents class of exposins for given problem geometry in addition to a given k2
        class_exposin k2_class;

        /// Solution trajectories
        std::vector<exposin> solv_exposins{};

        /// Solution boundary velocities
        std::vector<array3D> v1{}, v2{};

        /// Specified boundary positions
        array3D r1, r2;

        /// Maximum number of revolutions for solution set
        int multi_revs;

        /// True if solutions exist
        bool any_solutions;

        bool lw;
        double tof, mu, k2;
    public:
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
        lambert_exposin(const array3D &r1 = {1, 0, 0}, const array3D &r2 = {0, 1, 0}, const double &tof = 1,
                        const double &mu = 1, const bool &lw = false,
                        const int &multi_revs = -1, const double &k2 = 0.1) {
            this->r1 = r1;
            this->r2 = r2;
            this->tof = tof;
            this->mu = mu;
            this->lw = lw;
            this->multi_revs = multi_revs;
            this->k2 = k2;

            // k2-class uses 2D polar coordinates + revs, so compute the angle mod 2PI
            double r1_m = norm(r1);
            double r2_m = norm(r2);
            double angle = acos(dot(r1, r2) / r1_m / r2_m);
            if (lw) {
                angle = 2.0 * M_PI - angle;
            }

            k2_class = class_exposin(k2, r1_m, r2_m, angle, multi_revs);

            if (multi_revs == -1) { // Detect maximum revolutions
                // This strategy will iterate through Nmin <= n <= Nmax, where Nmax,min are heuristic limits that should
                // encompass all good solutions. Once the first the solution is found, we continue iterating up until a
                // failing solution or Nmax is hit and stop.
                bool found_one = false;
                int min_revs = get_min_revs_heuristic();
                int max_revs = get_max_revs_heuristic();
                for (int n = min_revs; n <= max_revs; n++) {
                    if (build_a_solution(n)) {
                        this->multi_revs = n;
                        found_one = true;
                    }
                    else if (found_one) break; // stop search when solutions stop appearing
                }
            }
            else { // Force the specified revolution number
                if (!build_a_solution(multi_revs)) this->multi_revs = -1; // signal no solutions
            }
            any_solutions = this->multi_revs != -1;
        }

        const std::vector<array3D> &get_v1() const {
            return v1;
        }

        const std::vector<array3D> &get_v2() const {
            return v2;
        }

        const array3D &get_r1() const {
            return r1;
        }

        const array3D &get_r2() const {
            return r2;
        }

        const double &get_tof() const {
            return tof;
        }

        const double &get_mu() const {
            return mu;
        }

        const int &get_revs() const {
            return multi_revs;
        }

        const bool &has_solutions() const {
            return any_solutions;
        }

        const std::vector<exposin> &get_exposins() const {
            return solv_exposins;
        }

        /*
         * Helper methods and misc.
         */
    private:
        /// Helper function to fully construct and save a solution for a given revolution; returns false if none found.
        bool build_a_solution(const int revs) {
            // Try to build a valid exposin
            exposin exps;
            k2_class.set_revs(revs);
            if (!k2_class.tof_to_exposin(exps, tof, mu, STOP_TOL)) {
                return false;
            }
            // Exposin is valid solution...

            array3D a, b;
            // Add 3D interpretation data to exposin
            exps.projection(r1, r2, lw, revs);

            // Get boundary velocities
            exps.v_vec(a, 0.0, mu);
            exps.v_vec(b, exps.get_psi(), mu);

            v1.push_back(a);
            v2.push_back(b);
            solv_exposins.push_back(exps);

            return true;
        }

        /// Helper function to approximate the maximum revolutions for a given transfer problem that are feasible.
        int get_max_revs_heuristic() const {
            double period1 = 2.0 * M_PI * sqrt(pow(norm(r1), 3.0) / mu); // circular orbit period
            double period2 = 2.0 * M_PI * sqrt(pow(norm(r2), 3.0) / mu);
            return ceil(tof / fmin(period1, period2));
        }

        /// Helper function to approximate the minimum revolutions for a given transfer problem that are feasible.
        int get_min_revs_heuristic() const {
            double period1 = 2.0 * M_PI * sqrt(pow(norm(r1), 3.0) / mu); // circular orbit period
            double period2 = 2.0 * M_PI * sqrt(pow(norm(r2), 3.0) / mu);
            return floor(tof / fmax(period1, period2));
        }

        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive &ar, const unsigned int) {
            ar &r1;
            ar &r2;
            ar &tof;
            ar &mu;
            ar &v1;
            ar &v2;
            ar &multi_revs;
            ar &lw;
            ar &k2;
            ar &k2_class;
            ar &solv_exposins;
        }

        friend std::ostream &operator<<(std::ostream &s, const lambert_exposin &lp) {
            s << std::setprecision(7) << "Lambert's problem (exponential sinusoid):" << std::endl;
            s << "mu = " << lp.mu << std::endl;
            s << "r1 = " << lp.r1 << std::endl;
            s << "r2 = " << lp.r2 << std::endl;
            s << "angle = " << lp.k2_class.get_angle() << std::endl;
            s << "time of flight: " << lp.tof << std::endl;
            s << "has solutions? " << lp.has_solutions() << std::endl;
            if (lp.has_solutions()) {
                s << "maximum number of revolutions: " << lp.multi_revs << std::endl;
                s << "solutions: " << (int) lp.solv_exposins.size() << std::endl;
                for (int i = 0; i < (int) lp.solv_exposins.size(); i++) {
                    s << "Rev := " << lp.solv_exposins[i].get_revs() << ", v1 := " << lp.get_v1()[i] << ", v2 := " <<
                    lp.get_v2()[i] << std::endl;
                }
            }
            return s;
        }
    };
}

#endif
