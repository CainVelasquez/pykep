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

#ifndef KEP_TOOLBOX_CLASS_EXPOSIN_H
#define KEP_TOOLBOX_CLASS_EXPOSIN_H

#include <cmath>
#include "exposin.h"
#include "../astro_constants.h"
#include "../serialization.h"
#include "../config.h"
#include "boost/math/special_functions/fpclassify.hpp"

#define NR_MAX_ITERS 1000
#define NR_EPSILON 1.0e-5

namespace kep_toolbox {
    /// Exponential Sinusoid k2-Class
    /**
     * Represents the set of exponential sinusoids with a fixed k2 and some geometrical constraints
     */
    class __KEP_TOOL_VISIBLE class_exposin {
    private:
        /// Problem geometry and k2
        double k2, r1_m, r2_m, psi, angle;
        int multi_revs;
    public:
        class_exposin(const double &k2 = 0, const double &r1_m = 0, const double &r2_m = 0, const double &angle = 0,
                      const int &multi_revs = 0) {
            this->k2 = k2;
            this->r1_m = r1_m;
            this->r2_m = r2_m;
            this->angle = angle;
            this->psi = 2.0 * M_PI * multi_revs + angle;
            this->multi_revs = multi_revs;
        }

        /// Get total traversed angular distance (radians)
        const double &get_psi() const {
            return psi;
        }

        /// Get angular distance between specified radii
        const double &get_angle() const {
            return angle;
        }

        /// Set required revolutions about central body between initial and final states
        void set_revs(const int &revs) {
            psi = 2.0 * M_PI * revs + angle;
            multi_revs = revs;
        }

        /// Calculate minimum and maximum permissible tan y1, returns true if range is valid
        bool tany1_range(double &tany1_l, double &tany1_u) {
            double logr1r2 = log(r1_m / r2_m);
            double cosk2O = cos(k2 * psi);
            double delta = 2.0 * (1.0 - cosk2O) / pow(k2, 4.0) - logr1r2 * logr1r2;
            double tany1min = k2 / 2.0 * (-logr1r2 / tan(k2 * psi / 2) - sqrt(delta));
            double tany1max = k2 / 2.0 * (-logr1r2 / tan(k2 * psi / 2) + sqrt(delta));
            tany1_l = tany1min;
            tany1_u = tany1max;
            return delta >= 0.0;
        }

        /// Build an exposin instance according to a given tan y1
        void create_exposin(exposin &expsn, const double &tany1) const {
            double logr1r2 = log(r1_m / r2_m);
            double sink2O = sin(k2 * psi);
            double cosk2O = cos(k2 * psi);

            double k1_sqr = pow((logr1r2 + tany1 / k2 * sink2O) / (1.0 - cosk2O), 2.0) + pow(tany1 / k2, 2.0);
            double k1_sign = (logr1r2 + tany1 / k2 * sink2O) / (1.0 - cosk2O);

            double k1;
            if (k1_sign < 0) {
                k1 = -sqrt(k1_sqr);
            }
            else {
                k1 = sqrt(k1_sqr);
            }
            double phi = acos(tany1 / k1 / k2);
            double k0 = r1_m / exp(k1 * sin(phi));
            expsn.set(k0, k1, k2, phi);
        }

        /// Find the tany1 that results in a given time of flight. Returns false if none exists in valid range.
        bool search_tany1(double &tany1, const double &dT, const double &mu, const double &stop_tol = 1.0e3) {
            // Uses Newton-Raphson with approximated derivatives to find tany1.
            exposin exps;
            double tof_guess;
            int iters = NR_MAX_ITERS;
            double tany1_lb, tany1_ub;
            if (!tany1_range(tany1_lb, tany1_ub)) return false;
            double tany1_guess = 0.5 * (tany1_lb + tany1_ub);
            do {
                // Calculate d_tof / d_tany1
                create_exposin(exps, tany1_guess + NR_EPSILON);
                double dy = exps.tof(psi, mu) - tof_guess;
                double dydx = dy / NR_EPSILON;
                // Update guess
                tany1_guess += -(tof_guess - dT) / dydx;
                // Break if NaNs come up
                if (!boost::math::isfinite(tany1_guess)) return false;
                // Get new f(guess)
                create_exposin(exps, tany1_guess);
                tof_guess = exps.tof(psi, mu);
            } while (fabs(tof_guess - dT) > stop_tol && (--iters));
            if (iters == 0) return false;
            if (tany1_guess < tany1_lb || tany1_guess > tany1_ub) return false;
            tany1 = tany1_guess;
            return true;
        }

        /// Build an exposin instance from a given tof. Returns false at failure, true otherwise.
        bool tof_to_exposin(exposin &exps, const double &tof, const double &mu, const double &stop_tol = 1.0e4) {
            double candidate_tany1;
            if (!search_tany1(candidate_tany1, tof, mu, stop_tol)) return false;
            create_exposin(exps, candidate_tany1);
            return true;
        }

        /*
         * Helper methods and misc.
         */

    private:
        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive &ar, const unsigned int) {
            ar &k2;
            ar &r1_m;
            ar &r2_m;
            ar &psi;
            ar &angle;
            ar &multi_revs;
        }
    };
}

#endif
