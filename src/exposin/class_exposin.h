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

#define REGULA_FALSI_ITERS 1000
#define MAX_TOF_LIMIT 100 * 365.0 * ASTRO_DAY2SEC
#define TANY1_HEURISTIC 0.95

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

        /// Find the tan y1 that results in a given time of flight using Regula Falsi.
        double search_tany1(const double &dT, const double &mu, const double &stop_tol = 1.0e3) {
            double tany1_a, tany1_b, tany1_c;
            if (!tany1_range(tany1_a, tany1_b)) return 1e20;
            tany1_a *= TANY1_HEURISTIC; // Boundaries have a few numerical issues, so avoid them
            tany1_b *= TANY1_HEURISTIC;
            double tof_a, tof_b, tof_c;

            exposin exps;

            // F(a)
            create_exposin(exps, tany1_a);
            tof_a = exps.tof(psi, mu);

            // F(b)
            create_exposin(exps, tany1_b);
            tof_b = exps.tof(psi, mu);

            // Max TOF can be problematically large, so impose a simple limit
            tof_b = fmin(tof_b, MAX_TOF_LIMIT);

            if (dT > tof_b || dT < tof_a) return 1e20; // Opportunistic prune

            double iters = REGULA_FALSI_ITERS;
            do {
                // Find c, then F(c)
                tany1_c = tany1_b - (tof_b - dT) * (tany1_b - tany1_a) / (tof_b - tof_a);
                create_exposin(exps, tany1_c);
                tof_c = exps.tof(psi, mu);

                // Revise bounds
                if ((tof_a > 0.0 && tof_c > 0.0) || (tof_a < 0.0 && tof_c < 0.0)) {
                    tany1_a = tany1_c;
                    tof_a = tof_c;
                }
                else {
                    tany1_b = tany1_c;
                    tof_b = tof_c;
                }
            }
            while (fabs(tof_c - dT) > stop_tol && (--iters));
            if (iters == 0) return 1e20;
            return tany1_c;
        }

        /// Build an exposin instance from a given tof. Returns false at failure, true otherwise.
        bool tof_to_exposin(exposin &exps, const double &tof, const double &mu, const double &stop_tol = 1.0e4) {
            double tany1_l, tany1_u;
            if (!tany1_range(tany1_l, tany1_u)) return false;
            double this_tany1 = search_tany1(tof, mu, stop_tol);
            if (this_tany1 > tany1_u || this_tany1 < tany1_l) return false;
            create_exposin(exps, this_tany1);
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
