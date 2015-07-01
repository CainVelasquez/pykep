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
#include "vectors.h"
#include "../serialization.h"
#include "../config.h"

namespace kep_toolbox {
    /// Lambert Solver for Exponential Sinusoid Trajectories
    /**
    * Represents lambert's boundary value problem for a class of trajectories
    * known as exponential sinusoids. These provide a sub-optimal but analytical
    * solution to the optimal control problem.
    *
    * See http://www.esa.int/gsp/ACT/doc/MAD/pub/ACT-TNT-MAD-LMSP01.pdf
    *
    * @author Chris Andre (chrisandre01 _AT_ gmail.com)
    */
    class __KEP_TOOL_VISIBLE lambert_exposin {
    private:
        class_exposin k2_class;
        std::vector<exposin> solv_exposins{};
        array3D r1, r2;
        std::vector<array3D> v1{}, v2{};
        int multi_revs;
        bool lw;
        double tof, mu, k2;
    public:
        lambert_exposin(const array3D &r1={0,0,0}, const array3D &r2={0,0,0}, const double &tof=0, const double &mu=0, const bool &lw=0,
                        const int &multi_revs = -1, const double &k2 = 0.1) {
            this->r1 = r1;
            this->r2 = r2;
            this->tof = tof;
            this->mu = mu;
            this->lw = lw;
            this->multi_revs = multi_revs;
            this->k2 = k2;

            double r1_m = norm(r1);
            double r2_m = norm(r2);
            double angle = acos(dot(r1, r2) / r1_m / r2_m);
            if (lw) {
                angle = 2.0 * M_PI - angle;
            }
            k2_class = class_exposin(k2, r1_m, r2_m, angle, multi_revs);

            if (multi_revs == -1) { // Detect maximum revolutions
                for (int n = 0; true; n++) {
                    if (build_a_solution(n)) {
                        this->multi_revs = n;
                    }
                    else break;
                }
            }
            else { // Force the specified revolution number
                if (!build_a_solution(this->multi_revs)) this->multi_revs = -1; // signal no solutions
            }
        }

        const std::vector<array3D>& get_v1() const {
            return v1;
        }

        const std::vector<array3D>& get_v2() const {
            return v2;
        }

        const array3D& get_r1() const {
            return r1;
        }

        const array3D& get_r2() const {
            return r2;
        }

        const double& get_tof() const {
            return tof;
        }

        const double& get_mu() const {
            return mu;
        }

        const int& get_revs() const {
            return multi_revs;
        }

        const std::vector<exposin>& get_exposins() const {
            return solv_exposins;
        }

        const double& get_traversal_final_mass(const double isp, const double m) const;

        const double& get_max_thrust(const double isp, const double m) const;

    private:
        bool build_a_solution(const int revs) {
            k2_class.set_revs(revs);
            exposin exps;
            if (!k2_class.tof_to_exposin(exps, tof, mu, 10*86400.0)) {
                return false;
            }
            array3D a, b;

            exps.projection(r1, r2, k2_class.get_angle(), lw, k2_class.get_psi(), revs);

            exps.v_vec(a, 0.0, mu);
            exps.v_vec(b, exps.get_psi(), mu);

            v1.push_back(a);
            v2.push_back(b);
            solv_exposins.push_back(exps);
            return true;
        }

        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
            ar & r1;
            ar & r2;
            ar & tof;
            ar & mu;
            ar & v1;
            ar & v2;
            ar & multi_revs;
            ar & lw;
            ar & k2;
            ar & k2_class;
            ar & solv_exposins;
        }
    public:
        friend std::ostream &operator<<(std::ostream &s, const lambert_exposin &lp) {
            s << std::setprecision(14) << "Lambert's problem (exponential sinusoid):" << std::endl;
            s << "mu = " << lp.mu << std::endl;
            s << "r1 = " << lp.r1 << std::endl;
            s << "r2 = " << lp.r2 << std::endl;
            s << "Time of flight: " << lp.tof <<std::endl<< std::endl;
            s << "Maximum number of revolutions: " << lp.multi_revs << std::endl;
            s << "Solutions: " << std::endl;
            for (int i = 0; i <= lp.multi_revs; i++){
                s << "Rev " << i << ": v1 := " << lp.get_v1()[i] << "; v2 := " << lp.get_v2()[i] << std::endl;
            }
            return s;
        }
    };
}

#endif
