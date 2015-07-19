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
#include "class_exposin.h"
#include "exposin.h"

#include "../astro_constants.h"
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
        static const array3D default_r1;
        static const array3D default_r2;
    public:
        friend std::ostream &operator<<(std::ostream &, const lambert_exposin &);
        lambert_exposin(const array3D &r1 = default_r1, const array3D &r2 = default_r2, const double &tof = 1,
                        const double &mu = 1, const bool &lw = false,
                        const int &revs = -1, const double &k2 = 0.1);
        const std::vector<array3D> &get_v1() const;
        const std::vector<array3D> &get_v2() const;
        const array3D &get_r1() const;
        const array3D &get_r2() const;
        const double &get_tof() const;
        const double &get_mu() const;
        const double &get_transfer_angle() const;
        const double &get_k2() const;
        const bool &has_solutions() const;
        const int &num_solutions() const;
        const int &min_revs() const;
        const int &max_revs() const;
        const std::vector<exposin> &get_exposins() const;
    private:
        bool build_a_solution(class_exposin &k2_class, const int revs);
        int get_max_revs_heuristic() const;
        int get_min_revs_heuristic() const;
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive &ar, const unsigned int) {
            ar & const_cast<array3D&> (m_r1);
            ar & const_cast<array3D&> (m_r2);
            ar & const_cast<double&> (m_tof);
            ar & const_cast<double&> (m_mu);
            ar & const_cast<double&> (m_k2);
            ar & const_cast<bool&> (m_lw);
            ar & const_cast<double&> (m_transfer_angle);
            ar & m_min_revs;
            ar & m_max_revs;
            ar & m_has_solutions;
            ar & m_num_solutions;
        }
        std::vector<exposin> m_solv_exposins;
        std::vector<array3D> m_v1, m_v2;
        const array3D m_r1, m_r2;
        int m_min_revs, m_max_revs;
        bool m_has_solutions;
        const double m_tof, m_mu, m_k2;
        const bool m_lw;
        const double m_transfer_angle;
        int m_num_solutions;
    };
    __KEP_TOOL_VISIBLE std::ostream &operator<<(std::ostream &, const lambert_exposin &);
}

#endif
