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

#include "exposin.h"

#include "../astro_constants.h"
#include "../serialization.h"
#include "../config.h"

#define NR_MAX_ITERS 20
#define NR_EPSILON 1.0e-5
#define DEFAULT_STOP_TOL 1.0e-3

namespace kep_toolbox {
    /// Exponential Sinusoid k2-Class
    /**
     * Represents the set of exponential sinusoids with a fixed k2 and some geometrical constraints
     *
     * @author Chris Andre (chrisandre01 _AT_ gmail.com)
     */
    class __KEP_TOOL_VISIBLE class_exposin {
    public:
        friend std::ostream &operator<<(std::ostream &, const class_exposin &);
        class_exposin(const double &k2 = 0, const double &R1 = 0, const double &R2 = 0, const double &transfer_angle = 0,
                      const int &revs = 0);
        const double &get_k2() const;
        const double &get_R1() const;
        const double &get_R2() const;
        const double &get_transfer_angle() const;
        const double &get_psi() const;
        void set_revs(const int &revs);
        bool tany1_range(double &tany1_lb, double &tany1_ub);
        void create_exposin(exposin &expsn, const double &tany1) const;
        bool search_tany1(double &tany1, const double &dT, const double &mu, int &iters, const double &stop_tol = DEFAULT_STOP_TOL);
        bool tof_to_exposin(exposin &exps, const double &tof, const double &mu, int &iters, const double &stop_tol = DEFAULT_STOP_TOL);
    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive &ar, const unsigned int) {
            ar & const_cast<double&> (m_k2);
            ar & const_cast<double&> (m_R1);
            ar & const_cast<double&> (m_R2);
            ar & const_cast<double&> (m_transfer_angle);
            ar & m_revs;
            ar & m_psi;
        }
        const double m_k2, m_R1, m_R2, m_transfer_angle;
        int m_revs;
        double m_psi;
    };
    __KEP_TOOL_VISIBLE std::ostream &operator<<(std::ostream &, const class_exposin &);
}

#endif
