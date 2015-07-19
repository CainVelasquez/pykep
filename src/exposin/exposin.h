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

#include "../astro_constants.h"
#include "../serialization.h"
#include "../config.h"

#define PROPAGATION_STEPS 10
#define DEFAULT_ABSCISSAS 10

namespace kep_toolbox {
    /// Exponential Sinusoid Trajectory Model
    /**
     * The exposin class represents a single trajectory in the form of an exponential sinusoid curve in 3D space.
     *
     * @author Chris Andre (chrisandre01 _AT_ gmail.com)
     */
    class __KEP_TOOL_VISIBLE exposin {
    public:
        friend std::ostream &operator<<(std::ostream &, const exposin &);
        exposin(const double &k0 = 0, const double &k1 = 0, const double &k2 = 0, const double &phi = 0);
        void set(const double &k0, const double &k1, const double &k2, const double &phi);
        void projection(const array3D &r1, const array3D &r2, const bool &lw, const int &revs);
        double tof(const double &psi, const double &mu, const int &abscissas = DEFAULT_ABSCISSAS) const;
        void r_vec(array3D &out, const double &theta) const;
        void v_vec(array3D &out, const double &theta, const double &mu) const;
        void a_vec(array3D &out, const double &theta, const double &mu) const;
        const double &get_psi() const;
        const int &get_revs() const;
        const double &get_transfer_angle() const;
        double get_final_mass(const double &mu, const double &isp, const double &m=1.0) const;
        double get_maximum_thrust(const double &mu, const double &isp, const double &m) const;
        double get_delta_v(const double &mu) const;
    private:
        double tany(const double &theta) const;
        double r(const double &theta) const;
        double r_dot(const double &theta, const double &mu) const;
        double theta_dot(const double &theta, const double &mu) const;
        double v(const double &theta, const double &mu) const;
        double local_a(const double &theta) const;
        double a(const double &theta, const double &mu) const;
        void rel_coord(array3D &out, const array3D &z_axis, const array3D &x_axis, const double &x_prime,
                       const double &y_prime) const;
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive &ar, const unsigned int) {
            ar & m_k0;
            ar & m_k1;
            ar & m_k2;
            ar & m_phi;
            ar & m_r1;
            ar & m_r2;
            ar & m_transfer_angle;
            ar & m_psi;
            ar & m_lw;
            ar & m_revs;
        }
        double m_k0, m_k1, m_k2, m_phi;
        array3D m_r1, m_r2;
        double m_transfer_angle, m_psi;
        bool m_lw;
        int m_revs;
    };
    __KEP_TOOL_VISIBLE std::ostream &operator<<(std::ostream &, const exposin &);
}

#endif
