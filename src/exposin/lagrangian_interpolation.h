#ifndef KEP_TOOLBOX_LAGRANGIAN_INTERPOLATION_H
#define KEP_TOOLBOX_LAGRANGIAN_INTERPOLATION_H

#include <vector>
#include <cmath>

namespace kep_toolbox {
    void poly_scale(std::vector<double> &out, std::vector<double> &poly, double scalar) {
        out.resize(poly.size());
        for (size_t i = 0; i < poly.size(); i++) {
            out[i] = poly[i] * scalar;
        }
    };

    void poly_multiply(std::vector<double> &out, std::vector<double> &a, std::vector<double> &b) {
        std::vector<double> tout(a.size() + b.size() - 1);
        for (size_t i = 0; i < a.size(); i++) {
            for (size_t j = 0; j < b.size(); j++) {
                int power = i + j;
                tout[power] += a[i] * b[j];
            }
        }
        out.resize(tout.size());
        for (size_t i = 0; i < tout.size(); i++) {
            out[i] = tout[i];
        }
    };

    void poly_add(std::vector<double> &out, std::vector<double> &a, std::vector<double> &b) {
        std::vector<double> large = a.size() > b.size() ? a : b;
        std::vector<double> small = a.size() > b.size() ? b : a;
        out.resize(large.size());
        for (size_t i = 0; i < small.size(); i++) {
            out[i] = small[i] + large[i];
        }
        for (size_t i = small.size(); i < large.size(); i++) {
            out[i] = large[i];
        }
    };

    double poly_eval(std::vector<double> &poly, double x) {
        double lin = 0.0;
        double x_i = 1.0;
        for (size_t i = 0; i < poly.size(); i++, x_i *= x) {
            lin += poly[i] * x_i;
        }
        return lin;
    };

    void poly_interpolate(std::vector<double> &out, std::vector<double> &x, std::vector<double> &y) {
        out.resize(1);
        out[0] = 0; // y(x) = 0
        for (size_t j = 0; j < x.size(); j++) {
            std::vector<double> numPoly = {1.0}; // y(x) = 1
            for (size_t i = 0; i < x.size(); i++) {
                if (i != j) {
                    std::vector<double> temp = {-x[i], 1}; // y(x) = x - x_m
                    poly_multiply(numPoly, numPoly, temp);
                }
            }
            double denomPoly = 1.0;
            for (size_t i = 0; i < x.size(); i++) {
                if (i != j) {
                    denomPoly *= x[j] - x[i];
                }
            }
            std::vector<double> scaled_basis;
            poly_scale(scaled_basis, numPoly, y[j] / denomPoly);
            poly_add(out, out, scaled_basis);
        }
    };
}

#endif
