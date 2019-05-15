/**
(C) Copyright 2019 DQ Robotics Developers

This file is part of DQ Robotics.

    DQ Robotics is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DQ Robotics is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with DQ Robotics.  If not, see <http://www.gnu.org/licenses/>.

Contributors:
- Murilo M. Marinho        (murilo@nml.t.u-tokyo.ac.jp)
*/

#ifndef DQ_CPLEX_INTERFACE_H
#define DQ_CPLEX_INTERFACE_H

#include <eigen3/Eigen/Dense>

using namespace Eigen;

class CplexInterface
{

public:
    /**
     * @brief solveCanonicalLinearProgramWithLinearConstraints
     * Solves the Canonical Linear Program in the form
     * min(g)  g'f
     * s.t.    Af=b;
     *         g \geq 0
     * using the dual-simplex algorithm.
     * @param f the n x 1 vector of the linear coeficients of the problem variables.
     * @param A the m x n matrix of equality constraints.
     * @param b the m x 1 value for the equality constraints.
     * @param silent
     * @return the optimal g
     */
    static VectorXd solveCanonicalLinearProgram(const MatrixXd &f, const MatrixXd& A, const MatrixXd &b, bool silent=true, bool* solved=nullptr);

    /**
     * @brief solveCanonicalQuadraticProgramWithLinearConstraints
     *   Solves the Canonical Quadratic Program With Linear Constraints in the form
     *   min(g)  g'Qg+c'g
     *   s.t.    Wg<=w;
     * @param Q the n x n matrix of the quadratic coeficitients of the decision variables.
     * @param c the n x 1 vector of the linear coeficients of the decision variables.
     * @param W the m x n matrix of inequality constraints.
     * @param w the m x 1 value for the inequality constraints.
     * @param silent
     * @return the optimal g
     */
    static VectorXd solveCanonicalQuadraticProgramWithLinearConstraints(const MatrixXd &Q, const MatrixXd &c, const MatrixXd &W, const MatrixXd &w, bool silent=true, bool* solved=nullptr);
};

#endif
