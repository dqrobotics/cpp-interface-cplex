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

#ifndef DQ_CPLEXSOLVER_H
#define DQ_CPLEXSOLVER_H

#include <dqrobotics/solvers/DQ_QuadraticProgrammingSolver.h>

#include <eigen3/Eigen/Dense>
#include <ilcplex/ilocplex.h>

using namespace Eigen;

namespace DQ_robotics
class DQ_CPLEXSolver: public DQ_QuadraticProgrammingSolver
{
private:
bool show_output_;


public:

DQ_CPLEXSolver()
{
    show_output_ = false;
}
~DQ_CPLEXSolver()=default;

/**
     * @brief solveCanonicalQuadraticProgramWithLinearinequality_constraints
     *   Solves the Canonical Quadratic Program With Linear inequality_constraints in the form
     *   min(x)  x'Hx+f'x
     *   s.t.    Ax<=b
     *           Aeqx=beq
     * @param H the n x n matrix of the quadratic coeficitients of the decision variables.
     * @param f the n x 1 vector of the linear coeficients of the decision variables.
     * @param A the m x n matrix of inequality inequality_constraints.
     * @param b the m x 1 value for the inequality inequality_constraints.
     * @param Aeq the m x n matrix of equality inequality_constraints.
     * @param beq the m x 1 value for the inequality inequality_constraints.
     * @return the optimal x
     */
VectorXd solve_quadratic_program(const MatrixXd &H, const MatrixXd &f, const MatrixXd &A, const MatrixXd &b, const MatrixXd& Aeq, const MatrixXd& beq)
{
    const int PROBLEM_SIZE = H.rows();
    const int INEQUALITY_CONSTRAINT_SIZE = b.size();
    const int EQUALITY_CONSTRAINT_SIZE = beq.size();

    ///Create optimization environment
    IloEnv env;

    ///Decision variables
    IloNumVarArray decision_variables(env);
    for(int i=0;i<PROBLEM_SIZE;i++)
        decision_variables.add(IloNumVar(env,-IloInfinity,+IloInfinity,ILOFLOAT));

    ///Objective function
    //QP in the form 0.5*x'*H*x+f'x
    IloObjective objective = IloMinimize(env); //The objective function

    IloNumArray linear_coeficient(env); //Each coeficient of f
    linear_coeficient.setSize(PROBLEM_SIZE);

    //f'x
    for(int i=0;i<PROBLEM_SIZE;i++)
        linear_coeficient[i] = c(i);
    objective.setLinearCoefs(decision_variables,linear_coeficient);

    //0.5*x'*H*x
    for(int i=0;i<PROBLEM_SIZE;i++)
    {
        for(int j=0;j<PROBLEM_SIZE;j++)
        {
            if(i!=j)
                objective.setQuadCoef(decision_variables[i],decision_variables[j],Q(i,j)+Q(j,i));
            else
                objective.setQuadCoef(decision_variables[i],decision_variables[j],Q(i,j));
        }
    }

    ///Linear inequalities
    IloRangeArray inequality_constraints(env);
    for(int i=0;i<CONSTRAINT_SIZE;i++)
        inequality_constraints.add(IloRange(env,-IloInfinity,b(i)));

    //Ax <= b
    IloNumArray cplex_A(env);
    cplex_A.setSize(PROBLEM_SIZE);
    for(int j=0;j<INEQUALITY_CONSTRAINT_SIZE;j++)
    {
        for(int k=0;k<PROBLEM_SIZE;k++)
            cplex_A[k] = A(j,k);
        inequality_constraints[j].setLinearCoefs(decision_variables,cplex_A);
    }

    ///Linear inequalities
    IloRangeArray equality_constraints(env);
    for(int i=0;i<EQUALITY_CONSTRAINT_SIZE;i++)
        equality_constraints.add(IloRange(env,beq(i),beq(i)));

    //Aeqx = beq
    IloNumArray cplex_Aeq(env);
    cplex_Aeq.setSize(PROBLEM_SIZE);
    for(int j=0;j<EQUALITY_CONSTRAINT_SIZE;j++)
    {
        for(int k=0;k<PROBLEM_SIZE;k++)
            cplex_Aeq[k] = Aeq(j,k);
        equality_constraints[j].setLinearCoefs(decision_variables,cplex_Aeq);
    }

    IloModel model(env);
    model.add(objective);
    model.add(inequality_constraints);
    model.add(equality_constraints);

    IloCplex cplex(model);

    ///Settings I found to give the fastest solving for the conditions I tested with.
    ///Single thread, deterministic mode (both mean the same thing basically)
    ///Dual simplex as solver
    cplex.setParam(IloCplex::RootAlg,IloCplex::Dual);
    cplex.setParam(IloCplex::Threads,1);
    cplex.setParam(IloCplex::ParallelMode,IloCplex::Deterministic);

    if(not show_output_)
    {
        cplex.setOut(env.getNullStream());
    }

    if(!cplex.solve())
    {
        throw std::runtime_error("Unable to solve quadratic program!!");
    }


    IloNumArray outvalues(env);
    if(show_output_){std::cout << "Solution status " << cplex.getStatus() << std::endl;}
    cplex.getValues(outvalues,decision_variables);

    VectorXd x(PROBLEM_SIZE);
    for(int j=0;j<PROBLEM_SIZE;j++)
    {
        x(j)=outvalues[j];
    }

    //End enviroment
    env.end();

    return g;
}
};

#endif
