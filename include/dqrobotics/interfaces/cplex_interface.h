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
#include <ilcplex/ilocplex.h>

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
    static VectorXd solveCanonicalLinearProgram(const MatrixXd &f, const MatrixXd& A, const MatrixXd &b, bool silent=true, bool* solved=nullptr)
    {
    
{
    const int PROBLEM_SIZE    = Q.rows();
    const int CONSTRAINT_SIZE = w.rows();

    ///Create optimization environment
    IloEnv env;

    ///Optimized variables
    IloNumVarArray decision_variables(env);
    for(int i=0;i<PROBLEM_SIZE;i++)
        decision_variables.add(IloNumVar(env,-IloInfinity,+IloInfinity,ILOFLOAT));

    ///Objective function
    ///QP in the form 0.5*x'*Q*x+c'x
    IloObjective objective = IloMinimize(env); //The objective function

    IloNumArray linear_coeficient(env); //Each coeficient of c
    linear_coeficient.setSize(PROBLEM_SIZE);

    //c'x
    for(int i=0;i<PROBLEM_SIZE;i++)
        linear_coeficient[i] = c(i);
    objective.setLinearCoefs(decision_variables,linear_coeficient);

    //0.5*x'*Q*x
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
    IloRangeArray constraints(env);
    for(int i=0;i<CONSTRAINT_SIZE;i++)
        constraints.add(IloRange(env,-IloInfinity,w(i)));

    //Wg <= w
    IloNumArray cplex_W(env);
    cplex_W.setSize(PROBLEM_SIZE);
    for(int j=0;j<CONSTRAINT_SIZE;j++)
    {
        for(int k=0;k<PROBLEM_SIZE;k++)
            cplex_W[k] = W(j,k);
        constraints[j].setLinearCoefs(decision_variables,cplex_W);
    }

    IloModel model(env);
    model.add(objective);
    model.add(constraints);

    IloCplex cplex(model);

    ///Settings I found to give the fastest solving for the conditions I tested with.
    ///Single thread, deterministic mode (both mean the same thing basically)
    ///Dual simplex as solver
    cplex.setParam(IloCplex::RootAlg,IloCplex::Dual);
    cplex.setParam(IloCplex::Threads,1);
    cplex.setParam(IloCplex::ParallelMode,IloCplex::Deterministic);

    if(silent)
    {
        cplex.setOut(env.getNullStream());
    }

    if(! cplex.solve())
    {
        if(solved!=nullptr){*solved = false;}
    }
    if(solved!=nullptr){*solved = true;}

    IloNumArray outvalues(env);
    if(!silent){std::cout << "Solution status " << cplex.getStatus() << std::endl;}
    cplex.getValues(outvalues,decision_variables);

    VectorXd g(PROBLEM_SIZE);
    for(int j=0;j<PROBLEM_SIZE;j++)
    {
        g(j)=outvalues[j];
    }

    //End enviroment
    env.end();

    return g;
}

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
    static VectorXd solveCanonicalQuadraticProgramWithLinearConstraints(const MatrixXd &Q, const MatrixXd &c, const MatrixXd &W, const MatrixXd &w, bool silent=true, bool* solved=nullptr)
    {
        const int PROBLEM_SIZE    = Q.rows();

    ///Create optimization environment
    IloEnv env;

    ///Optimized variables
    IloNumVarArray decision_variables(env);
    for(int i=0;i<PROBLEM_SIZE;i++)
        decision_variables.add(IloNumVar(env,-IloInfinity,+IloInfinity,ILOFLOAT));

    ///Objective function
    ///QP in the form 0.5*x'*Q*x+c'x
    IloObjective objective = IloMinimize(env); //The objective function

    IloNumArray linear_coeficient(env); //Each coeficient of c
    linear_coeficient.setSize(PROBLEM_SIZE);

    //c'x
    for(int i=0;i<PROBLEM_SIZE;i++)
        linear_coeficient[i] = c(i);
    objective.setLinearCoefs(decision_variables,linear_coeficient);

    //0.5*x'*Q*x
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

    IloModel model(env);
    model.add(objective);

    IloCplex cplex(model);

    ///Settings I found to give the fastest solving for the conditions I tested with.
    ///Single thread, deterministic mode (both mean the same thing basically)
    ///Dual simplex as solver
    cplex.setParam(IloCplex::RootAlg,IloCplex::Dual);
    cplex.setParam(IloCplex::Threads,1);
    cplex.setParam(IloCplex::ParallelMode,IloCplex::Deterministic);

    if(silent)
    {
        cplex.setOut(env.getNullStream());
    }

    if(!cplex.solve())
    {
        if(solved!=nullptr){*solved = false;}
        else
        {
            throw int(0);
        }
    }
    if(solved!=nullptr){*solved = true;}

    IloNumArray outvalues(env);
    if(!silent){std::cout << "Solution status " << cplex.getStatus() << std::endl;}
    cplex.getValues(outvalues,decision_variables);

    VectorXd g(PROBLEM_SIZE);
    for(int j=0;j<PROBLEM_SIZE;j++)
    {
        g(j)=outvalues[j];
    }

    //End enviroment
    env.end();

    return g;
    }
};

#endif
