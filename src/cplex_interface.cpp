#include<dqrobotics/interfaces/cplex_interface.h>

//Include CPLEX
#include <ilcplex/ilocplex.h>

VectorXd solveCanonicalLinearProgram(const MatrixXd &f, const MatrixXd &A, const MatrixXd &b, bool silent, bool* solved)
{
    const int PROBLEM_SIZE    = f.cols();
    const int CONSTRAINT_SIZE = b.cols();

    IloEnv env;

    //Each variable of g
    IloNumVarArray cplex_x(env);
    for(int i=0;i<PROBLEM_SIZE;i++)
        cplex_x.add(IloNumVar(env,0.0,IloInfinity));

    //Minizing function f
    IloObjective cplex_objective = IloMinimize(env); //The objective function
    IloNumArray cplex_c(env); //Each coeficient of f
    cplex_c.setSize(PROBLEM_SIZE);

    //Constraints for (A and b)
    IloRangeArray cplex_range(env); //Each constraint
    IloNumArray cplex_a(env);
    cplex_a.setSize(PROBLEM_SIZE);
    for(int i=0;i<CONSTRAINT_SIZE;i++)
        cplex_range.add(IloRange(env,0,0));

    //CPLEX_f
    for(int j=0;j<PROBLEM_SIZE;j++)
        cplex_c[j] = f(j);
    cplex_objective.setLinearCoefs(cplex_x,cplex_c);

    //CPLEX A and B
    for(int j=0;j<CONSTRAINT_SIZE;j++)
    {
        for(int k=0;k<PROBLEM_SIZE;k++)
            cplex_a[k] = A(j,k);
        cplex_range[j].setLinearCoefs(cplex_x,cplex_a);
        cplex_range[j].setBounds(b(j),b(j));
    }

    IloModel model(env); //Almost everything uses this
    model.add(cplex_objective);
    model.add(cplex_range);
    IloCplex cplex(model);
    if(silent)
    {
        cplex.setOut(env.getNullStream());
    }
    cplex.setParam(IloCplex::RootAlg,IloCplex::Dual);
    //cplex.setParam(IloCplex::Threads,1);

    if(! cplex.solve())
    {
        *solved = false;
    }
    *solved = true;
    IloNumArray outvalues(env);
    if(!silent){std::cout << "Solution status " << cplex.getStatus() << std::endl;}
    cplex.getValues(outvalues,cplex_x);

    VectorXd g(PROBLEM_SIZE);
    for(int j=0;j<PROBLEM_SIZE;j++)
    {
        g(j)=outvalues[j];
    }

    return g;
}


VectorXd solveCanonicalQuadraticProgramWithLinearConstraints(const MatrixXd &Q, const MatrixXd &c, const MatrixXd &W, const MatrixXd &w, bool silent, bool* solved)
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

