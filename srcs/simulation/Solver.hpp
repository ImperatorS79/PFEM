#pragma once
#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED

#include <map>
#include <memory>
#include <vector>

#include "utility/Clock.hpp"
#include "utility/SolTable.hpp"

#include "simulation_defines.h"

class Problem;
class Mesh;
class Node;
class Equation;

/**
 * \class Solver
 * \brief Represents a solver (sets of equations) for a problem.
 *
 * The solver class is responsible for setting up the different equations needed by a
 * problem (based on the problem id and solver id) and to solve them step by step.
 */
class SIMULATION_API Solver
{
    public:
        Solver(Problem* pProblem, Mesh* pMesh, std::vector<SolTable> problemParams);
        Solver(const Solver& solver)             = delete;
        Solver& operator=(const Solver& solver)  = delete;
        Solver(Solver&& solver)                  = delete;
        Solver& operator=(Solver&& solver)       = delete;
        virtual ~Solver();

        /// \brief Display general parameters of the solver.
        virtual void displayParams() const;

        /// \brief Display time statistics.
        void displayTimeStats() const;

        bool checkBC(SolTable bcParam, unsigned int n, const Node& node, std::string bcString, unsigned int expectedBCSize);

        /// \return The id of the solver (child class have to set m_id).
        std::string getID() const noexcept;

        /// \brief Solve the equations for one time step and update the mesh if needed.
        /// \return true if the solver succeeded, false otherwise.
        virtual bool solveOneTimeStep();

        /// \brief Compute the next time step which will be used at the next call of solveOneTimeStep.
        virtual void computeNextDT();

        /// \return Get the number of additional states required by the solver with respect to the problem.
        virtual std::size_t getAdditionalStateCount() const;

        inline double getTimeStep() const noexcept;

    protected:
        std::string m_id;   /**<  The id of the solver (should be set by child class). */

        std::vector<SolTable> m_solverParams;   /**<  sol::table wrapper for the solver parameters (1 per thread). */
        double m_timeStep;                      /**< The current time step used. */

        std::vector<std::unique_ptr<Equation>> m_pEquations;    /**< Smart pointers to the equations */
        bool m_solveSucceed;    /**< Did the solveOneTimeStep function succeeded ? */

        Mesh* m_pMesh;          /**< Pointer to the mesh used. */
        Problem* m_pProblem;    /**< Pointer to the underlying problem. */

        std::map<std::string, double> m_accumalatedTimes;  /**< Accumulated time in each part of the solver*/
        Clock m_clock;

        double m_nextTimeToRemesh;
        void m_computeNextRemeshTime(bool force);
        void m_conditionalRemesh(std::size_t speedIndex);
};

#include "Solver.inl"

#endif // SOLVER_HPP_INCLUDED
