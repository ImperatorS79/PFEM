#pragma once
#ifndef CONTEQWCOMPNEWTON_HPP_INCLUDED
#define CONTEQWCOMPNEWTON_HPP_INCLUDED

#include "../../Equation.hpp"

class Problem;
class Mesh;

template<unsigned short dim>
class ContEqWCompNewton : public Equation
{
    public:
        ContEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                          std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                          const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        ContEqWCompNewton(const ContEqWCompNewton& equation)             = delete;
        ContEqWCompNewton& operator=(const ContEqWCompNewton& equation)  = delete;
        ContEqWCompNewton(ContEqWCompNewton&& equation)                  = delete;
        ContEqWCompNewton& operator=(ContEqWCompNewton&& equation)       = delete;
        ~ContEqWCompNewton() override;

        void displayParams() const override;

        double getSquaredSpeedEquiv(const Node& node) const override;
        bool solve() override;
        void preCompute() override;

    private:
        double m_K0;
        double m_K0p;
        double m_rhoStar;
        double m_mu;
        double m_alpha;
        double m_Tr;
        Eigen::VectorXd m_bodyForce;
        bool m_strongContinuity;
        bool m_enableStab;
        double m_keStab;

        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder;
        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder2; // TO DO: change mat builder to allow multiple M functions

        Eigen::DiagonalMatrix<double,Eigen::Dynamic> m_invM;
        Eigen::VectorXd m_F0;

        void m_buildF0();
        void m_buildSystem();
        void m_applyBC();
        Eigen::VectorXd m_getPFromRhoTaitMurnagham(const Eigen::VectorXd& rho);

        void m_buildSystemFIC();
        void m_applyBCFIC();
        double m_computeTauFIC(const Element& element) const;
};

#include "ContEquation.inl"

#endif // CONTEQWCOMPNEWTON_HPP_INCLUDED
