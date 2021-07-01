#pragma once
#ifndef HEATEQWCOMPNEWTON_HPP_INCLUDED
#define HEATEQWCOMPNEWTON_HPP_INCLUDED

#include "../../Equation.hpp"

class Problem;
class Mesh;

template<unsigned short dim>
class HeatEqWCompNewton : public Equation
{
    public:
        HeatEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                          std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                          const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        HeatEqWCompNewton(const HeatEqWCompNewton& equation)             = delete;
        HeatEqWCompNewton& operator=(const HeatEqWCompNewton& equation)  = delete;
        HeatEqWCompNewton(HeatEqWCompNewton&& equation)                  = delete;
        HeatEqWCompNewton& operator=(HeatEqWCompNewton&& equation)       = delete;
        ~HeatEqWCompNewton() override;

        void displayParams() const override;

        double getDiffusionParam(const Node& node) const override;
        bool solve() override;

    private:
        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder;
        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder2;

        double m_k;
        double m_cv;
        double m_Tm;
        double m_DT;
        double m_Lm;
        double m_h;
        double m_Tinf;
        double m_epsRad;
        bool m_phaseChange;

        Eigen::DiagonalMatrix<double,Eigen::Dynamic> m_invM;
        Eigen::VectorXd m_F;

        void m_buildSystem();
        void m_applyBC();

        double m_getDflDT(double T);
};

#include "HeatEquation.inl"

#endif // HEATEQWCOMPNEWTON_HPP_INCLUDED
