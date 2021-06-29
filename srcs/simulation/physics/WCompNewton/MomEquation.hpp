#pragma once
#ifndef MOMEQWCOMPNEWTON_HPP_INCLUDED
#define MOMEQWCOMPNEWTON_HPP_INCLUDED

#include "../../Equation.hpp"

class Problem;
class Mesh;

template<unsigned short dim>
class MomEqWCompNewton : public Equation
{
    public:
        MomEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                          std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                          const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        MomEqWCompNewton(const MomEqWCompNewton& equation)             = delete;
        MomEqWCompNewton& operator=(const MomEqWCompNewton& equation)  = delete;
        MomEqWCompNewton(MomEqWCompNewton&& equation)                  = delete;
        MomEqWCompNewton& operator=(MomEqWCompNewton&& equation)       = delete;
        ~MomEqWCompNewton() override;

        void displayParams() const override;

        double getSquaredSpeedEquiv(const Node& node) const override;
        double getDiffusionParam(const Node& node) const override;
        bool solve() override;
        void setQVhalf(Eigen::VectorXd qV1half);

    private:
        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder;
        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder2;

        double m_mu;
        double m_gamma;
        double m_alpha;
        double m_Tr;
        double m_DgammaDT;

        bool m_phaseChange;
        double m_C;
        double m_eps;
        double m_Tm;
        double m_DT;

        Eigen::Matrix<double, dim, 1> m_bodyForce;

        Eigen::DiagonalMatrix<double,Eigen::Dynamic> m_invM;
        Eigen::VectorXd m_F;

        void m_buildSystem();
        void m_applyBC();

        double m_getFl(double T);
};

#include "MomEquation.inl"

#endif // MOMEQWCOMPNEWTON_HPP_INCLUDED
