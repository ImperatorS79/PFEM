#pragma once
#ifndef MOMCONTEQINCOMPNEWTON_HPP_INCLUDED
#define MOMCONTEQINCOMPNEWTON_HPP_INCLUDED

#include <Eigen/IterativeLinearSolvers>
#ifdef EIGEN_USE_MKL_ALL
    #include <Eigen/PardisoSupport>
    typedef Eigen::PardisoLU<Eigen::SparseMatrix<double>> EigenSparseSolver;
#else
    #include <Eigen/SparseLU>
    typedef Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> EigenSparseSolver;
#endif

#include "../../Equation.hpp"
#include "../../nonLinearAlgo/PicardAlgo.hpp"

class Problem;
class Mesh;

template<unsigned short dim>
class MomContEqIncompNewton : public Equation
{
    public:
        MomContEqIncompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                              std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                              const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        MomContEqIncompNewton(const MomContEqIncompNewton& equation)             = delete;
        MomContEqIncompNewton& operator=(const MomContEqIncompNewton& equation)  = delete;
        MomContEqIncompNewton(MomContEqIncompNewton&& equation)                  = delete;
        MomContEqIncompNewton& operator=(MomContEqIncompNewton&& equation)       = delete;
        ~MomContEqIncompNewton() override;

        void displayParams() const override;

        bool solve() override;

    private:
        enum class Res
        {
            U,
            U_P,
            Ax_f
        };

        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder; /**< Class responsible of building the required matrices. */
        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder2; /**< Class responsible of building the required matrices. */
        std::unique_ptr<PicardAlgo> m_pPicardAlgo;
        EigenSparseSolver m_solver;
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> m_solverIt;
        Res m_residual;

        double m_rho;
        double m_mu;
        double m_gamma;
        double m_DgammaDT;
        double m_alpha;
        double m_Tr;
        double m_tau0;
        double m_mReg;

        bool m_phaseChange;
        double m_C;
        double m_eps;
        double m_Tm;
        double m_DT;

        Eigen::Matrix<double, dim, 1> m_bodyForce;

        //PSPG
        Eigen::SparseMatrix<double> m_A;
        Eigen::VectorXd m_b;

        void m_setupPicardPSPG(unsigned int maxIter, double minRes);

        void m_buildAbPSPG(const Eigen::VectorXd& qPrev);
        void m_applyBCPSPG(const Eigen::VectorXd& qPrev);

        //Fractionnal Step
        double m_gammaFS;

        Eigen::SparseMatrix<double> m_M;
        Eigen::SparseMatrix<double> m_MK_dt;
        Eigen::SparseMatrix<double> m_L;
        std::vector<Eigen::Matrix<double, dim*(dim + 1), dim + 1>> m_DTelm; //One per element
        std::vector<Eigen::Matrix<double, dim + 1, dim + 1>> m_Lelm; //One per element
        Eigen::VectorXd m_bVAppStep;
        Eigen::VectorXd m_bPcorrStep;
        Eigen::VectorXd m_bVStep;

        void m_setupPicardFracStep(unsigned int maxIter, double minRes);

        void m_buildMatFracStep(const std::vector<Eigen::VectorXd>& qPrev);
        void m_buildMatPcorrStep(const Eigen::VectorXd& qVTilde, const Eigen::VectorXd& qPprev);
        void m_buildMatVStep(const Eigen::VectorXd& qDeltaP);
        void m_applyBCVAppStep(const Eigen::VectorXd& qPrev);
        void m_applyBCPCorrStep();
        void m_applyBCVStep();

        double m_computeTauPSPG(const Element& element) const;
        double m_getFl(double T);
};

#include "MomContEquation.inl"
#include "MomContEquationPSPG.inl"
#include "MomContEquationFracStep.inl"

#endif // MOMCONTEQINCOMPNEWTON_HPP_INCLUDED
