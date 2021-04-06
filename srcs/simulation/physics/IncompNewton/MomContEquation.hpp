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

class SIMULATION_API MomContEqIncompNewton : public Equation
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
        std::unique_ptr<PicardAlgo> m_pPicardAlgo;
        EigenSparseSolver m_solver;
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> m_solverIt;

        double m_rho;
        double m_mu;
        double m_gamma;
        double m_alpha;
        double m_Tr;

        Eigen::VectorXd m_bodyForce;

        //PSPG
        Eigen::SparseMatrix<double> m_A;
        Eigen::VectorXd m_b;

        void m_setupPicardPSPG(unsigned int maxIter, double minRes);

        void m_buildAbPSPG(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::VectorXd& qPrev);
        void m_applyBCPSPG(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::VectorXd& qPrev);

        //Fractionnal Step
        double m_gammaFS;

        Eigen::SparseMatrix<double> m_M;
        Eigen::SparseMatrix<double> m_MK_dt;
        Eigen::SparseMatrix<double> m_L;
        std::vector<Eigen::MatrixXd> m_DTelm; //One per element
        std::vector<Eigen::MatrixXd> m_Lelm; //One per element
        Eigen::VectorXd m_bVAppStep;
        Eigen::VectorXd m_bPcorrStep;
        Eigen::VectorXd m_bVStep;

        void m_setupPicardFracStep(unsigned int maxIter, double minRes);

        void m_buildMatFracStep(Eigen::SparseMatrix<double>& M,
                                Eigen::SparseMatrix<double>& MK_dt,
                                Eigen::SparseMatrix<double>& dtL,
                                std::vector<Eigen::MatrixXd>& DTelm,
                                std::vector<Eigen::MatrixXd>& dtLelm, Eigen::VectorXd& bVAppStep,
                                const std::vector<Eigen::VectorXd>& qPrev);
        void m_buildMatPcorrStep(Eigen::VectorXd& bPcorrStep, const Eigen::VectorXd& qVTilde, const Eigen::VectorXd& qPprev);
        void m_buildMatVStep(Eigen::VectorXd& bVStep, const Eigen::VectorXd& qDeltaP);
        void m_applyBCVAppStep(Eigen::SparseMatrix<double>& MK_dt, Eigen::VectorXd& bVAppStep, const Eigen::VectorXd& qPrev);
        void m_applyBCPCorrStep(Eigen::SparseMatrix<double>& L, Eigen::VectorXd& bPcorrStep);
        void m_applyBCVStep(Eigen::SparseMatrix<double>& M, Eigen::VectorXd& bVAppStep);

        double m_computeTauPSPG(const Element& element) const;
};

#endif // MOMCONTEQINCOMPNEWTON_HPP_INCLUDED
