#pragma once
#ifndef HEATEQINCOMPNEWTON_HPP_INCLUDED
#define HEATEQINCOMPNEWTON_HPP_INCLUDED

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
class HeatEqIncompNewton : public Equation
{
    public:
        HeatEqIncompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                          std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                          const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex);
        HeatEqIncompNewton(const HeatEqIncompNewton& equation)             = delete;
        HeatEqIncompNewton& operator=(const HeatEqIncompNewton& equation)  = delete;
        HeatEqIncompNewton(HeatEqIncompNewton&& equation)                  = delete;
        HeatEqIncompNewton& operator=(HeatEqIncompNewton&& equation)       = delete;
        ~HeatEqIncompNewton() override;

        void displayParams() const override;

        bool solve() override;

    private:
        enum class Res
        {
            T,
            Ax_f
        };

        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder; /**< Class responsible of building the required matrices. */
        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder2; /**< Class responsible of building the required matrices. */
        std::unique_ptr<PicardAlgo> m_pPicardAlgo;
        Eigen::SparseMatrix<double> m_A;
        Eigen::VectorXd m_b;
        EigenSparseSolver m_solver;
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> m_solverIt;

        Res m_residual;
        double m_k;
        double m_rho;
        double m_cv;
        double m_Tm;
        double m_Lm;
        double m_DT;
        double m_h;
        double m_Tinf;
        double m_epsRad;
        bool m_phaseChange;

        void m_buildAb(const Eigen::VectorXd& qPrev);
        void m_applyBC(const Eigen::VectorXd& qPrev);

        double m_getDflDT(double T);
};

#include "HeatEquation.inl"

#endif // HEATEQINCOMPNEWTON_HPP_INCLUDED
