#pragma once
#ifndef NONLINEARALGO_HPP_INCLUDED
#define NONLINEARALGO_HPP_INCLUDED

#include <functional>
#include <Eigen/Sparse>
#include <Eigen/Dense>

class Mesh;

/**
 * \class NonLinearAlgo
 * \brief Represents a non-linear algorithm to solve a Ax = b system of equations.
 */
class NonLinearAlgo
{
    public:
        NonLinearAlgo(std::function<void(const std::vector<Eigen::VectorXd>& /** qPrevVec **/)> prepare,
                      std::function<bool(std::vector<Eigen::VectorXd>& /** qIterVec **/,
                                         const std::vector<Eigen::VectorXd>& /** qPrevVec **/)> solve,
                      std::function<double(const std::vector<Eigen::VectorXd>& /** qIterVec **/,
                                           const std::vector<Eigen::VectorXd>& /** qIterPrevVec **/)> computeRes);

        virtual ~NonLinearAlgo();

        virtual void displayParams();
        virtual bool solve(Mesh* pMesh, const std::vector<Eigen::VectorXd>& qPrevVec, bool verboseOutput);
        void runOnlyOnce(bool runOnce)
        {
            m_runOnce = runOnce;
        }

    protected:
        bool m_runOnce = false;
        std::function<void(const std::vector<Eigen::VectorXd>&)> m_prepare;

        std::function<bool(std::vector<Eigen::VectorXd>&,
                           const std::vector<Eigen::VectorXd>&)> m_solve;

        std::function<double(const std::vector<Eigen::VectorXd>&,
                             const std::vector<Eigen::VectorXd>&)> m_computeRes;
};

#endif // NONLINEARALGO_HPP_INCLUDED
