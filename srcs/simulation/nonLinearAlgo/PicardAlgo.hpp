#pragma once
#ifndef PICARDALGO_HPP_INCLUDED
#define PICARDALGO_HPP_INCLUDED

#include <functional>

#include "NonLinearAlgo.hpp"


/**
 * \class PicardAlgo
 * \brief Represents a Picard non-linear algorithm to solve a Ax = b system of equations.
 */
class PicardAlgo : public NonLinearAlgo
{
    public:
        PicardAlgo(std::function<void(const std::vector<Eigen::VectorXd>& /** qPrevVec **/)> prepare,
                   std::function<bool(std::vector<Eigen::VectorXd>& /** qIterVec **/,
                                      const std::vector<Eigen::VectorXd>& /** qPrevVec **/)> solve,
                   std::function<double(const std::vector<Eigen::VectorXd>& /** qIterVec **/,
                                        const std::vector<Eigen::VectorXd>& /** qIterPrevVec **/)> computeRes,
                   unsigned int maxIter, double minRes);

        ~PicardAlgo() override;

        void displayParams() override;
        bool solve(Mesh* pMesh, const std::vector<Eigen::VectorXd>& qPrevVec, bool verboseOutput) override;
    private:
        unsigned int m_maxIter;
        double m_minRes;
};

#endif // PICARDALGO_HPP_INCLUDED
