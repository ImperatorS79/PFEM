#include "NonLinearAlgo.hpp"

#include <stdexcept>

NonLinearAlgo::NonLinearAlgo(std::function<void(const std::vector<Eigen::VectorXd>& /** qPrevVec **/)> prepare,
                             std::function<bool(std::vector<Eigen::VectorXd>& /** qIterVec **/,
                                                const std::vector<Eigen::VectorXd>& /** qPrevVec **/)> solve,
                             std::function<double(const std::vector<Eigen::VectorXd>& /** qIterVec **/,
                                                  const std::vector<Eigen::VectorXd>& /** qIterPrevVec **/)> computeRes):
m_prepare(prepare),
m_solve(solve),
m_computeRes(computeRes)
{

}

NonLinearAlgo::~NonLinearAlgo()
{

}

void NonLinearAlgo::displayParams()
{
    throw std::runtime_error("Unimplemented function by the child class -> NonLinearAlgo::displayParams()");
}

bool NonLinearAlgo::solve(Mesh* pMesh, const std::vector<Eigen::VectorXd>& qPrevVec, bool verboseOutput)
{
    (void)pMesh;
    (void)qPrevVec;
    (void)verboseOutput;
    throw std::runtime_error("Unimplemented function by the child class -> NonLinearAlgo::solve(pMesh, qPrev)");
}

