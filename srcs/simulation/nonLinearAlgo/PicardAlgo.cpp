#include "PicardAlgo.hpp"

#include <iostream>

#include "../../mesh/Mesh.hpp"

PicardAlgo::PicardAlgo(std::function<void(const std::vector<Eigen::VectorXd>& /** qPrevVec **/)> prepare,
                       std::function<bool(std::vector<Eigen::VectorXd>& /** qIterVec **/,
                                          const std::vector<Eigen::VectorXd>& /** qPrevVec **/)> solve,
                       std::function<double(const std::vector<Eigen::VectorXd>& /** qIterVec **/,
                                            const std::vector<Eigen::VectorXd>& /** qIterPrevVec **/)> computeRes,
                       unsigned int maxIter, double minRes):
NonLinearAlgo(prepare, solve, computeRes),
m_maxIter(maxIter),
m_minRes(minRes)
{

}

PicardAlgo::~PicardAlgo()
{

}

void PicardAlgo::displayParams()
{
    std::cout << " * Maximum residual: " << m_minRes << "\n"
              << " * Maximum iteration count: " << m_maxIter << std::endl;
}

bool PicardAlgo::solve(Mesh* pMesh, const std::vector<Eigen::VectorXd>& qPrevVec, bool verboseOutput)
{
    m_prepare(qPrevVec);

    std::vector<Eigen::VectorXd> qIterVec(qPrevVec.size());
    std::vector<Eigen::VectorXd> qIterPrevVec(qPrevVec.size());

	for(std::size_t i = 0 ; i < qIterVec.size() ; ++i)
	{
		qIterVec[i].resize(qPrevVec[i].rows()); qIterVec[i].setZero();
		qIterPrevVec[i] = qPrevVec[i];
	}

	unsigned int iterCount = 0;
    double res = std::numeric_limits<double>::max();

    while(res > m_minRes)
    {
        if(verboseOutput && !m_runOnce)
        {
            std::cout << " - Picard algorithm (mesh position) - iteration ("
                      << iterCount << ")" << std::endl;
        }

        if(iterCount > m_maxIter)
        {
            if(verboseOutput && !m_runOnce)
            {
                std::cout << "\t * Iteration count " << iterCount
                      << " greater than maximum: " << m_maxIter << std::endl;
            }

            pMesh->restoreNodesList();
            return false;
        }

        if(!m_solve(qIterVec, qPrevVec))
            return false;

        res = m_computeRes(qIterVec, qIterPrevVec);

        if(verboseOutput && !m_runOnce)
            std::cout << "\t * Relative 2-norm of q: " << res << " vs "
                      << m_minRes << std::endl;

        for(std::size_t i = 0 ; i < qIterVec.size() ; ++i)
			qIterPrevVec[i] = qIterVec[i];

        if(std::isnan(res))
        {
            if(verboseOutput)
                std::cout << "\t * NaN residual" << std::endl;

            pMesh->restoreNodesList();
            return false;
        }

        iterCount++;

        if(m_runOnce)
            break;
    }

    return true;
}
