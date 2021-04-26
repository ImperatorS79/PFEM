#include "MomContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

void MomContEqIncompNewton::m_buildAbPSPG(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    m_clock.start();
    unsigned int dim = m_pMesh->getDim();
    unsigned int noPerEl = m_pMesh->getNodesPerElm();
    unsigned int tripletPerElm = (dim + 1)*noPerEl*(dim + 1)*noPerEl;
    unsigned int doubletPerElm = (dim + 1)*noPerEl;
    std::size_t nElm = m_pMesh->getElementsCount();
    std::size_t nNodes = m_pMesh->getNodesCount();
    double dt = m_pSolver->getTimeStep();

    std::vector<Eigen::Triplet<double>> indexA(tripletPerElm*nElm);
    std::vector<std::pair<std::size_t, double>> indexb(doubletPerElm*nElm); b.setZero();
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);
        Eigen::MatrixXd Ae((dim + 1)*noPerEl, (dim + 1)*noPerEl);
        Eigen::VectorXd be((dim + 1)*noPerEl);

        double tau = m_computeTauPSPG(element);
        Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
        Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);
        Eigen::MatrixXd Me_dt = (1/dt)*m_pMatBuilder->getM(element);
        Me_dt = MatrixBuilder::diagBlock(Me_dt, dim);
        Eigen::MatrixXd Ke = m_pMatBuilder->getK(element, Be);
        Eigen::MatrixXd De = m_pMatBuilder->getD(element, Be);
        Eigen::MatrixXd Ce_dt = (tau/dt)*m_pMatBuilder->getC(element, Be, gradNe);
        Eigen::MatrixXd Le = tau*m_pMatBuilder->getL(element, Be, gradNe);
        Eigen::VectorXd Fe = m_pMatBuilder->getF(element, m_bodyForce, Be);
        Eigen::VectorXd He = tau*m_pMatBuilder->getH(element, m_bodyForce, Be, gradNe);

        Ae << Me_dt + Ke, -De.transpose(), Ce_dt + De, Le;

        Eigen::VectorXd vPrev(noPerEl*dim);
        if(dim == 2)
            vPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)],
                     qPrev[element.getNodeIndex(0) + nNodes], qPrev[element.getNodeIndex(1) + nNodes], qPrev[element.getNodeIndex(2) + nNodes];
        else
            vPrev << qPrev[element.getNodeIndex(0)], qPrev[element.getNodeIndex(1)], qPrev[element.getNodeIndex(2)], qPrev[element.getNodeIndex(3)],
                     qPrev[element.getNodeIndex(0) + nNodes], qPrev[element.getNodeIndex(1) + nNodes], qPrev[element.getNodeIndex(2) + nNodes], qPrev[element.getNodeIndex(3) + nNodes],
                     qPrev[element.getNodeIndex(0) + 2*nNodes], qPrev[element.getNodeIndex(1) + 2*nNodes], qPrev[element.getNodeIndex(2) + 2*nNodes], qPrev[element.getNodeIndex(3) + 2*nNodes];

        be << Fe + Me_dt*vPrev, He + Ce_dt*vPrev;

        std::size_t countA = 0;
        std::size_t countb = 0;

        for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            const Node& ni = m_pMesh->getNode(element.getNodeIndex(i));

            for(unsigned short j = 0 ; j < noPerEl ; ++j)
            {
                for(unsigned short d1 = 0 ; d1 < dim ; ++d1)
                {
                    for(unsigned short d2 = 0 ; d2 <= dim ; ++d2)
                    {
                        if(!(ni.isBound() || ni.isFree()))
                        {
                            indexA[tripletPerElm*elm + countA] =
                                Eigen::Triplet<double>(element.getNodeIndex(i) + d1*nNodes,
                                                       element.getNodeIndex(j) + d2*nNodes,
                                                       Ae(i + d1*noPerEl, j + d2*noPerEl));
                        }
                        countA++;
                    }
                }

                for(unsigned short d2 = 0 ; d2 <= dim ; ++d2)
                {
                    if(!ni.isFree())
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i) + dim*nNodes,
                                                   element.getNodeIndex(j) + d2*nNodes,
                                                   Ae(i + dim*noPerEl, j + d2*noPerEl));
                    }

                    countA++;
                }
            }

            for(unsigned short d = 0 ; d <= dim ; ++d)
            {
                indexb[doubletPerElm*elm + countb] = std::make_pair(element.getNodeIndex(i) + d*nNodes, be(i + d*noPerEl));
                countb++;
            }
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    //Best would be to know the number of nodes in which case :/
    //This can still be fasten using OpenMP but will never be as good as using []
    //with preallocated memory
    m_clock.start();
    for(std::size_t n = 0 ; n < nNodes ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree())
        {
            indexA.push_back(Eigen::Triplet<double>(n + dim*nNodes,
                                                    n + dim*nNodes,
                                                    1));
        }

        if(node.isBound() || node.isFree())
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                indexA.push_back(Eigen::Triplet<double>(n + d*nNodes,
                                                        n + d*nNodes,
                                                        1));
            }
        }
    }
    m_accumalatedTimes["Push back (n, n, 1)"] += m_clock.end();

    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    m_clock.start();
    A.setFromTriplets(indexA.begin(), indexA.end());
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();

    m_clock.start();
    for(const auto& doublet : indexb)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        b[doublet.first] += doublet.second;
    }
    m_accumalatedTimes["Assemble vector"] += m_clock.end();
}

void MomContEqIncompNewton::m_applyBCPSPG(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::VectorXd& qPrev)
{
    const uint8_t dim = m_pMesh->getDim();
    const std::size_t nodesCount = m_pMesh->getNodesCount();

    const std::size_t facetsCount = m_pMesh->getFacetsCount();
    const std::size_t noPerFacet = m_pMesh->getNodesPerFacet();

    for(std::size_t f = 0 ; f < facetsCount ; ++f)
    {
        if(m_gamma < 1e-15)
            continue;

        const Facet& facet = m_pMesh->getFacet(f);

        bool onFS = true;
        for(unsigned short n = 0 ; n < noPerFacet ; ++n)
        {
            if(!facet.getNode(n).isOnFreeSurface())
            {
                onFS = false;
                break;
            }
        }
        if(!onFS)
            continue;

        Eigen::MatrixXd MGamma = m_pMatBuilder->getMGamma(facet);
        MGamma = MatrixBuilder::diagBlock(MGamma, m_pMesh->getDim());

        Eigen::VectorXd kappa_n(dim*noPerFacet);
        for(uint8_t n = 0 ; n < noPerFacet; ++n)
        {
            double curvature = m_pMesh->getFreeSurfaceCurvature(facet.getNodeIndex(n));
            std::array<double, 3> normal = m_pMesh->getBoundFSNormal(facet.getNodeIndex(n));

            for(uint8_t d = 0 ; d < dim ; ++d)
                kappa_n(n + d*noPerFacet) = curvature*normal[d];
        }

        Eigen::VectorXd Ff = MGamma*kappa_n;

        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
                b(facet.getNodeIndex(i) + d*nodesCount) += Ff[d*noPerFacet + i];
        }
    }

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree())
        {
            b(n + dim*nodesCount) = 0;

            if(!node.isBound())
            {
                for(uint8_t d = 0 ; d < dim  ; ++d)
                {
                    b(n + d*nodesCount) = qPrev(n + d*nodesCount) + m_pSolver->getTimeStep()*m_bodyForce[d];
                }
            }
        }

        if(node.isBound())
        {
            if(node.getFlag(m_bcFlags[0]))
            {
                std::vector<double> result;
                result = m_bcParams[0].call<std::vector<double>>(m_pMesh->getNodeType(n) + "V",
                                                        node.getPosition(),
                                                        m_pMesh->getBoundNodeInitPos(n),
                                                        m_pProblem->getCurrentSimTime() +
                                                        m_pSolver->getTimeStep());

                for(uint8_t d = 0 ; d < dim ; ++d)
                {
                    b(n + d*nodesCount) = result[d];
                    for(Eigen::SparseMatrix<double>::InnerIterator it(A, n + d*nodesCount); it; ++it)
                    {
                        Eigen::Index row = it.row();
                        if(row == it.col())
                            continue;

                        double value = it.value();
                        b(row) -= value*result[d];
                        it.valueRef() = 0;
                    }
                }
            }
        }
    }

    A.makeCompressed();
}

double MomContEqIncompNewton::m_computeTauPSPG(const Element& element) const
{
    const double h = std::sqrt(m_pMesh->getRefElementSize(m_pMesh->getDim())*element.getDetJ()/M_PI);

    double U = 0;
    for (unsigned short n = 0 ; n < m_pMesh->getNodesPerElm() ; ++n)
    {
        const Node& node = m_pMesh->getNode(element.getNodeIndex(n));

        double nodeU = 0;
        for (unsigned short d = 0 ; d < m_pMesh->getDim() ; ++d)
        {
            nodeU += node.getState(d)*node.getState(d);
        }
        U += std::sqrt(nodeU);
    }
    U /= (m_pMesh->getDim() + 1);

    return 1/std::sqrt((2/m_pSolver->getTimeStep())*(2/m_pSolver->getTimeStep()) + (2*U/h)*(2*U/h)
                        + 9*(4*m_mu/(h*h*m_rho))*(4*m_mu/(h*h*m_rho)));
}

void MomContEqIncompNewton::m_setupPicardPSPG(unsigned int maxIter, double minRes)
{
    m_pPicardAlgo = std::make_unique<PicardAlgo>([&](const auto& qPrevVec){
        m_clock.start();
        m_A.resize(qPrevVec[0].rows(), qPrevVec[0].rows());
        m_b.resize(qPrevVec[0].rows()); m_b.setZero();
        m_accumalatedTimes["Prepare Picard algorithm"] += m_clock.end();
        m_clock.start();
        m_pMesh->saveNodesList();
        m_accumalatedTimes["Save/restore nodelist"] += m_clock.end();
    },
    [&](auto& qIterVec, const auto& qPrevVec){
        m_buildAbPSPG(m_A, m_b, qPrevVec[0]);
        m_clock.start();
        m_applyBCPSPG(m_A, m_b, qPrevVec[0]);
        m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();

        m_clock.start();
        m_solver.compute(m_A);
        m_accumalatedTimes["Compute A matrix"] += m_clock.end();

        if(m_solver.info() == Eigen::Success)
        {
            m_clock.start();
            qIterVec[0] = m_solver.solve(m_b);
            m_accumalatedTimes["Solve system"] += m_clock.end();
            m_clock.start();
            setNodesStatesfromQ(m_pMesh, qIterVec[0], m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim());
            Eigen::VectorXd deltaPos = qIterVec[0]*m_pSolver->getTimeStep();
            m_pMesh->updateNodesPositionFromSave(deltaPos);
            m_accumalatedTimes["Update solutions"] += m_clock.end();
            return true;
        }
        else
        {
            if(m_pProblem->isOutputVerbose())
                std::cout << "\t * The Eigen::SparseLU solver failed to factorize the A matrix!" << std::endl;
            m_clock.start();
            m_pMesh->restoreNodesList();
            m_accumalatedTimes["Save/restore nodelist"] += m_clock.end();
            return false;
        }
    },
    [&](const auto& qIterVec, const auto& qIterPrevVec) -> double {
        m_clock.start();
        double num = 0, den = 0;
        Mesh* p_Mesh = this->m_pMesh;

        for(std::size_t n = 0 ; n < p_Mesh->getNodesCount() ; ++n)
        {
            const Node& node = p_Mesh->getNode(n);

            if(!node.isFree())
            {
                for(unsigned short d = 0 ; d < p_Mesh->getDim() ; ++d)
                {
                    num += (qIterVec[0](n + d*p_Mesh->getNodesCount()) - qIterPrevVec[0](n + d*p_Mesh->getNodesCount()))*(qIterVec[0](n + d*p_Mesh->getNodesCount()) - qIterPrevVec[0](n + d*p_Mesh->getNodesCount()));
                    den += qIterPrevVec[0](n + d*p_Mesh->getNodesCount())*qIterPrevVec[0](n + d*p_Mesh->getNodesCount());
                }
            }
        }

        double res;
        if(den == 0)
            res = std::numeric_limits<double>::max();
        else
            res = std::sqrt(num/den);
        m_accumalatedTimes["Compute Picard Algo residual"] += m_clock.end();
        return res;
    }, maxIter, minRes);
}
