#include "MomContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_buildAbPSPG(const Eigen::VectorXd& qPrev)
{
    m_clock.start();
    constexpr unsigned short nodPerEl = dim + 1;
    const unsigned int tripletPerElm = (dim + 1)*nodPerEl*(dim + 1)*nodPerEl;
    const unsigned int doubletPerElm = (dim + 1)*nodPerEl;
    const std::size_t nElm = m_pMesh->getElementsCount();
    const std::size_t nNodes = m_pMesh->getNodesCount();
    const double dt = m_pSolver->getTimeStep();

    std::vector<Eigen::Triplet<double>> indexA(tripletPerElm*nElm);
    std::vector<std::pair<std::size_t, double>> indexb(doubletPerElm*nElm); m_b.setZero();
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);
        Eigen::Matrix<double, (dim + 1)*nodPerEl, (dim + 1)*nodPerEl> Ae;
        Eigen::Matrix<double, (dim + 1)*nodPerEl, 1> be;

        double tau = m_computeTauPSPG(element);
        GradNmatType<dim> gradNe = m_pMatBuilder->getGradN(element);
        BmatType<dim> Be = m_pMatBuilder->getB(gradNe);
        Eigen::Matrix<double, nodPerEl, nodPerEl> Me_dt_s = (1/dt)*m_pMatBuilder->getM(element);
        Eigen::Matrix<double, dim*nodPerEl, dim*nodPerEl> Me_dt = MatrixBuilder<dim>::diagBlock(Me_dt_s);
        Eigen::Matrix<double, dim*nodPerEl, dim*nodPerEl> Ke = m_pMatBuilder->getK(element, Be);
        Eigen::Matrix<double, nodPerEl, dim*nodPerEl> De = m_pMatBuilder->getD(element, Be);
        Eigen::Matrix<double, nodPerEl, dim*nodPerEl> Ce_dt = (tau/dt)*m_pMatBuilder->getC(element, Be, gradNe);
        Eigen::Matrix<double, nodPerEl, nodPerEl> Le = tau*m_pMatBuilder->getL(element, Be, gradNe);
        Eigen::Matrix<double, dim*nodPerEl, 1> Fe = m_pMatBuilder->getF(element, m_bodyForce, Be);
        Eigen::Matrix<double, nodPerEl, 1> He = tau*m_pMatBuilder->getH(element, m_bodyForce, Be, gradNe);

        Ae << Me_dt + Ke, -De.transpose(), Ce_dt + De, Le;

        Eigen::Matrix<double, dim*nodPerEl, 1> vPrev = getElementVecState<dim>(qPrev, element, 0, nNodes);

        be << Fe + Me_dt*vPrev, He + Ce_dt*vPrev;

        std::size_t countA = 0;
        std::size_t countb = 0;

        for(unsigned short i = 0 ; i < nodPerEl ; ++i)
        {
            const Node& ni = m_pMesh->getNode(element.getNodeIndex(i));

            for(unsigned short j = 0 ; j < nodPerEl ; ++j)
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
                                                       Ae(i + d1*nodPerEl, j + d2*nodPerEl));
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
                                                   Ae(i + dim*nodPerEl, j + d2*nodPerEl));
                    }

                    countA++;
                }
            }

            for(unsigned short d = 0 ; d <= dim ; ++d)
            {
                indexb[doubletPerElm*elm + countb] = std::make_pair(element.getNodeIndex(i) + d*nNodes, be(i + d*nodPerEl));
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
    m_A.setFromTriplets(indexA.begin(), indexA.end());
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();

    m_clock.start();
    for(const auto& doublet : indexb)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        m_b[doublet.first] += doublet.second;
    }
    m_accumalatedTimes["Assemble vector"] += m_clock.end();
}

template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_applyBCPSPG(const Eigen::VectorXd& qPrev)
{
    const std::size_t nodesCount = m_pMesh->getNodesCount();
    const std::size_t facetsCount = m_pMesh->getFacetsCount();
    constexpr unsigned short noPerFacet = dim;

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

        auto MGamma_s = m_pMatBuilder->getMGamma(facet);
        auto MGamma = MatrixBuilder<dim>::diagBlock(MGamma_s);

        Eigen::Matrix<double, dim*noPerFacet, 1> nVec;
        for(uint8_t n = 0 ; n < noPerFacet; ++n)
        {
            std::array<double, 3> normal = m_pMesh->getBoundFSNormal(facet.getNodeIndex(n));

            for(uint8_t d = 0 ; d < dim ; ++d)
                nVec(n + d*noPerFacet) = normal[d];
        }

        Eigen::Matrix<double, dim*noPerFacet, 1> Ff = MGamma*nVec;

        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
                m_b(facet.getNodeIndex(i) + d*nodesCount) += Ff[d*noPerFacet + i];
        }
    }

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree())
        {
            m_b(n + dim*nodesCount) = 0;

            if(!node.isBound())
            {
                for(uint8_t d = 0 ; d < dim  ; ++d)
                {
                    m_b(n + d*nodesCount) = qPrev(n + d*nodesCount) + m_pSolver->getTimeStep()*m_bodyForce[d];
                }
            }
        }

        if(node.isBound())
        {
            if(node.getFlag(m_bcFlags[0]))
            {
                std::array<double, dim> result;
                result = m_bcParams[0].call<std::array<double, dim>>(m_pMesh->getNodeType(n) + "V",
                                                        node.getPosition(),
                                                        m_pMesh->getBoundNodeInitPos(n),
                                                        m_pProblem->getCurrentSimTime() +
                                                        m_pSolver->getTimeStep());

                for(uint8_t d = 0 ; d < dim ; ++d)
                {
                    m_b(n + d*nodesCount) = result[d];
                    for(Eigen::SparseMatrix<double>::InnerIterator it(m_A, n + d*nodesCount); it; ++it)
                    {
                        Eigen::Index row = it.row();
                        if(row == it.col())
                            continue;

                        double value = it.value();
                        m_b(row) -= value*result[d];
                        it.valueRef() = 0;
                    }
                }
            }
        }
    }

    m_A.makeCompressed();
}

template<unsigned short dim>
double MomContEqIncompNewton<dim>::m_computeTauPSPG(const Element& element) const
{
    const double h = std::sqrt(m_pMesh->getRefElementSize(dim)*element.getDetJ()/M_PI);
    constexpr unsigned short nodPerEl = dim + 1;

    double U = 0;
    for (unsigned short n = 0 ; n < nodPerEl ; ++n)
    {
        const Node& node = m_pMesh->getNode(element.getNodeIndex(n));

        double nodeU = 0;
        for (unsigned short d = 0 ; d < dim ; ++d)
        {
            nodeU += node.getState(d)*node.getState(d);
        }
        U += std::sqrt(nodeU);
    }
    U /= (m_pMesh->getDim() + 1);

    return 1/std::sqrt((2/m_pSolver->getTimeStep())*(2/m_pSolver->getTimeStep()) + (2*U/h)*(2*U/h)
                        + 9*(4*m_mu/(h*h*m_rho))*(4*m_mu/(h*h*m_rho)));
}

template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_setupPicardPSPG(unsigned int maxIter, double minRes)
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
        m_buildAbPSPG(qPrevVec[0]);
        m_clock.start();
        m_applyBCPSPG(qPrevVec[0]);
        m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();

        m_clock.start();
        m_solver.analyzePattern(m_A);
        m_accumalatedTimes["Analyse pattern of A matrix"] += m_clock.end();
        m_clock.start();
        m_solver.factorize(m_A);
        m_accumalatedTimes["Factorize A matrix"] += m_clock.end();

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
        Mesh* p_Mesh = this->m_pMesh;
        const std::size_t nNodes = p_Mesh->getNodesCount();
        double resV = 0, resP = 0;

        double num = 0, den = 0;
        for(std::size_t n = 0 ; n < nNodes ; ++n)
        {
            const Node& node = p_Mesh->getNode(n);

            if(!node.isFree())
            {
                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    num += (qIterVec[0](n + d*nNodes) - qIterPrevVec[0](n + d*nNodes))*(qIterVec[0](n + d*nNodes) - qIterPrevVec[0](n + d*nNodes));
                    den += qIterPrevVec[0](n + d*nNodes)*qIterPrevVec[0](n + d*nNodes);
                }
            }
        }

        if(den == 0)
            resV = std::numeric_limits<double>::max();
        else
            resV = std::sqrt(num/den);

        if(m_computePres)
        {
            num = 0, den = 0;
            for(std::size_t n = 0 ; n < nNodes ; ++n)
            {
                const Node& node = p_Mesh->getNode(n);

                if(!node.isFree())
                {
                    num += (qIterVec[0](n + dim*nNodes) - qIterPrevVec[0](n + dim*nNodes))*(qIterVec[0](n + dim*nNodes) - qIterPrevVec[0](n + dim*nNodes));
                    den += qIterPrevVec[0](n + dim*nNodes)*qIterPrevVec[0](n + dim*nNodes);
                }
            }

            if(den == 0)
                resP = std::numeric_limits<double>::max();
            else
                resP = std::sqrt(num/den);
        }

        m_accumalatedTimes["Compute Picard Algo residual"] += m_clock.end();

        return std::max(resV, resP);
    }, maxIter, minRes);
}
