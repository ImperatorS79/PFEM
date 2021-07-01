#include "HeatEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

template<unsigned short dim>
HeatEqIncompNewton<dim>::HeatEqIncompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "HeatEq")
{
    unsigned int nGPHD = 0;
    unsigned int nGPLD = 0;
    if constexpr (dim == 2)
    {
        nGPHD = 3;
        nGPLD = 3;
    }
    else if constexpr (dim == 3)
    {
        nGPHD = 4;
        nGPLD = 3;
    }
    m_pMatBuilder = std::make_unique<MatrixBuilder<dim>>(*pMesh, nGPHD, nGPLD);
    m_pMatBuilder2 = std::make_unique<MatrixBuilder<dim>>(*pMesh, nGPHD, nGPLD);

    m_k = m_materialParams[0].checkAndGet<double>("k");
    m_rho = m_materialParams[0].checkAndGet<double>("rho");
    m_cv = m_materialParams[0].checkAndGet<double>("cv");
    m_h = m_materialParams[0].checkAndGet<double>("h");
    m_Tinf = m_materialParams[0].checkAndGet<double>("Tinf");
    m_epsRad = m_materialParams[0].checkAndGet<double>("epsRad");

    m_phaseChange = false;
    if(m_materialParams[0].doesVarExist("Tm")  || m_materialParams[0].doesVarExist("C") ||
       m_materialParams[0].doesVarExist("eps") || m_materialParams[0].doesVarExist("DT") ||
       m_materialParams[0].doesVarExist("Lm"))
    {
        m_phaseChange = true;

        m_Tm = m_materialParams[0].checkAndGet<double>("Tm");
        m_DT = m_materialParams[0].checkAndGet<double>("DT");
        m_Lm = m_materialParams[0].checkAndGet<double>("Lm");
    }

    if(bcFlags.size() != 2)
        throw std::runtime_error("the " + getID() + " equation require two flags for two possible boundary conditiosn!");

    if(statesIndex.size() != 1)
        throw std::runtime_error("the " + getID() + " equation require one state index describing the T state !");

    if(!m_phaseChange)
    {
        m_pMatBuilder->setMcomputeFactor([&](const Element& /** element **/,
                                             const NmatTypeHD<dim>& /** N **/) -> double {
            return m_cv*m_rho;
        });
    }
    else
    {
        m_pMatBuilder->setMcomputeFactor([&](const Element& element,
                                             const NmatTypeHD<dim>& N) -> double {
            double T = (N*getElementState<dim>(element, m_statesIndex[0])).value();
            return (m_cv + m_getDflDT(T)*m_Lm)*m_rho;
        });
    }

    m_pMatBuilder->setLcomputeFactor([&](const Element& /** element **/,
                                         const NmatTypeHD<dim>& /** N **/,
                                         const BmatType<dim>& /** B **/) -> double {
        return m_k;
    });

    m_pMatBuilder->setQFunc([&](const Facet& facet, const std::array<double, 3>& gp) -> Eigen::Matrix<double, dim, 1> {
        std::array<double, 3> pos = facet.getPosFromGP(gp);

        std::string nodeType = facet.isOnFreeSurface() ? "FreeSurface" : m_pMesh->getNodeType(facet.getNodeIndex(0));

        std::array<double, dim> result;
            result = m_bcParams[0].call<std::array<double, dim>>(nodeType + "Q",
                                                             pos,
                                                             m_pProblem->getCurrentSimTime() +
                                                             m_pSolver->getTimeStep());

        return Eigen::Matrix<double, dim, 1>::Map(result.data(), result.size());
    });

    m_pMatBuilder->setSGammacomputeFactor([&](const Facet& facet,
                                              const NmatTypeLD<dim>& N) -> double {
        double T = (N*getFacetState<dim>(facet, m_statesIndex[0])).value();
        return m_h*(T - m_Tinf);
    });

    m_pMatBuilder2->setSGammacomputeFactor([&](const Facet& facet,
                                               const NmatTypeLD<dim>& N) -> double {
        double T = (N*getFacetState<dim>(facet, m_statesIndex[0])).value();
        constexpr double sigma = 5.670374419e-8;
        return sigma*m_epsRad*(T*T*T*T - m_Tinf*m_Tinf*m_Tinf*m_Tinf);
    });

    unsigned int maxIter = m_equationParams[0].checkAndGet<unsigned int>("maxIter");
    double minRes = m_equationParams[0].checkAndGet<double>("minRes");
    std::string residual = m_equationParams[0].checkAndGet<std::string>("residual");
    if(residual == "T")
        m_residual = Res::T;
    else if(residual == "Ax_f")
        m_residual = Res::Ax_f;
    else
        throw std::runtime_error("unknown residual type: " + residual);

    m_pPicardAlgo = std::make_unique<PicardAlgo>([&](const auto& qPrevVec){
        m_clock.start();
        m_A.resize(qPrevVec[0].rows(), qPrevVec[0].rows());
        m_b.resize(qPrevVec[0].rows()); m_b.setZero();
        m_accumalatedTimes["Prepare Picard Algo"] += m_clock.end();

        m_clock.start();
        m_pMesh->saveNodesList();
        m_accumalatedTimes["Save/restore nodeslist"] += m_clock.end();

        m_buildAb(qPrevVec[0]);

        m_clock.start();
        m_applyBC(qPrevVec[0]);
        m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();
    },
    [&](auto& qIterVec, const auto& qPrevVec){
        m_clock.start();
        m_solverIt.compute(m_A);
        m_accumalatedTimes["Compute matrix"] += m_clock.end();

        if(m_solverIt.info() == Eigen::Success)
        {
            m_clock.start();
            qIterVec[0] = m_solverIt.solveWithGuess(m_b, qPrevVec[0]);
            m_accumalatedTimes["Solve system"] += m_clock.end();
            m_clock.start();
            setNodesStatesfromQ(m_pMesh, qIterVec[0], m_statesIndex[0], m_statesIndex[0]);
            m_accumalatedTimes["Update solution"] += m_clock.end();
            m_buildAb(qPrevVec[0]);

            m_clock.start();
            m_applyBC(qPrevVec[0]);
            m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();
            return true;
        }
        else
        {
            if(m_pProblem->isOutputVerbose())
                std::cout << "\t * The Eigen::SparseLU solver failed to factorize the A matrix!" << std::endl;
            m_clock.start();
            m_pMesh->restoreNodesList();
            m_accumalatedTimes["Save/restore nodeslist"] += m_clock.end();
            return false;
        }
    },
    [&](const auto& qIterVec, const auto& qIterPrevVec) -> double {
        m_clock.start();
        Mesh* p_Mesh = this->m_pMesh;

        if(m_residual == Res::T)
        {
            double num = 0, den = 0;
            for(std::size_t n = 0 ; n < p_Mesh->getNodesCount() ; ++n)
            {
                const Node& node = p_Mesh->getNode(n);

                if(!node.isFree())
                {
                    num += (qIterVec[0](n) - qIterPrevVec[0](n))*(qIterVec[0](n) - qIterPrevVec[0](n));
                    den += qIterPrevVec[0](n)*qIterPrevVec[0](n);
                }
            }

            double res;
            if(den == 0)
                res = std::numeric_limits<double>::max();
            else
                res = std::sqrt(num/den);
            m_accumalatedTimes["Compute Picard Algo residual"] += m_clock.end();
            return res;
        }
        else
        {
            return(m_A*qIterVec[0] - m_b).norm();
        }

    }, maxIter, minRes);

    if(!m_phaseChange)
        m_pPicardAlgo->runOnlyOnce(true); //When k, cv independant of T, no need of Picard

    m_needNormalCurv = true;
}

template<unsigned short dim>
HeatEqIncompNewton<dim>::~HeatEqIncompNewton()
{

}

template<unsigned short dim>
void HeatEqIncompNewton<dim>::displayParams() const
{
    std::cout << "Heat equation parameters:\n"
              << " * Specific heat capacity: " << m_cv << " J/(kg K)\n"
              << " * Heat conduction: " << m_k << " W/(mK)" << std::endl;

    m_pPicardAlgo->displayParams();
}

template<unsigned short dim>
bool HeatEqIncompNewton<dim>::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Heat Equation" << std::endl;

    std::vector<Eigen::VectorXd> qPrevVec = {getQFromNodesStates(m_pMesh, m_statesIndex[0], m_statesIndex[0])};
    return m_pPicardAlgo->solve(m_pMesh, qPrevVec, m_pProblem->isOutputVerbose());
}

template<unsigned short dim>
void HeatEqIncompNewton<dim>::m_buildAb(const Eigen::VectorXd& qPrev)
{
    constexpr unsigned short noPerEl = dim + 1;
    const unsigned int tripletPerElm = (dim*noPerEl*noPerEl + dim*noPerEl*dim*noPerEl + 3*noPerEl*dim*noPerEl + noPerEl*noPerEl);
    const unsigned int doubletPerElm = (2*dim*noPerEl + 2*noPerEl);
    const std::size_t nElm = m_pMesh->getElementsCount();
    const std::size_t nNodes = m_pMesh->getNodesCount();
    const double dt = m_pSolver->getTimeStep();

    std::vector<Eigen::Triplet<double>> indexA(tripletPerElm*nElm);
    std::vector<std::pair<std::size_t, double>> indexb(doubletPerElm*nElm); m_b.setZero();

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        auto gradNe = m_pMatBuilder->getGradN(element);
        auto Be = m_pMatBuilder->getB(gradNe);
        auto Me_dt = (1/dt)*m_pMatBuilder->getM(element);
        auto Le = m_pMatBuilder->getL(element, Be, gradNe);

        auto thetaPrev = getElementState<dim>(qPrev, element, 0, nNodes);

        auto MThetaPreve_dt = Me_dt*thetaPrev;

        std::size_t countA = 0;
        std::size_t countb = 0;

         for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            const Node& ni = element.getNode(i);

            for(unsigned short j = 0 ; j < noPerEl ; ++j)
            {
                if(!m_pSolver->getBcTagFlags(ni.getTag(), m_bcFlags[0]))
                {
                    if(!ni.isFree())
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i),
                                                   element.getNodeIndex(j),
                                                   Me_dt(i, j));
                    }
                    countA++;

                    if(!ni.isFree())
                    {
                        indexA[tripletPerElm*elm + countA] =
                            Eigen::Triplet<double>(element.getNodeIndex(i),
                                                   element.getNodeIndex(j),
                                                   Le(i, j));
                    }
                    countA++;
                }
            }

            /************************************************************************
                                            Build h
            ************************************************************************/
            indexb[noPerEl*elm + countb] = std::make_pair(element.getNodeIndex(i), MThetaPreve_dt[i]);

            countb++;
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());

    //Best would be to know the number of nodes in which case :/
    //This can still be fasten using OpenMP but will never be as good as using []
    //with preallocated memory
    for(std::size_t n = 0 ; n < nNodes ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(m_pSolver->getBcTagFlags(node.getTag(), m_bcFlags[0]) || node.isFree())
        {
            indexA.push_back(Eigen::Triplet<double>(n, n, 1));
        }
    }

    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    m_A.setFromTriplets(indexA.begin(), indexA.end());

    for(const auto& doublet : indexb)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        m_b[doublet.first] += doublet.second;
    }
}

template<unsigned short dim>
void HeatEqIncompNewton<dim>::m_applyBC(const Eigen::VectorXd& qPrev)
{
    const std::size_t nodesCount = m_pMesh->getNodesCount();
    const std::size_t facetsCount = m_pMesh->getFacetsCount();
    constexpr std::size_t noPerFacet = dim;

    for(std::size_t f = 0 ; f < facetsCount ; ++f)
    {
        const Facet& facet = m_pMesh->getFacet(f);
        bool boundaryQ = true;
        bool boundaryQh = true;
        bool boundaryQr = true;
        bool fs = facet.isOnFreeSurface();
        for(uint8_t n = 0 ; n < noPerFacet ; ++n)
        {
            const Node& node = facet.getNode(n);

            int tag = node.getTag();
            if(fs)
                tag = -2;

            if(!m_pSolver->getBcTagFlags(tag, m_bcFlags[1]))
                boundaryQ = false;

            if(!m_pSolver->getBcTagFlags(tag, m_bcFlags[2]))
                boundaryQh = false;

            if(!m_pSolver->getBcTagFlags(tag, m_bcFlags[3]))
                boundaryQr = false;
        }

        if(!boundaryQ && !boundaryQh && !boundaryQr)
            continue;

        if(boundaryQ)
        {
            auto qn = m_pMatBuilder->getQN(facet);

            for(unsigned short i = 0 ; i < noPerFacet ; ++i)
            {
                m_b(facet.getNodeIndex(i)) -= qn[i];
            }
        }

        if(boundaryQh)
        {
            auto SGammah = m_pMatBuilder->getSGamma(facet);

            for(unsigned short i = 0 ; i < noPerFacet ; ++i)
            {
                m_b(facet.getNodeIndex(i)) -= SGammah[i];
            }
        }

        if(boundaryQr)
        {
            auto SGammar = m_pMatBuilder2->getSGamma(facet);

            for(unsigned short i = 0 ; i < noPerFacet ; ++i)
            {
                m_b(facet.getNodeIndex(i)) -= SGammar[i];
            }
        }

    }

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree() && !m_pSolver->getBcTagFlags(node.getTag(), m_bcFlags[0]))
        {
            m_b(n) = qPrev[n];
        }
        else if(m_pSolver->getBcTagFlags(node.getTag(), m_bcFlags[0]))
        {
            std::array<double, 1> result;
            result = m_bcParams[0].call<std::array<double, 1>>(m_pMesh->getNodeType(n) + "T",
                                                             node.getPosition(),
                                                             m_pProblem->getCurrentSimTime() +
                                                             m_pSolver->getTimeStep());
            m_b(n) = result[0];
            for(Eigen::SparseMatrix<double>::InnerIterator it(m_A, n); it; ++it)
            {
                Eigen::Index row = it.row();
                if(row == it.col())
                    continue;

                double value = it.value();
                m_b(row) -= value*result[0];
                it.valueRef() = 0;
            }
        }
    }

    m_A.makeCompressed();
}

template<unsigned short dim>
double HeatEqIncompNewton<dim>::m_getDflDT(double T)
{
    return static_cast<double>(T > m_Tm - m_DT/2)*static_cast<double>(T < m_Tm + m_DT/2)*(-2*3/(m_DT*m_DT*m_DT)*T*T + 6*2*m_Tm/(m_DT*m_DT*m_DT)*T + (3/(2*m_DT) - 6*m_Tm*m_Tm/(m_DT*m_DT*m_DT)));
}
