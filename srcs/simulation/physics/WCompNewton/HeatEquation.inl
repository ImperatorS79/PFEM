#include "HeatEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

template<unsigned short dim>
HeatEqWCompNewton<dim>::HeatEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
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

    if(bcFlags.size() != 4)
        throw std::runtime_error("the " + getID() + " equation require four flags for four possible boundary conditions!");

    if(statesIndex.size() != 2)
        throw std::runtime_error("the " + getID() + " equation require two states index describing the T and rho state !");


    if(!m_phaseChange)
    {
        m_pMatBuilder->setMcomputeFactor([&](const Element& element,
                                             const NmatTypeHD<dim>& N) -> double {
            double rho = (N*getElementState<dim>(element, m_statesIndex[1])).value();
            return m_cv*rho;
        });
    }
    else
    {
        m_pMatBuilder->setMcomputeFactor([&](const Element& element,
                                             const NmatTypeHD<dim>& N) -> double {
            double rho = (N*getElementState<dim>(element, m_statesIndex[1])).value();
            double T = (N*getElementState<dim>(element, m_statesIndex[0])).value();
            return (m_cv + m_getDflDT(T)*m_Lm)*rho;
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
                                                             m_pProblem->getCurrentSimTime());

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

    m_needNormalCurv = true;
}

template<unsigned short dim>
HeatEqWCompNewton<dim>::~HeatEqWCompNewton()
{

}

template<unsigned short dim>
void HeatEqWCompNewton<dim>::displayParams() const
{
    std::cout << "Heat equation parameters:\n"
              << " * Specific heat capacity: " << m_cv << " J/(kg K)\n"
              << " * Heat conduction: " << m_k << " W/(mK)" << std::endl;
    if(m_phaseChange)
    {
        std::cout << " * Melting temperature: " << m_Tm << " K\n"
                  << " * DT: " << m_DT << " K" << std::endl;
    }
}

template<unsigned short dim>
double HeatEqWCompNewton<dim>::getDiffusionParam(const Node& node) const
{
    return m_k/(m_cv*node.getState(m_statesIndex[1]));
}

template<unsigned short dim>
bool HeatEqWCompNewton<dim>::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Heat Equation" << std::endl;

    m_buildSystem();

    m_clock.start();
    m_applyBC();
    m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();

    m_clock.start();
    Eigen::VectorXd qT = m_invM*m_F;
    m_accumalatedTimes["Solve system"] += m_clock.end();

    m_clock.start();
    setNodesStatesfromQ(m_pMesh, qT, m_statesIndex[0], m_statesIndex[0]);
    m_accumalatedTimes["Update solutions"] += m_clock.end();

    return true;
}

template<unsigned short dim>
void HeatEqWCompNewton<dim>::m_applyBC()
{
    auto& invMDiag = m_invM.diagonal();

    const std::size_t nodesCount = m_pMesh->getNodesCount();
    const std::size_t facetsCount = m_pMesh->getFacetsCount();
    constexpr std::size_t noPerFacet = dim;

    const double dt = m_pSolver->getTimeStep();

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
                m_F(facet.getNodeIndex(i)) -= dt*qn[i];
            }
        }

        if(boundaryQh)
        {
            auto SGammah = m_pMatBuilder->getSGamma(facet);

            for(unsigned short i = 0 ; i < noPerFacet ; ++i)
            {
                m_F(facet.getNodeIndex(i)) -= dt*SGammah[i];
            }
        }

        if(boundaryQr)
        {
            auto SGammar = m_pMatBuilder2->getSGamma(facet);

            for(unsigned short i = 0 ; i < noPerFacet ; ++i)
            {
                m_F(facet.getNodeIndex(i)) -= dt*SGammar[i];
            }
        }

    }

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree() && !m_pSolver->getBcTagFlags(node.getTag(), m_bcFlags[0]))
        {
            m_F(n) = node.getState(m_statesIndex[0]);
            invMDiag[n] = 1;
        }
        else if(m_pSolver->getBcTagFlags(node.getTag(), m_bcFlags[0]))
        {
            std::array<double, 1> result;
            result = m_bcParams[0].call<std::array<double, 1>>(m_pMesh->getNodeType(n) + "T",
                                                             node.getPosition(),
                                                             m_pProblem->getCurrentSimTime() +
                                                             m_pSolver->getTimeStep());
            m_F(n) = result[0];
            invMDiag[n] = 1;
        }
    }
}

template<unsigned short dim>
void HeatEqWCompNewton<dim>::m_buildSystem()
{
    m_clock.start();
    constexpr unsigned short nodPerEl = dim + 1;
    const std::size_t elementsCount = m_pMesh->getElementsCount();
    const std::size_t nodesCount = m_pMesh->getNodesCount();

    m_invM.resize(nodesCount); m_invM.setZero();
    std::vector<Eigen::Matrix<double, nodPerEl, nodPerEl>> Me(elementsCount);

    m_F.resize(nodesCount); m_F.setZero();
    std::vector<Eigen::Matrix<double, nodPerEl, 1>> FTote(elementsCount);
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < elementsCount ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::Matrix<double, nodPerEl, 1> T = getElementState<dim>(element, m_statesIndex[0]);

        GradNmatType<dim> gradNe = m_pMatBuilder->getGradN(element);
        BmatType<dim> Be = m_pMatBuilder->getB(gradNe);
        Me[elm] = m_pMatBuilder->getM(element);
        MatrixBuilder<dim>:: template lump<nodPerEl>(Me[elm]);
        Eigen::Matrix<double, nodPerEl, nodPerEl> Le = m_pMatBuilder->getL(element, Be, gradNe);

        FTote[elm] = - m_pSolver->getTimeStep()*Le*T + Me[elm]*T;
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    m_clock.start();
    auto& invMDiag = m_invM.diagonal();

    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            invMDiag[element.getNodeIndex(i)] += Me[elm].diagonal()[i];

            m_F(element.getNodeIndex(i)) += FTote[elm](i);
        }
    }

    MatrixBuilder<dim>::inverse(m_invM);
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();
}

template<unsigned short dim>
double HeatEqWCompNewton<dim>::m_getDflDT(double T)
{
    return static_cast<double>(T > m_Tm - m_DT/2)*static_cast<double>(T < m_Tm + m_DT/2)*(-2*3/(m_DT*m_DT*m_DT)*T*T + 6*2*m_Tm/(m_DT*m_DT*m_DT)*T + (3/(2*m_DT) - 6*m_Tm*m_Tm/(m_DT*m_DT*m_DT)));
}
