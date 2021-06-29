#include "ContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

template<unsigned short dim>
ContEqWCompNewton<dim>::ContEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "ContEq")
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

    m_K0 = m_materialParams[0].checkAndGet<double>("K0");
    m_K0p = m_materialParams[0].checkAndGet<double>("K0p");
    m_rhoStar = m_materialParams[0].checkAndGet<double>("rhoStar");
    m_mu = m_materialParams[0].checkAndGet<double>("mu");
    std::string stab = m_equationParams[0].checkAndGet<std::string>("stabilization");

    if(stab == "None")
        m_stabilization = Stab::None;
    else if (stab == "Meduri")
        m_stabilization = Stab::Meduri;
    else
        throw std::runtime_error("unknown stabilization: " + stab);

    if(m_pSolver->getID() == "CDS_rho")
        m_version = EqType::Rho;
    else if(m_pSolver->getID() == "CDS_drhodt")
        m_version = EqType::DRhoDt;
    else if(m_pSolver->getID() == "CDS_dpdt")
        m_version = EqType::DPDt;

    if(m_pProblem->getID() == "BoussinesqWC")
    {
        m_alpha = m_materialParams[0].checkAndGet<double>("alpha");
        m_Tr = m_materialParams[0].checkAndGet<double>("Tr");
    }

    mVecType<dim> m;
    if constexpr (dim == 2)
    {
        m << 1, 1, 0;
    }
    else if constexpr (dim == 3)
    {
        m << 1, 1, 1, 0, 0, 0;
    }
    m_pMatBuilder->setm(m);

    if(m_pProblem->getID() == "BoussinesqWC")
    {
        if(m_statesIndex.size()!= 4)
            throw std::runtime_error("the " + getID() + " equation requires 4 statesIndex: one index for the p unknown, one for the rho unknown, one for the beginning of the (u,v,w) unknown, and the T unknown!");

    }
    else
    {
        if(m_statesIndex.size()!= 3)
            throw std::runtime_error("the " + getID() + " equation requires 3 statesIndex: one index for the p unknown, one for the rho unknown, and one for the beginning of the (u,v,w) unknown!");
    }

    m_pMatBuilder->setMcomputeFactor([&](const Element& /** element **/,
                                         const NmatTypeHD<dim>& /** N **/) -> double {
        return 1;
    });

    if(m_version == EqType::DPDt)
    {
        m_pMatBuilder->setDcomputeFactor([&](const Element& element,
                                         const NmatTypeHD<dim>& N,
                                         const BmatType<dim>& /** B **/) -> double {
            return m_K0 + m_K0p*(N*getElementState<dim>(element, m_statesIndex[0])).value();
        });
    }
    else
    {
        m_pMatBuilder->setDcomputeFactor([&](const Element& element,
                                         const NmatTypeHD<dim>& N,
                                         const BmatType<dim>& /** B **/) -> double {
            return (N*getElementState<dim>(element, m_statesIndex[1])).value();
        });
    }

    m_needNormalCurv = false;
}

template<unsigned short dim>
ContEqWCompNewton<dim>::~ContEqWCompNewton()
{

}

template<unsigned short dim>
void ContEqWCompNewton<dim>::displayParams() const
{
    std::cout << "Continuity equation parameters:\n"
              << " * K0: " << m_K0 << " Pa\n"
              << " * K0': " << m_K0p << " \n"
              << " * Density at zero p: " << m_rhoStar << " Kg/m^3\n"
              << " * Version: " << static_cast<int>(m_version) << std::endl;
}

template<unsigned short dim>
double ContEqWCompNewton<dim>::getSquaredSpeedEquiv(const Node& node) const
{
    return (m_K0 + m_K0p*node.getState(dim))/node.getState(dim + 1);
}

template<unsigned short dim>
bool ContEqWCompNewton<dim>::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Continuity Equation" << std::endl;

    m_clock.start();
    Eigen::VectorXd qRho(m_pMesh->getNodesCount());
    Eigen::VectorXd qP(m_pMesh->getNodesCount());
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    if(m_version == EqType::DPDt)
    {
        m_buildSystemdpdt();
        m_clock.start();
        m_applyBCdpdt();
        m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();
        qP = m_invM*m_F0;
        m_accumalatedTimes["Solve system"] += m_clock.end();

        m_clock.start();
        setNodesStatesfromQ(m_pMesh, qP, m_statesIndex[0], m_statesIndex[0]);

        qRho = m_getRhoFromPTaitMurnagham(qP);
        setNodesStatesfromQ(m_pMesh, qRho, m_statesIndex[1], m_statesIndex[1]);
        m_accumalatedTimes["Update solutions"] += m_clock.end();
    }
    else
    {
        m_buildSystem();
        m_clock.start();
        m_applyBC();
        m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();
        m_clock.start();
        qRho = m_invM*m_F0;
        m_accumalatedTimes["Solve system"] += m_clock.end();

        m_clock.start();
        setNodesStatesfromQ(m_pMesh, qRho, m_statesIndex[1], m_statesIndex[1]);

        qP = m_getPFromRhoTaitMurnagham(qRho);
        setNodesStatesfromQ(m_pMesh, qP, m_statesIndex[0], m_statesIndex[0]);
        m_accumalatedTimes["Update solutions"] += m_clock.end();
    }

    return true;
}

template<unsigned short dim>
void ContEqWCompNewton<dim>::preCompute()
{
    if(m_version == EqType::Rho)
        m_buildF0();
}

template<unsigned short dim>
void ContEqWCompNewton<dim>::m_applyBC()
{
    auto& invMDiag = m_invM.diagonal();

    for (std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree() || node.isOnFreeSurface())
        {
            m_F0(n) = m_rhoStar;

            invMDiag[n] = 1;
        }
    }
}

template<unsigned short dim>
void ContEqWCompNewton<dim>::m_buildF0()
{
    m_clock.start();
    constexpr unsigned short nodPerEl = dim + 1;

    m_F0.resize(m_pMesh->getNodesCount()); m_F0.setZero();

    std::vector<Eigen::Matrix<double, nodPerEl, 1>> F0e(m_pMesh->getElementsCount());
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::Matrix<double, nodPerEl, 1> Rho = getElementState<dim>(element, m_statesIndex[1]);

        Eigen::Matrix<double, nodPerEl, nodPerEl> Mrhoe = m_pMatBuilder->getM(element);

        F0e[elm] = Mrhoe*Rho;
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    m_clock.start();
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(uint8_t i = 0 ; i < m_pMesh->getNodesPerElm() ; ++i)
            m_F0(element.getNodeIndex(i)) += F0e[elm](i);
    }
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();
}

template<unsigned short dim>
void ContEqWCompNewton<dim>::m_buildSystem()
{
    m_clock.start();
    constexpr unsigned short nodPerEl = dim + 1;

    m_invM.resize(m_pMesh->getNodesCount()); m_invM.setZero();

    std::vector<Eigen::DiagonalMatrix<double, nodPerEl>> MeLumped(m_pMesh->getElementsCount());
    std::vector<Eigen::Matrix<double, nodPerEl, 1>> F0e;

    if(m_version == EqType::DRhoDt)
    {
        m_F0.resize(m_pMesh->getNodesCount()); m_F0.setZero();
        F0e.resize(m_pMesh->getElementsCount());
    }
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::Matrix<double, nodPerEl, nodPerEl> Me = m_pMatBuilder->getM(element);
        MeLumped[elm] = MatrixBuilder<dim>:: template lump2<nodPerEl>(Me);

        if(m_version == EqType::DRhoDt)
        {
            Eigen::Matrix<double, nodPerEl, 1> Rho = getElementState<dim>(element, m_statesIndex[1]);
            Eigen::Matrix<double, dim*nodPerEl, 1> V = getElementVecState<dim>(element, m_statesIndex[2]);

            GradNmatType<dim> gradNe = m_pMatBuilder->getGradN(element);
            BmatType<dim> Be = m_pMatBuilder->getB(gradNe);
            Eigen::Matrix<double, nodPerEl, dim*nodPerEl> Drhoe = m_pMatBuilder->getD(element, Be);

            F0e[elm] = - m_pSolver->getTimeStep()*Drhoe*V;

            if(m_stabilization != Stab::None)
                F0e[elm] += Me*Rho;
            else
                F0e[elm] += MeLumped[elm]*Rho;
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    m_clock.start();
    auto& invMDiag = m_invM.diagonal();

    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(unsigned short i = 0 ; i < nodPerEl ; ++i)
        {
            invMDiag[element.getNodeIndex(i)] += MeLumped[elm].diagonal()[i];

            if(m_version == EqType::DRhoDt)
            {
                m_F0(element.getNodeIndex(i)) += F0e[elm](i);
            }
        }
    }

    MatrixBuilder<dim>::inverse(m_invM);
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();
}

template<unsigned short dim>
Eigen::VectorXd ContEqWCompNewton<dim>::m_getPFromRhoTaitMurnagham(const Eigen::VectorXd& qRho)
{
    Eigen::VectorXd qP(m_pMesh->getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
    {
        double rho = qRho[n];
        qP[n] = (m_K0/m_K0p)*(std::pow(rho/m_rhoStar, m_K0p) - 1);
    }

    return qP;
}

template<unsigned short dim>
Eigen::VectorXd ContEqWCompNewton<dim>::m_getRhoFromPTaitMurnagham(const Eigen::VectorXd& qP)
{
    Eigen::VectorXd qRho(m_pMesh->getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
    {
        double p = qP[n];
        qRho[n] = std::pow((m_K0p/m_K0)*p + 1, 1/m_K0p)*m_rhoStar;
    }

    return qRho;
}

template<unsigned short dim>
void ContEqWCompNewton<dim>::m_applyBCdpdt()
{
    auto& invMDiag = m_invM.diagonal();

    #pragma omp parallel for default(shared) schedule(static)
    for (std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree())
        {
            m_F0(n) = 0;

            invMDiag[n] = 1;
        }
    }
}

template<unsigned short dim>
void ContEqWCompNewton<dim>::m_buildSystemdpdt()
{
    m_clock.start();
    constexpr unsigned short nodPerEl = dim + 1;

    m_invM.resize(m_pMesh->getNodesCount()); m_invM.setZero();
    m_MeLumped.resize(m_pMesh->getElementsCount());

    m_F0.resize(m_pMesh->getNodesCount()); m_F0.setZero();
    m_F0e.resize(m_pMesh->getElementsCount());

    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    double dt = m_pSolver->getTimeStep();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::Matrix<double, nodPerEl, nodPerEl> Me = m_pMatBuilder->getM(element);
        m_MeLumped[elm] = MatrixBuilder<dim>:: template lump2<nodPerEl>(Me);

        Eigen::Matrix<double, nodPerEl, 1> P = getElementState<dim>(element, m_statesIndex[0]);
        Eigen::Matrix<double, dim*nodPerEl, 1> V = getElementVecState<dim>(element, m_statesIndex[2]);

        GradNmatType<dim> gradNe = m_pMatBuilder->getGradN(element);
        BmatType<dim> Be = m_pMatBuilder->getB(gradNe);
        Eigen::Matrix<double, nodPerEl, dim*nodPerEl> Drhoe = m_pMatBuilder->getD(element, Be);

        m_F0e[elm] = - dt*Drhoe*V;

        if(m_stabilization == Stab::Meduri)
            m_F0e[elm] += Me*P;
        else
            m_F0e[elm] += m_MeLumped[elm]*P;
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    m_clock.start();
    auto& invMDiag = m_invM.diagonal();

    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(unsigned short i = 0 ; i < nodPerEl ; ++i)
        {
            invMDiag[element.getNodeIndex(i)] += m_MeLumped[elm].diagonal()[i];

            m_F0(element.getNodeIndex(i)) += m_F0e[elm](i);
        }
    }

    MatrixBuilder<dim>::inverse(m_invM);

    m_accumalatedTimes["Assemble matrix"] += m_clock.end();
}
