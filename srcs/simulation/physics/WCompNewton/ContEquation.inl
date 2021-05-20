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
    m_pMatBuilder2 = std::make_unique<MatrixBuilder<dim>>(*pMesh, nGPHD, nGPLD);

    m_K0 = m_materialParams[0].checkAndGet<double>("K0");
    m_K0p = m_materialParams[0].checkAndGet<double>("K0p");
    m_rhoStar = m_materialParams[0].checkAndGet<double>("rhoStar");
    m_mu = m_materialParams[0].checkAndGet<double>("mu");
    m_strongContinuity = m_equationParams[0].checkAndGet<bool>("strongContinuity");
    m_enableStab = m_equationParams[0].checkAndGet<bool>("enableStab");
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

    m_pMatBuilder->setDcomputeFactor([&](const Element& element,
                                         const NmatTypeHD<dim>& N,
                                         const BmatType<dim>& /** B **/) -> double {
        return (N*getElementState<dim>(element, m_statesIndex[1])).value();
    });

    if(m_pSolver->getID() == "CDS_FIC")
    {
        m_pMatBuilder2->setMcomputeFactor([&](const Element& element,
                                              const NmatTypeHD<dim>& N) -> double {
            return (N*getElementState<dim>(element, m_statesIndex[1])).value();
        });

        m_pMatBuilder2->setLcomputeFactor([&](const Element& element,
                                              const NmatTypeHD<dim>& N,
                                              const BmatType<dim>& /** B **/) -> double {
            return (N*getElementState<dim>(element, m_statesIndex[1])).value();
        });

        m_pMatBuilder2->setHcomputeFactor([&](const Element& element,
                                              const NmatTypeHD<dim>& N,
                                              const BmatType<dim>& /** B **/) -> double {
            double rho = (N*getElementState<dim>(element, m_statesIndex[1])).value();
            return rho*rho;
        });

//        m_pMatBuilder2->setCcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N, const Eigen::MatrixXd& /** B **/) -> double {
//            return (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
//        });
//
//        m_pMatBuilder2->setLcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
//            return 1;
//        });
//
//        m_pMatBuilder2->setHcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N, const Eigen::MatrixXd& /** B **/) -> double {
//            return (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
//        });
    }

    auto bodyForce = m_equationParams[0].checkAndGet<std::vector<double>>("bodyForce");
    if(bodyForce.size() != m_pMesh->getDim())
        throw std::runtime_error("the body force vector has not the right dimension!");
    m_bodyForce = Eigen::Map<Eigen::Matrix<double, dim, 1>>(bodyForce.data(), bodyForce.size());

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
              << " * Strong continuity: " << std::boolalpha << m_strongContinuity << std::endl;
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
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    if(m_pSolver->getID() == "CDS_FIC")
    {
        m_buildSystemFIC();
        m_clock.start();
        m_applyBC();
        m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();
        std::size_t indexPrevRho = m_pMesh->getNode(0).getStates().size() - 1;
        m_clock.start();
        qRho = m_invM*m_F0;
        m_accumalatedTimes["Solve system"] += m_clock.end();

        m_clock.start();
        Eigen::VectorXd qPrevRho = getQFromNodesStates(m_pMesh, m_statesIndex[1], m_statesIndex[1]);
        setNodesStatesfromQ(m_pMesh, qPrevRho, indexPrevRho, indexPrevRho);
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
    }

    m_clock.start();
    setNodesStatesfromQ(m_pMesh, qRho, m_statesIndex[1], m_statesIndex[1]);

    Eigen::VectorXd qP = m_getPFromRhoTaitMurnagham(qRho);
    setNodesStatesfromQ(m_pMesh, qP, m_statesIndex[0], m_statesIndex[0]);
    m_accumalatedTimes["Update solutions"] += m_clock.end();

    return true;
}

template<unsigned short dim>
void ContEqWCompNewton<dim>::preCompute()
{
    if(m_strongContinuity)
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

    if(!m_strongContinuity)
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

        if(!m_strongContinuity)
        {
            Eigen::Matrix<double, nodPerEl, 1> Rho = getElementState<dim>(element, m_statesIndex[1]);
            Eigen::Matrix<double, dim*nodPerEl, 1> V = getElementVecState<dim>(element, m_statesIndex[2]);

            GradNmatType<dim> gradNe = m_pMatBuilder->getGradN(element);
            BmatType<dim> Be = m_pMatBuilder->getB(gradNe);
            Eigen::Matrix<double, nodPerEl, dim*nodPerEl> Drhoe = m_pMatBuilder->getD(element, Be);

            F0e[elm] = - m_pSolver->getTimeStep()*Drhoe*V;

            if(m_enableStab)
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

            if(!m_strongContinuity)
            {
                m_F0(element.getNodeIndex(i)) += F0e[elm](i);
            }
        }
    }

    MatrixBuilder<dim>::inverse(m_invM);
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();
}

template<unsigned short dim>
void ContEqWCompNewton<dim>::m_buildSystemFIC()
{
    m_clock.start();
    constexpr unsigned short nodPerEl = dim + 1;

    m_invM.resize(m_pMesh->getNodesCount()); m_invM.setZero();

    std::vector<Eigen::DiagonalMatrix<double, nodPerEl>> MeLumped(m_pMesh->getElementsCount());
    std::vector<Eigen::Matrix<double, nodPerEl, 1>> F0e;
    m_F0.resize(m_pMesh->getNodesCount()); m_F0.setZero();
    F0e.resize(m_pMesh->getElementsCount());

    std::size_t indexPrevRho = m_pMesh->getNode(0).getStates().size() - 1;
    double dt = m_pSolver->getTimeStep();
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::Matrix<double, nodPerEl, 1> Rho = getElementState<dim>(element, m_statesIndex[1]);
        Eigen::Matrix<double, nodPerEl, 1> prevRho = getElementState<dim>(element, indexPrevRho);
        Eigen::Matrix<double, nodPerEl, 1> P = getElementState<dim>(element, m_statesIndex[0]);
        Eigen::Matrix<double, dim*nodPerEl, 1> V = getElementVecState<dim>(element, m_statesIndex[2]);
        //Eigen::Matrix<double, dim*nodPerEl, 1> A = getElementVecState<dim>(element, dim + 2);

        GradNmatType<dim> gradNe = m_pMatBuilder->getGradN(element);
        BmatType<dim> Be = m_pMatBuilder->getB(gradNe);

        double tau = m_computeTauFIC(element);
        Eigen::Matrix<double, nodPerEl, nodPerEl> Me = m_pMatBuilder->getM(element);
        MatrixBuilder<dim>:: template lump<nodPerEl>(Me);
        Eigen::Matrix<double, nodPerEl, dim*nodPerEl> Drhoe = m_pMatBuilder->getD(element, Be);

        Eigen::Matrix<double, nodPerEl, nodPerEl> Me_tau_rho = tau*m_pMatBuilder2->getM(element);
        MatrixBuilder<dim>:: template lump<nodPerEl>(Me_tau_rho);
        Eigen::Matrix<double, nodPerEl, nodPerEl> Le_tau_rho = tau*m_pMatBuilder2->getL(element, Be, gradNe);
        Eigen::Matrix<double, nodPerEl, 1> He_tau_rho_2 = tau*m_pMatBuilder2->getH(element, m_bodyForce, Be, gradNe);

        F0e[elm] = Me*Rho + dt*(-Drhoe*V - Le_tau_rho*P + He_tau_rho_2) + (1/dt)*Me_tau_rho*(2*Rho - prevRho);
        MeLumped[elm] = MatrixBuilder<dim>:: template lump2<nodPerEl>(Me + (1/dt)*Me_tau_rho);

//        Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
//        Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);
//
//        double tau = m_computeTauFIC(element);
//        Eigen::MatrixXd Me = m_pMatBuilder->getM(element);
//        MatrixBuilder::lump(Me);
//        Eigen::MatrixXd Drhoe = m_pMatBuilder->getD(element, Be);
//
//        Eigen::MatrixXd Ce = tau*m_pMatBuilder2->getC(element, Be, gradNe);
//        Eigen::MatrixXd Le = tau*m_pMatBuilder2->getL(element, Be, gradNe);
//        Eigen::MatrixXd He = tau*m_pMatBuilder2->getH(element, m_bodyForce, Be, gradNe);
//
//        F0e[elm] = Me*Rho + dt*(-Drhoe*V + Ce*A + Le*P - He);
//        MeLumped[elm] = MatrixBuilder::lump2(Me);
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
            m_F0(element.getNodeIndex(i)) += F0e[elm](i);
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
double ContEqWCompNewton<dim>::m_computeTauFIC(const Element& element) const
{
    const double h = std::sqrt(m_pMesh->getRefElementSize(m_pMesh->getDim())*element.getDetJ()/M_PI);

    double rho = 0;
    for (unsigned short n = 0 ; n < m_pMesh->getNodesPerElm() ; ++n)
    {
        const Node& node = m_pMesh->getNode(element.getNodeIndex(n));
        rho += node.getState(m_statesIndex[1]);
    }
    rho /= (m_pMesh->getDim() + 1);

    return 1/(8*m_mu/(h*h) + 2*rho/m_pSolver->getTimeStep());

//    const double h = std::sqrt(m_pMesh->getRefElementSize(m_pMesh->getDim())*element.getDetJ()/M_PI);
//
//    double U = 0;
//    double Rho = 0;
//    for (unsigned short n = 0 ; n < m_pMesh->getNodesPerElm() ; ++n)
//    {
//        const Node& node = m_pMesh->getNode(element.getNodeIndex(n));
//
//        double nodeU = 0;
//        double nodeRho = node.getState(m_statesIndex[1]);
//        for (unsigned short d = 0 ; d < m_pMesh->getDim() ; ++d)
//        {
//            nodeU += node.getState(d)*node.getState(d);
//        }
//        U += std::sqrt(nodeU);
//        Rho += nodeRho;
//    }
//    U /= (m_pMesh->getDim() + 1);
//    Rho /= (m_pMesh->getDim() + 1);
//
//    return 1/std::sqrt((2/m_pSolver->getTimeStep())*(2/m_pSolver->getTimeStep()) + (2*U/h)*(2*U/h)
//                        + 9*(4*m_mu/(h*h*Rho))*(4*m_mu/(h*h*Rho)));
}
