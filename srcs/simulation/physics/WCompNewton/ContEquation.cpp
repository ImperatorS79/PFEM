#include "ContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

ContEqWCompNewton::ContEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "ContEq")
{
    unsigned int nGPHD = 0;
    unsigned int nGPLD = 0;
    if(m_pMesh->getDim() == 2)
    {
        nGPHD = 3;
        nGPLD = 3;
    }
    else
    {
        nGPHD = 4;
        nGPLD = 3;
    }
    m_pMatBuilder2 = std::make_unique<MatrixBuilder>(*pMesh, nGPHD, nGPLD);

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

    Eigen::VectorXd m;
    if(m_pMesh->getDim() == 2)
    {
        m.resize(3);
        m << 1, 1, 0;
    }
    else
    {
        m.resize(6);
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

    m_pMatBuilder->setMcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/) -> double {
        return 1;
    });

    m_pMatBuilder->setDcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N, const Eigen::MatrixXd& /** B **/) -> double {
        return (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
    });

    if(m_pSolver->getID() == "CDS_FIC")
    {
        m_pMatBuilder2->setMcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N) -> double {
            return (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
        });

        m_pMatBuilder2->setLcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N, const Eigen::MatrixXd& /** B **/) -> double {
            return (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
        });

        m_pMatBuilder2->setHcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N, const Eigen::MatrixXd& /** B **/) -> double {
            double rho = (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
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
    m_bodyForce = Eigen::Map<Eigen::VectorXd>(bodyForce.data(), bodyForce.size());

    m_needNormalCurv = false;
}

ContEqWCompNewton::~ContEqWCompNewton()
{

}

void ContEqWCompNewton::displayParams() const
{
    std::cout << "Continuity equation parameters:\n"
              << " * K0: " << m_K0 << " Pa\n"
              << " * K0': " << m_K0p << " \n"
              << " * Density at zero p: " << m_rhoStar << " Kg/m^3\n"
              << " * Strong continuity: " << std::boolalpha << m_strongContinuity << std::endl;
}

double ContEqWCompNewton::getSpeedEquiv(double /** he **/, const Node& node)
{
    return (m_K0 + m_K0p*node.getState(m_pMesh->getDim()))/node.getState(m_pMesh->getDim() + 1);
}

bool ContEqWCompNewton::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Continuity Equation" << std::endl;

    m_clock.start();
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> invM; //The mass matrix of the continuity.
    Eigen::VectorXd qRho(m_pMesh->getNodesCount());
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    if(m_pSolver->getID() == "CDS_FIC")
    {
        m_buildSystemFIC(invM, m_F0);
        m_clock.start();
        m_applyBC(invM, m_F0);
        m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();
        std::size_t indexPrevRho = m_pMesh->getNode(0).getStates().size() - 1;
        m_clock.start();
        qRho = invM*m_F0;
        m_accumalatedTimes["Solve system"] += m_clock.end();

        m_clock.start();
        Eigen::VectorXd qPrevRho = getQFromNodesStates(m_pMesh, m_statesIndex[1], m_statesIndex[1]);
        setNodesStatesfromQ(m_pMesh, qPrevRho, indexPrevRho, indexPrevRho);
        m_accumalatedTimes["Update solutions"] += m_clock.end();

    }
    else
    {
        m_buildSystem(invM, m_F0);
        m_clock.start();
        m_applyBC(invM, m_F0);
        m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();
        m_clock.start();
        qRho = invM*m_F0;
        m_accumalatedTimes["Solve system"] += m_clock.end();
    }

    m_clock.start();
    setNodesStatesfromQ(m_pMesh, qRho, m_statesIndex[1], m_statesIndex[1]);

    Eigen::VectorXd qP = m_getPFromRhoTaitMurnagham(qRho);
    setNodesStatesfromQ(m_pMesh, qP, m_statesIndex[0], m_statesIndex[0]);
    m_accumalatedTimes["Update solutions"] += m_clock.end();

    return true;
}

void ContEqWCompNewton::preCompute()
{
    if(m_strongContinuity)
        m_buildF0(m_F0);
}

void ContEqWCompNewton::m_applyBC(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F0)
{
    auto& invMDiag = invM.diagonal();

    for (std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree())
        {
            F0(n) = m_rhoStar;

            invMDiag[n] = 1;
        }
    }
}

void ContEqWCompNewton::m_buildF0(Eigen::VectorXd& F0)
{
    m_clock.start();
    F0.resize(m_pMesh->getNodesCount()); F0.setZero();

    std::vector<Eigen::VectorXd> F0e(m_pMesh->getElementsCount());
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::VectorXd Rho = getElementState(m_pMesh, element, m_statesIndex[1]);

        Eigen::MatrixXd Mrhoe = m_pMatBuilder->getM(element);

        F0e[elm] = Mrhoe*Rho;
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    m_clock.start();
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(uint8_t i = 0 ; i < m_pMesh->getNodesPerElm() ; ++i)
            F0(element.getNodeIndex(i)) += F0e[elm](i);
    }
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();
}

void ContEqWCompNewton::m_buildSystem(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F0)
{
    m_clock.start();
    const unsigned short dim = m_pMesh->getDim();
    const unsigned short noPerEl = m_pMesh->getNodesPerElm();

    invM.resize(m_pMesh->getNodesCount()); invM.setZero();

    std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> MeLumped(m_pMesh->getElementsCount());
    std::vector<Eigen::VectorXd> F0e;

    if(!m_strongContinuity)
    {
        F0.resize(m_pMesh->getNodesCount()); F0.setZero();
        F0e.resize(m_pMesh->getElementsCount());
    }
    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::MatrixXd Me = m_pMatBuilder->getM(element);
        MeLumped[elm] = MatrixBuilder::lump2(Me);

        if(!m_strongContinuity)
        {
            Eigen::VectorXd Rho = getElementState(m_pMesh, element, m_statesIndex[1]);
            Eigen::VectorXd V(noPerEl*dim);
            if(dim == 2)
            {
                V << getElementState(m_pMesh, element, m_statesIndex[2]),
                     getElementState(m_pMesh, element, m_statesIndex[2] + 1);
            }
            else
            {
                 V << getElementState(m_pMesh, element, m_statesIndex[2]),
                      getElementState(m_pMesh, element, m_statesIndex[2] + 1),
                      getElementState(m_pMesh, element, m_statesIndex[2] + 2);
            }

            Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
            Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);
            Eigen::MatrixXd Drhoe = m_pMatBuilder->getD(element, Be);

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
    auto& invMDiag = invM.diagonal();

    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            invMDiag[element.getNodeIndex(i)] += MeLumped[elm].diagonal()[i];

            if(!m_strongContinuity)
            {
                F0(element.getNodeIndex(i)) += F0e[elm](i);
            }
        }
    }

    MatrixBuilder::inverse(invM);
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();
}

void ContEqWCompNewton::m_buildSystemFIC(Eigen::DiagonalMatrix<double,Eigen::Dynamic>& invM, Eigen::VectorXd& F0)
{
    m_clock.start();
    const unsigned short dim = m_pMesh->getDim();
    const unsigned short noPerEl = m_pMesh->getNodesPerElm();

    invM.resize(m_pMesh->getNodesCount()); invM.setZero();

    std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> MeLumped(m_pMesh->getElementsCount());
    std::vector<Eigen::VectorXd> F0e;
    F0.resize(m_pMesh->getNodesCount()); F0.setZero();
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

        Eigen::VectorXd Rho = getElementState(m_pMesh, element, m_statesIndex[1]);
        Eigen::VectorXd prevRho = getElementState(m_pMesh, element, indexPrevRho);
        Eigen::VectorXd P = getElementState(m_pMesh, element, m_statesIndex[0]);
        Eigen::VectorXd V(noPerEl*dim);
        Eigen::VectorXd A(noPerEl*dim);
        if(dim == 2)
        {
            V << getElementState(m_pMesh, element, m_statesIndex[2]),
                 getElementState(m_pMesh, element, m_statesIndex[2] + 1);

            A << getElementState(m_pMesh, element, dim + 2),
                 getElementState(m_pMesh, element, dim + 3);
        }
        else
        {
             V << getElementState(m_pMesh, element, m_statesIndex[2]),
                  getElementState(m_pMesh, element, m_statesIndex[2] + 1),
                  getElementState(m_pMesh, element, m_statesIndex[2] + 2);

            A << getElementState(m_pMesh, element, dim + 2),
                  getElementState(m_pMesh, element, dim + 3),
                  getElementState(m_pMesh, element, dim + 4);
        }
        Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
        Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);

        double tau = m_computeTauFIC(element);
        Eigen::MatrixXd Me = m_pMatBuilder->getM(element);
        MatrixBuilder::lump(Me);
        Eigen::MatrixXd Drhoe = m_pMatBuilder->getD(element, Be);

        Eigen::MatrixXd Me_tau_rho = tau*m_pMatBuilder2->getM(element);
        MatrixBuilder::lump(Me_tau_rho);
        Eigen::MatrixXd Le_tau_rho = tau*m_pMatBuilder2->getL(element, Be, gradNe);
        Eigen::MatrixXd He_tau_rho_2 = tau*m_pMatBuilder2->getH(element, m_bodyForce, Be, gradNe);

        F0e[elm] = Me*Rho + dt*(-Drhoe*V - Le_tau_rho*P + He_tau_rho_2) + (1/dt)*Me_tau_rho*(2*Rho - prevRho);
        MeLumped[elm] = MatrixBuilder::lump2(Me + (1/dt)*Me_tau_rho);

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
    auto& invMDiag = invM.diagonal();

    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            invMDiag[element.getNodeIndex(i)] += MeLumped[elm].diagonal()[i];
            F0(element.getNodeIndex(i)) += F0e[elm](i);
        }
    }

    MatrixBuilder::inverse(invM);
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();
}

Eigen::VectorXd ContEqWCompNewton::m_getPFromRhoTaitMurnagham(const Eigen::VectorXd& qRho)
{
    Eigen::VectorXd qP(m_pMesh->getNodesCount());

    if(m_pProblem->getID() == "BoussinesqWC")
    {
        #pragma omp parallel for default(shared)
        for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
        {
            double rho = qRho[n];
            double T = m_pMesh->getNode(n).getState(m_statesIndex[3]);
            qP[n] = (m_K0/m_K0p)*(std::pow(rho/(m_rhoStar*(1 - m_alpha*(T - m_Tr))), m_K0p) - 1);
        }
    }
    else
    {
        #pragma omp parallel for default(shared)
        for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
        {
            double rho = qRho[n];
            qP[n] = (m_K0/m_K0p)*(std::pow(rho/m_rhoStar, m_K0p) - 1);
        }
    }

    return qP;
}

double ContEqWCompNewton::m_computeTauFIC(const Element& element) const
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
