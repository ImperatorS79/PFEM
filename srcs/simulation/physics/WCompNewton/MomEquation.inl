#include "MomEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"
#include "../../utility/Clock.hpp"

template<unsigned short dim>
MomEqWCompNewton<dim>::MomEqWCompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "MomEq")
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

    m_mu = m_materialParams[0].checkAndGet<double>("mu");
    m_gamma = m_materialParams[0].checkAndGet<double>("gamma");

    if(bcFlags.size() != 1)
        throw std::runtime_error("the " + getID() + " only required one flag for one possible boundary condition!");

    if(m_pProblem->getID() == "BoussinesqWC")
    {
        m_alpha = m_materialParams[0].checkAndGet<double>("alpha");
        m_Tr = m_materialParams[0].checkAndGet<double>("Tr");

        if(statesIndex.size() != 5)
            throw std::runtime_error("the " + getID() + " equation requires 5 statesIndex: beginning of (u,v,w), beginning of (ax, ay, az), p and rho and T");
    }
    else
    {
        if(statesIndex.size() != 4)
            throw std::runtime_error("the " + getID() + "  equation requires 4 statesIndex: beginning of (u,v,w), beginning of (ax, ay, az), p and rho");
    }

    DdevMatType<dim> ddev;
    mVecType<dim> m;
    if constexpr (dim == 2)
    {
        ddev <<  4.0/3, -2.0/3, 0,
                -2.0/3,  4.0/3, 0,
                     0,      0, 1;

        m << 1, 1, 0;
    }
    else if constexpr (dim == 3)
    {
        ddev <<  4.0/3, -2.0/3, -2.0/3, 0, 0, 0,
                -2.0/3,  4.0/3, -2.0/3, 0, 0, 0,
                -2.0/3, -2.0/3,  4.0/3, 0, 0, 0,
                     0,      0,      0, 1, 0, 0,
                     0,      0,      0, 0, 1, 0,
                     0,      0,      0, 0, 0, 1;

        m << 1, 1, 1, 0, 0, 0;
    }
    m_pMatBuilder->setddev(ddev);
    m_pMatBuilder->setm(m);

    m_pMatBuilder->setMcomputeFactor([&](const Element& element,
                                         const NmatTypeHD<dim>& N) -> double {
        return (N*getElementState<dim>(element, m_statesIndex[3])).value();
    });

    m_pMatBuilder->setKcomputeFactor([&](const Element& /** element **/,
                                         const NmatTypeHD<dim>& /** N **/,
                                         const BmatType<dim>& /** B **/,
                                         const DdevMatType<dim>& /** ddev **/) -> double {
        return m_mu;
    });

    m_pMatBuilder->setDcomputeFactor([&](const Element& /** element **/,
                                         const NmatTypeHD<dim>& /** N **/,
                                         const BmatType<dim>& /** B **/) -> double {
        return 1;
    });

    if(m_pProblem->getID() == "BoussinesqWC")
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& element,
                                             const NmatTypeHD<dim>&  N,
                                             const BmatType<dim>& /** B **/) -> double {
            double rho = (N*getElementState<dim>(element, m_statesIndex[3])).value();
            double T = (N*getElementState<dim>(element, m_statesIndex[4])).value();
            return rho*(1 - m_alpha*(T - m_Tr));
        });
    }
    else
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& element,
                                             const NmatTypeHD<dim>&  N,
                                             const BmatType<dim>& /** B **/) -> double {
            return (N*getElementState<dim>(element, m_statesIndex[3])).value();
        });
    }

    auto bodyForce = m_equationParams[0].checkAndGet<std::vector<double>>("bodyForce");
    if(bodyForce.size() != m_pMesh->getDim())
        throw std::runtime_error("the body force vector has not the right dimension!");

    m_bodyForce = Eigen::Map<Eigen::Matrix<double, dim, 1>>(bodyForce.data(), bodyForce.size());

    m_needNormalCurv = (m_gamma < 1e-15) ? false : true;
}

template<unsigned short dim>
MomEqWCompNewton<dim>::~MomEqWCompNewton()
{

}

template<unsigned short dim>
void MomEqWCompNewton<dim>::displayParams() const
{
    std::cout << "Momentum equation parameters:\n"
              << " * Viscosity: " << m_mu << " Pa s\n"
              << " * Surface Tension: " << m_gamma << " Nm" << std::endl;

    if(m_pMesh->getDim() == 2)
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << ")" << std::endl;
    else
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << "," << m_bodyForce[2] << ")" << std::endl;
}

template<unsigned short dim>
double MomEqWCompNewton<dim>::getSquaredSpeedEquiv(const Node& node) const
{
    double nodeU2 = 0;
    for(unsigned int d = 0 ; d < dim ; ++d)
    {
        nodeU2 += node.getState(d)*node.getState(d);
    }
    return nodeU2;
}

template<unsigned short dim>
double MomEqWCompNewton<dim>::getDiffusionParam(const Node& node) const
{
    return m_mu/node.getState(m_statesIndex[3]);
}

template<unsigned short dim>
bool MomEqWCompNewton<dim>::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Momentum Equation" << std::endl;

    m_clock.start();
    Eigen::VectorXd qV1half = getQFromNodesStates(m_pMesh, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim() - 1);
    m_accumalatedTimes["Update solutions"] += m_clock.end();

    m_buildSystem();
    m_clock.start();
    m_applyBC();
    m_accumalatedTimes["Apply boundary conditions"] += m_clock.end();
    m_clock.start();
    Eigen::VectorXd qAcc = m_invM*m_F;
    m_accumalatedTimes["Solve system"] += m_clock.end();

    m_clock.start();
    Eigen::VectorXd qV = qV1half + 0.5*m_pSolver->getTimeStep()*qAcc;
    setNodesStatesfromQ(m_pMesh, qV, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim() - 1);
    setNodesStatesfromQ(m_pMesh, qAcc, m_statesIndex[1], m_statesIndex[1] + m_pMesh->getDim() - 1);
    m_accumalatedTimes["Update solutions"] += m_clock.end();

    return true;
}

template<unsigned short dim>
void MomEqWCompNewton<dim>::m_buildSystem()
{
    m_clock.start();
    const std::size_t elementsCount = m_pMesh->getElementsCount();
    const std::size_t nodesCount = m_pMesh->getNodesCount();
    constexpr unsigned short nodPerEl = dim + 1;

    m_invM.resize(dim*nodesCount); m_invM.setZero();
    std::vector<Eigen::Matrix<double, dim*nodPerEl, dim*nodPerEl>> Me(elementsCount);

    m_F.resize(dim*nodesCount); m_F.setZero();
    std::vector<Eigen::Matrix<double, dim*nodPerEl, 1>> FTote(elementsCount);

    m_accumalatedTimes["Prepare matrix assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < elementsCount ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::Matrix<double, dim*nodPerEl, 1> V = getElementVecState<dim>(element, m_statesIndex[0]);
        Eigen::Matrix<double, nodPerEl, 1> P = getElementState<dim>(element, m_statesIndex[2]);

        GradNmatType<dim> gradNe = m_pMatBuilder->getGradN(element);
        BmatType<dim> Be = m_pMatBuilder->getB(gradNe);

        Eigen::Matrix<double, nodPerEl, nodPerEl> MeTemp = m_pMatBuilder->getM(element);
        Me[elm] = MatrixBuilder<dim>::diagBlock(MeTemp);
        MatrixBuilder<dim>:: template lump<dim*nodPerEl>(Me[elm]);
        Eigen::Matrix<double, dim*nodPerEl, dim*nodPerEl> Ke = m_pMatBuilder->getK(element, Be);
        Eigen::Matrix<double, nodPerEl, dim*nodPerEl> De = m_pMatBuilder->getD(element, Be);
        Eigen::Matrix<double, dim*nodPerEl, 1> Fe = m_pMatBuilder->getF(element, m_bodyForce, Be);

        FTote[elm] = -Ke*V + De.transpose()*P + Fe;
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    m_clock.start();
    auto& invMDiag = m_invM.diagonal();

    for(std::size_t elm = 0 ; elm < elementsCount ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                /********************************************************************
                                             Build M
                ********************************************************************/
                invMDiag[element.getNodeIndex(i) + d*nodesCount] += Me[elm](i + d*nodPerEl, i + d*nodPerEl);

                /************************************************************************
                                                Build f
                ************************************************************************/
                m_F(element.getNodeIndex(i) + d*nodesCount) += FTote[elm](i + d*nodPerEl);
            }
        }
    }

    MatrixBuilder<dim>::inverse(m_invM);
    m_accumalatedTimes["Assemble matrix"] += m_clock.end();
}

template<unsigned short dim>
void MomEqWCompNewton<dim>::m_applyBC()
{
    assert(m_pMesh->getNodesCount() != 0);

    const std::size_t nodesCount = m_pMesh->getNodesCount();

    auto& invMDiag = m_invM.diagonal();

    //Do not parallelize this
    #pragma omp parallel for default(shared) schedule(dynamic)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        int threadIndex = omp_get_thread_num();
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree() && !node.isBound())
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                m_F(n + d*nodesCount) = m_bodyForce[d];
                invMDiag[n + d*nodesCount] = 1;
            }
        }
        else if(node.isBound())
        {
            if(node.getFlag(m_bcFlags[0]))
            {
                std::array<double, dim> result;
                result = m_bcParams[threadIndex].call<std::array<double, dim>>(m_pMesh->getNodeType(n) + "V",
                                                        node.getPosition(),
                                                        m_pMesh->getBoundNodeInitPos(n),
                                                        m_pProblem->getCurrentSimTime() +
                                                        m_pSolver->getTimeStep());

                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    m_F(n + d*nodesCount) = result[d];
                    invMDiag[n + d*nodesCount] = 1;
                }
            }
        }
    }
}
