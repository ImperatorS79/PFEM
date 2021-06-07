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

    m_k = m_materialParams[0].checkAndGet<double>("k");
    m_cv = m_materialParams[0].checkAndGet<double>("cv");

    if(bcFlags.size() != 2)
        throw std::runtime_error("the " + getID() + " equation require two flags for two possible boundary conditions!");

    if(statesIndex.size() != 2)
        throw std::runtime_error("the " + getID() + " equation require two states index describing the T and rho state !");


    m_pMatBuilder->setMcomputeFactor([&](const Element& element,
                                         const NmatTypeHD<dim>& N) -> double {
        double rho = (N*getElementState<dim>(element, m_statesIndex[1])).value();
        return m_cv*rho;
    });

    m_pMatBuilder->setLcomputeFactor([&](const Element& /** element **/,
                                         const NmatTypeHD<dim>& /** N **/,
                                         const BmatType<dim>& /** B **/) -> double {
        return m_k;
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
    //const std::size_t facetsCount = m_pMesh->getFacetsCount();
    //constexpr unsigned short noPerFacet = dim;

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree())
        {
            m_F(n) = node.getState(m_statesIndex[0]);
            invMDiag[n] = 1;
        }
        else if(node.getFlag(m_bcFlags[0]))
        {
            std::array<double, 1> result;
            result = m_bcParams[0].call<std::array<double, 1>>(m_pMesh->getNodeType(n) + "T",
                                                             node.getPosition(),
                                                             m_pMesh->getBoundNodeInitPos(n),
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
