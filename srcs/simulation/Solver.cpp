#include "Solver.hpp"

#include <iomanip>

#include "../mesh/Mesh.hpp"
#include "Problem.hpp"
#include "Equation.hpp"

Solver::Solver(Problem* pProblem, Mesh* pMesh, std::vector<SolTable> m_problemParams):
m_timeStep(0),
m_pMesh(pMesh),
m_pProblem(pProblem)
{
    m_solverParams.resize(m_pProblem->getThreadCount());
    for(std::size_t i = 0 ; i < m_solverParams.size() ; ++i)
        m_solverParams[i] = SolTable("Solver", m_problemParams[i]);

    m_id = m_solverParams[0].checkAndGet<std::string>("id");

    m_nextTimeToRemesh = 10;
}

Solver::~Solver()
{

}

void Solver::displayParams() const
{
    throw std::runtime_error("Unimplemented function by the child class -> Solver::displayParams()");
}

void Solver::displayTimeStats() const
{
    for(auto& keyPair : m_accumalatedTimes)
    {
        std::cout << std::defaultfloat << std::setprecision(7) << std::setw(40) << std::left << keyPair.first << ": " << std::setw(10) << std::right << keyPair.second << " s" << std::endl;
    }

    for(auto& pEquation : m_pEquations)
    {
        std::cout << "--------------------------------------" << std::endl;
        pEquation->displayTimeStats();
    }
}

std::string Solver::getID() const noexcept
{
    return m_id;
}

bool Solver::solveOneTimeStep()
{
    throw std::runtime_error("Unimplemented function by the child class -> Solver::solveOneTimeStep()");
}

void Solver::computeNextDT()
{
    throw std::runtime_error("Unimplemented function by the child class -> Solver::computeNextDT()");
}

bool Solver::checkBC(SolTable bcParam, unsigned int n, const Node& node, std::string bcString, unsigned int expectedBCSize)
{
    bool res = bcParam.checkCallNoThrow(m_pMesh->getNodeType(n) + bcString,
                                        node.getPosition(),
                                        m_pMesh->getBoundNodeInitPos(n),
                                        node.getStates(), 0);

    if(res)
    {
        std::vector<double> result = bcParam.call<std::vector<double>>(
            m_pMesh->getNodeType(n) + bcString,
            node.getPosition(),
            m_pMesh->getBoundNodeInitPos(n),
            node.getStates(),
            0);

        if(result.size() != expectedBCSize)
        {
            throw std::runtime_error(" the boundary condition " + m_pMesh->getNodeType(n) + bcString +
                                     " does not return the right number of states: " +
                                     std::to_string(result.size()) + " vs " +  std::to_string(expectedBCSize));
        }
    }

    return res;
}

std::size_t Solver::getAdditionalStateCount() const
{
    throw std::runtime_error("Unimplemented function by the child class -> Solver::getAdditionalStateCount()");
}

void Solver::m_conditionalRemesh(std::size_t speedIndex)
{
    bool force = false;
    if(m_pProblem->getCurrentSimTime() > m_nextTimeToRemesh)
    {
        force = true;
        m_pMesh->remesh(m_pProblem->isOutputVerbose());
        //std::cout << "Remeshing at: " << std::fixed << m_pProblem->getCurrentSimTime() << std::endl;
    }

    double dt = std::numeric_limits<double>::max();
    for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);
        double he = 2*element.getRin();

        double maxSpeed = 0;
        for(std::size_t n = 0 ; n < m_pMesh->getNodesPerElm() ; ++n)
        {
            const Node& node = element.getNode(n);

            double u = 0;
            for(unsigned int d = 0 ; d < m_pMesh->getDim() ; ++d)
            {
                u += node.getState(speedIndex + d)*node.getState(speedIndex + d);
            }
            u = std::sqrt(u);

            maxSpeed = std::max(u, maxSpeed);
        }
        dt = std::min(dt, he/maxSpeed);
    }

    if(force || m_pProblem->getCurrentSimTime() + dt < m_nextTimeToRemesh)
        m_nextTimeToRemesh = m_pProblem->getCurrentSimTime() + dt;
}

