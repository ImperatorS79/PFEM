#include "Equation.hpp"

#include <iomanip>

#include "Problem.hpp"

Equation::Equation(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                 std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                 const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex,
                 const std::string& id):
m_id(id),
m_materialParams(materialParams),
m_bcFlags(bcFlags),
m_statesIndex(statesIndex),
m_pProblem(pProblem),
m_pSolver(pSolver),
m_pMesh(pMesh)
{
    m_equationParams.resize(m_pProblem->getThreadCount());
    m_bcParams.resize(m_pProblem->getThreadCount());
    for(std::size_t i = 0 ; i < solverParams.size() ; ++i)
    {
        m_equationParams[i] = SolTable(getID(), solverParams[i]);
        m_bcParams[i] = SolTable("BC", m_equationParams[i]);
    }
}

Equation::~Equation()
{

}

void Equation::displayParams() const
{
    throw std::runtime_error("Unimplemented function by the child class -> Equation::displayParams()");
}

void Equation::displayTimeStats() const
{
    std::cout << getID() << " stats" << std::endl;
    std::cout << "--------------------------------------" << std::endl;

    for(auto& keyPair : m_accumalatedTimes)
    {
        std::cout << std::defaultfloat << std::setprecision(7) <<  std::setw(40) << std::left << keyPair.first << ": " << std::setw(10) << std::right << keyPair.second << " s" << std::endl;
    }
}

SolTable Equation::getBCParam(unsigned int thread) const
{
    return m_bcParams[thread];
}

double Equation::getDiffusionParam(const Node& /** node **/) const
{
    throw std::runtime_error("Unimplemented function by the child class -> Equation::getDiffusionParam()!");
}

std::string Equation::getID() const noexcept
{
    return m_id;
}

bool Equation::isNormalCurvNeeded() const noexcept
{
    return m_needNormalCurv;
}

double Equation::getSquaredSpeedEquiv(const Node& /** node **/) const
{
    throw std::runtime_error("Unimplemented function by the child class -> Equation::getSquaredSpeedEquiv()!");
}

void Equation::preCompute()
{
    throw std::runtime_error("Unimplemented function by the child class -> Equation::preCompute()!");
}

bool Equation::solve()
{
    throw std::runtime_error("Unimplemented function by the child class -> Equation::solve()!");
}
