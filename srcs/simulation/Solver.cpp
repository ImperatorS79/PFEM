#include "Solver.hpp"

#include <iomanip>

#include "../mesh/Mesh.hpp"
#include "Problem.hpp"
#include "Equation.hpp"
#include "meshSmoother/PSmoother.hpp"
#include "meshSmoother/GETMe.hpp"


Solver::Solver(Problem* pProblem, Mesh* pMesh, std::vector<SolTable> m_problemParams):
m_timeStep(0),
m_pMesh(pMesh),
m_pProblem(pProblem)
{
    m_solverParams.resize(m_pProblem->getThreadCount());
    for(std::size_t i = 0 ; i < m_solverParams.size() ; ++i)
        m_solverParams[i] = SolTable("Solver", m_problemParams[i]);

    m_id = m_solverParams[0].checkAndGet<std::string>("id");

    m_adaptDT = m_solverParams[0].checkAndGet<bool>("adaptDT");
    m_maxDT = m_solverParams[0].checkAndGet<double>("maxDT");
//    m_maxRemeshDT = m_solverParams[0].checkAndGet<double>("maxRemeshDT");
    m_initialDT = m_solverParams[0].checkAndGet<double>("initialDT");

    if(m_solverParams[0].doesVarExist("MeshSmoother"))
    {
        SolTable meshSmoothers("MeshSmoother", m_solverParams[0]);
        meshSmoothers.for_each([&](sol::object /*key*/, sol::object value){
            SolTable smoother = SolTable(value);

            std::string kind = smoother.checkAndGet<std::string>("kind");

            if(kind == "PSmoother")
            {
                double a = smoother.checkAndGet<double>("a");
                double epsADRTol = smoother.checkAndGet<double>("epsADRTol");
                double beta = smoother.checkAndGet<double>("betaInit");

                if(m_pMesh->getDim() == 2)
                    m_pMeshSmoothers.push_back(std::make_unique<PSmoother<2>>(pProblem, *m_pMesh, a, epsADRTol, beta));
                else
                    m_pMeshSmoothers.push_back(std::make_unique<PSmoother<3>>(pProblem, *m_pMesh, a, epsADRTol, beta));
            }
            else if (kind == "GETMe")
            {
                double epsTol = smoother.checkAndGet<double>("epsTol");
                double maxIter = smoother.checkAndGet<double>("maxIter");

                if(m_pMesh->getDim() == 2)
                    m_pMeshSmoothers.push_back(std::make_unique<GETMe<2>>(pProblem, *m_pMesh, maxIter, epsTol));
                else
                    m_pMeshSmoothers.push_back(std::make_unique<GETMe<3>>(pProblem, *m_pMesh, maxIter, epsTol));
            }
            else
                throw std::runtime_error("unknown mesh smoother " + kind);

        });
    }

//    if(m_maxRemeshDT != -1)
//        m_nextTimeToRemesh = m_maxRemeshDT;
//    else
//        m_nextTimeToRemesh = m_maxDT;
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
                                        node.getPosition(), 0);

    if(res)
    {
        std::vector<double> result = bcParam.call<std::vector<double>>(
            m_pMesh->getNodeType(n) + bcString,
            node.getPosition(), 0);

        if(result.size() != expectedBCSize)
        {
            throw std::runtime_error(" the boundary condition " + m_pMesh->getNodeType(n) + bcString +
                                     " does not return the right number of states: " +
                                     std::to_string(result.size()) + " vs " +  std::to_string(expectedBCSize));
        }
    }

    return res;
}

bool Solver::checkFreeSurfaceBC(SolTable bcParam, const Node& node, std::string bcString, unsigned int expectedBCSize)
{
    bool res = bcParam.checkCallNoThrow("FreeSurface" + bcString,
                                        node.getPosition(), 0);

    if(res)
    {
        std::vector<double> result = bcParam.call<std::vector<double>>(
            "FreeSurface" + bcString,
            node.getPosition(), 0);

        if(result.size() != expectedBCSize)
        {
            throw std::runtime_error(" the boundary condition FreeSurface" + bcString +
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

void Solver::m_conditionalRemesh()
{
    bool force = false;
    if(m_pProblem->getCurrentSimTime() > m_nextTimeToRemesh)
    {
        force = true;
        m_pMesh->remesh(m_pProblem->isOutputVerbose());
        std::cout << "Remeshing at: " << std::fixed << m_pProblem->getCurrentSimTime() << std::endl;
    }

    if(force || m_pProblem->getCurrentSimTime() + m_remeshTimeStep < m_nextTimeToRemesh)
        m_nextTimeToRemesh = m_pProblem->getCurrentSimTime() + m_remeshTimeStep;
}

