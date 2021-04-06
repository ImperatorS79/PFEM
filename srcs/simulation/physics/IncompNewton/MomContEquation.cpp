#include "MomContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

MomContEqIncompNewton::MomContEqIncompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "MomContEq")
{
    m_rho = m_materialParams[0].checkAndGet<double>("rho");
    m_mu = m_materialParams[0].checkAndGet<double>("mu");
    m_gamma = m_materialParams[0].checkAndGet<double>("gamma");

    if(m_pProblem->getID() == "Boussinesq")
    {
        m_alpha = m_materialParams[0].checkAndGet<double>("alpha");
        m_Tr = m_materialParams[0].checkAndGet<double>("Tr");

        if(statesIndex.size() != 2)
            throw std::runtime_error("the " + getID() + " equation require two statesIndex describing the beginning of the states span and the temperature state!");
    }
    else
    {
        if(statesIndex.size() != 1)
            throw std::runtime_error("the " + getID() + " equation require one statesIndex describing the beginning of the states span!");
    }

    if(bcFlags.size() != 1)
        throw std::runtime_error("the " + getID() + " equation require one flag for one possible boundary condition!");



    Eigen::MatrixXd ddev;
    if(m_pMesh->getDim() == 2)
    {
        ddev.resize(3, 3);
        ddev << 2, 0, 0,
                0, 2, 0,
                0, 0, 1;
    }
    else
    {
        ddev.resize(6, 6);
        ddev << 2, 0, 0, 0, 0, 0,
                0, 2, 0, 0, 0, 0,
                0, 0, 2, 0, 0, 0,
                0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1;
    }
    m_pMatBuilder->setddev(ddev);

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

    m_pMatBuilder->setMcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/) -> double {
        return m_rho;
    });

    m_pMatBuilder->setMGammacomputeFactor([&](const Facet& /** facet **/, const Eigen::MatrixXd& /** N **/) -> double {
        return m_gamma;
    });

    m_pMatBuilder->setKcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return m_mu;
    });

    m_pMatBuilder->setDcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
        return 1;
    });

    if(m_pSolver->getID() == "PSPG")
    {
        m_pMatBuilder->setCcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
            return 1;
        });

        m_pMatBuilder->setLcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
            return 1/m_rho;
        });
    }
    else if(m_pSolver->getID() == "FracStep")
    {
        m_gammaFS = m_equationParams[0].checkAndGet<double>("gammaFS");
        if(m_gammaFS < 0 || m_gammaFS > 1)
            throw std::runtime_error("gammaFS should be between 0 and 1");

        m_pMatBuilder->setLcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
            return 1;
        });
    }

    if(m_pProblem->getID() == "Boussinesq")
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& element, const Eigen::MatrixXd&  N, const Eigen::MatrixXd& /** B **/) -> double {
            double T = (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
            return m_rho*(1 - m_alpha*(T - m_Tr));
        });

        if(m_pSolver->getID() == "PSPG")
        {
            m_pMatBuilder->setHcomputeFactor([&](const Element& element, const Eigen::MatrixXd& N, const Eigen::MatrixXd& /** B **/) -> double {
                double T = (N*getElementState(m_pMesh, element, m_statesIndex[1])).value();
                return (1 - m_alpha*(T - m_Tr));
            });
        }
    }
    else
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
            return m_rho;
        });

        if(m_pSolver->getID() == "PSPG")
        {
            m_pMatBuilder->setHcomputeFactor([&](const Element& /** element **/, const Eigen::MatrixXd& /** N **/, const Eigen::MatrixXd& /** B **/) -> double {
                return 1;
            });
        }
    }

    unsigned int maxIter = m_equationParams[0].checkAndGet<unsigned int>("maxIter");
    double minRes = m_equationParams[0].checkAndGet<double>("minRes");

    auto bodyForce = m_equationParams[0].checkAndGet<std::vector<double>>("bodyForce");
    if(bodyForce.size() != m_pMesh->getDim())
        throw std::runtime_error("the body force vector has not the right dimension!");

    m_bodyForce = Eigen::Map<Eigen::VectorXd>(bodyForce.data(), bodyForce.size());

    if(m_pSolver->getID() == "PSPG")
    {
        m_setupPicardPSPG(maxIter, minRes);
    }
    else if(m_pSolver->getID() == "FracStep")
    {
        m_setupPicardFracStep(maxIter, minRes);
    }


    m_needNormalCurv = (m_gamma < 1e-15) ? false : true;
}

MomContEqIncompNewton::~MomContEqIncompNewton()
{

}

void MomContEqIncompNewton::displayParams() const
{
    std::cout << "Momentum Continuity equation parameters: \n"
              << " * Density: " << m_rho << " kg/m^3\n"
              << " * Viscosity: " << m_mu << " Pa s\n"
              << " * Surface Tension: " << m_gamma << " Nm" << std::endl;

    if(m_pProblem->getID() == "Boussinesq")
    {
        std::cout << " * Thermal expansion coefficient: " << m_alpha << " 1/K\n"
                  << " * Reference temperature: " << m_Tr << " K" << std::endl;
    }

    if(m_pMesh->getDim() == 2)
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << ")" << std::endl;
    else
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << "," << m_bodyForce[2] << ")" << std::endl;

    m_pPicardAlgo->displayParams();
}

bool MomContEqIncompNewton::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Momentum Continuity Equation" << std::endl;

    if(m_pSolver->getID() == "PSPG")
    {
        std::vector<Eigen::VectorXd> qPrev = {getQFromNodesStates(m_pMesh, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim())};
        return m_pPicardAlgo->solve(m_pMesh, qPrev, m_pProblem->isOutputVerbose());
    }
    else if(m_pSolver->getID() == "FracStep")
    {
        std::vector<Eigen::VectorXd> qPrev = {getQFromNodesStates(m_pMesh, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim() - 1),
                                              getQFromNodesStates(m_pMesh, m_pMesh->getDim(), m_pMesh->getDim())};
        return m_pPicardAlgo->solve(m_pMesh, qPrev, m_pProblem->isOutputVerbose());
    }

    return false;

}
