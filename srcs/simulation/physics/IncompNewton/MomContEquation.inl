#include "MomContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"

template<unsigned short dim>
MomContEqIncompNewton<dim>::MomContEqIncompNewton(Problem* pProblem, Solver* pSolver, Mesh* pMesh,
                                     std::vector<SolTable> solverParams, std::vector<SolTable> materialParams,
                                     const std::vector<unsigned short>& bcFlags, const std::vector<unsigned int>& statesIndex) :
Equation(pProblem, pSolver, pMesh, solverParams, materialParams, bcFlags, statesIndex, "MomContEq")
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

    m_rho = m_materialParams[0].checkAndGet<double>("rho");
    m_mu = m_materialParams[0].checkAndGet<double>("mu");
    m_gamma = m_materialParams[0].checkAndGet<double>("gamma");

    m_phaseChange = false;

    if(m_pProblem->getID() == "Boussinesq")
    {
        m_alpha = m_materialParams[0].checkAndGet<double>("alpha");
        m_Tr = m_materialParams[0].checkAndGet<double>("Tr");
        m_DgammaDT = m_materialParams[0].checkAndGet<double>("DgammaDT");

        if(m_materialParams[0].doesVarExist("Tm")  || m_materialParams[0].doesVarExist("C") ||
           m_materialParams[0].doesVarExist("eps") || m_materialParams[0].doesVarExist("DT"))
        {
            m_phaseChange = true;

            m_C = m_materialParams[0].doesVarExist("C");
            m_eps = m_materialParams[0].doesVarExist("eps");
            m_Tm = m_materialParams[0].doesVarExist("Tm");
            m_DT = m_materialParams[0].doesVarExist("DT");
        }

        if(statesIndex.size() != 2)
            throw std::runtime_error("the " + getID() + " equation require two statesIndex describing the beginning of the states span and the temperature state!");
    }
    else if(m_pProblem->getID() == "Bingham")
    {
        m_tau0 = m_materialParams[0].checkAndGet<double>("tau0");
        m_mReg = m_materialParams[0].checkAndGet<double>("mReg");

        if(statesIndex.size() != 1)
            throw std::runtime_error("the " + getID() + " equation require one statesIndex describing the beginning of the states span!");
    }
    else
    {
        if(statesIndex.size() != 1)
            throw std::runtime_error("the " + getID() + " equation require one statesIndex describing the beginning of the states span!");
    }

    if(bcFlags.size() != 1)
        throw std::runtime_error("the " + getID() + " equation require one flag for one possible boundary condition!");



    DdevMatType<dim> ddev;
    mVecType<dim> m;
    if constexpr (dim == 2)
    {
        ddev << 2, 0, 0,
                0, 2, 0,
                0, 0, 1;

        m << 1, 1, 0;
    }
    else if constexpr (dim == 3)
    {
        ddev << 2, 0, 0, 0, 0, 0,
                0, 2, 0, 0, 0, 0,
                0, 0, 2, 0, 0, 0,
                0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1;

        m << 1, 1, 1, 0, 0, 0;
    }
    m_pMatBuilder->setddev(ddev);
    m_pMatBuilder->setm(m);

    m_pMatBuilder->setMcomputeFactor([&](const Element& /** element **/,
                                         const NmatTypeHD<dim>& /** N **/) -> double {
        return m_rho;
    });

    if(m_pProblem->getID() == "Bingham")
    {
        m_pMatBuilder->setKcomputeFactor([&](const Element& element,
                                             const NmatTypeHD<dim>& /** N **/,
                                             const BmatType<dim>& B,
                                             const DdevMatType<dim>& ddevParam) -> double {

            Eigen::Matrix<double, dim*(dim + 1), 1> V = getElementVecState<dim>(element,  m_statesIndex[0]);

            double gammaDot = std::sqrt(V.transpose()*B.transpose()*ddevParam*B*V);
            double muEq = m_tau0;
            if(gammaDot < 1e-15)
                muEq *= m_mReg;
            else
                muEq *=(1 - std::exp(- m_mReg*gammaDot))/gammaDot;

            return m_mu + muEq;
        });
    }
    else
    {
        m_pMatBuilder->setKcomputeFactor([&](const Element& /** element **/,
                                             const NmatTypeHD<dim>& /** N **/,
                                             const BmatType<dim>& /** B **/,
                                             const DdevMatType<dim>& /** ddev **/) -> double {
            return m_mu;
        });
    }

    m_pMatBuilder->setDcomputeFactor([&](const Element& /** element **/,
                                         const NmatTypeHD<dim>& /** N **/,
                                         const BmatType<dim>& /** B **/) -> double {
        return 1;
    });

    if(m_pSolver->getID() == "PSPG")
    {
        m_pMatBuilder->setCcomputeFactor([&](const Element& /** element **/,
                                             const NmatTypeHD<dim>& /** N **/,
                                             const BmatType<dim>& /** B **/) -> double {
            return 1;
        });

        m_pMatBuilder->setLcomputeFactor([&](const Element& /** element **/,
                                             const NmatTypeHD<dim>& /** N **/,
                                             const BmatType<dim>& /** B **/) -> double {
            return 1/m_rho;
        });
    }
    else if(m_pSolver->getID() == "FracStep")
    {
        m_gammaFS = m_equationParams[0].checkAndGet<double>("gammaFS");
        if(m_gammaFS < 0 || m_gammaFS > 1)
            throw std::runtime_error("gammaFS should be between 0 and 1");

        m_pMatBuilder->setLcomputeFactor([&](const Element& /** element **/,
                                             const NmatTypeHD<dim>& /** N **/,
                                             const BmatType<dim>& /** B **/) -> double {
            return 1;
        });
    }

    if(m_pProblem->getID() == "Boussinesq")
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& element,
                                             const NmatTypeHD<dim>&  N,
                                             const BmatType<dim>& /** B **/) -> double {
            double T = (N*getElementState<dim>(element, m_statesIndex[1])).value();
            return m_rho*(1 - m_alpha*(T - m_Tr));
        });

        m_pMatBuilder->setFSTcomputeFactor([&](const Facet& facet,
                                               const NmatTypeLD<dim>& N,
                                               const NmatTildeTypeLD<dim>&  /** Ntilde **/,
                                               const GradNmatType<dim>& /** gradNe **/) -> double {
            double T = (N*getFacetState<dim>(facet, m_statesIndex[1])).value();
            return m_gamma + m_DgammaDT*(T - m_Tr);
        });

        if(m_phaseChange)
        {
            m_pMatBuilder2->setMcomputeFactor([&](const Element& element,
                                                  const NmatTypeHD<dim>& N) -> double {
                double T = (N*getElementState<dim>(element, m_statesIndex[1])).value();
                double fl = m_getFl(T);
                return m_C*(1 - fl)*(1 - fl)/(fl*fl*fl + m_eps);
            });
        }

        if(m_pSolver->getID() == "PSPG")
        {
            m_pMatBuilder->setHcomputeFactor([&](const Element& element,
                                                 const NmatTypeHD<dim>& N,
                                                 const BmatType<dim>& /** B **/) -> double {
                double T = (N*getElementState<dim>(element, m_statesIndex[1])).value();
                return (1 - m_alpha*(T - m_Tr));
            });
        }
    }
    else
    {
        m_pMatBuilder->setFcomputeFactor([&](const Element& /** element **/,
                                             const NmatTypeHD<dim>& /** N **/,
                                             const BmatType<dim>& /** B **/) -> double {
            return m_rho;
        });

        m_pMatBuilder->setFSTcomputeFactor([&](const Facet& /** facet **/,
                                               const NmatTypeLD<dim>& /** N **/,
                                               const NmatTildeTypeLD<dim>&  /** Ntilde **/,
                                               const GradNmatType<dim>& /** gradNe **/) -> double {
            return m_gamma;
        });

        if(m_pSolver->getID() == "PSPG")
        {
            m_pMatBuilder->setHcomputeFactor([&](const Element& /** element **/,
                                                 const NmatTypeHD<dim>& /** N **/,
                                                 const BmatType<dim>& /** B **/) -> double {
                return 1;
            });
        }
    }

    unsigned int maxIter = m_equationParams[0].checkAndGet<unsigned int>("maxIter");
    double minRes = m_equationParams[0].checkAndGet<double>("minRes");
    std::string residual = m_equationParams[0].checkAndGet<std::string>("residual");
    if(residual == "U")
        m_residual = Res::U;
    else if(residual == "U_P")
        m_residual = Res::U_P;
    else if(residual == "Ax_f")
        m_residual = Res::Ax_f;
    else
        throw std::runtime_error("unknown residual type: " + residual);


    auto bodyForce = m_equationParams[0].checkAndGet<std::vector<double>>("bodyForce");
    if(bodyForce.size() != m_pMesh->getDim())
        throw std::runtime_error("the body force vector has not the right dimension!");

    m_bodyForce = Eigen::Map<Eigen::Matrix<double, dim, 1>>(bodyForce.data(), bodyForce.size());

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

template<unsigned short dim>
MomContEqIncompNewton<dim>::~MomContEqIncompNewton()
{

}

template<unsigned short dim>
void MomContEqIncompNewton<dim>::displayParams() const
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

    if constexpr (dim == 2)
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << ")" << std::endl;
    else if constexpr (dim == 3)
        std::cout << " * Body force: (" << m_bodyForce[0] << ", " << m_bodyForce[1] << "," << m_bodyForce[2] << ")" << std::endl;

    m_pPicardAlgo->displayParams();
}

template<unsigned short dim>
bool MomContEqIncompNewton<dim>::solve()
{
    if(m_pProblem->isOutputVerbose())
        std::cout << "Momentum Continuity Equation" << std::endl;

    if(m_pSolver->getID() == "PSPG")
    {
        m_clock.start();
        std::vector<Eigen::VectorXd> qPrev = {getQFromNodesStates(m_pMesh, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim())};
        m_accumalatedTimes["Update solutions"] += m_clock.end();
        return m_pPicardAlgo->solve(m_pMesh, qPrev, m_pProblem->isOutputVerbose());
    }
    else if(m_pSolver->getID() == "FracStep")
    {
        m_clock.start();
        std::vector<Eigen::VectorXd> qPrev = {getQFromNodesStates(m_pMesh, m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim() - 1),
                                              getQFromNodesStates(m_pMesh, m_pMesh->getDim(), m_pMesh->getDim())};
        m_accumalatedTimes["Update solutions"] += m_clock.end();
        return m_pPicardAlgo->solve(m_pMesh, qPrev, m_pProblem->isOutputVerbose());
    }

    return false;
}

template<unsigned short dim>
double MomContEqIncompNewton<dim>::m_getFl(double T)
{
    return static_cast<double>(T > m_Tm - m_DT/2)*static_cast<double>(T < m_Tm + m_DT/2)*(-2/(m_DT*m_DT*m_DT)*T*T*T + 6*m_Tm/(m_DT*m_DT*m_DT)*T*T
           + (3/(2*m_DT) - 6*m_Tm*m_Tm/(m_DT*m_DT*m_DT))*T + (0.5 - 1.5*m_Tm/m_DT + 2*m_Tm*m_Tm*m_Tm/(m_DT*m_DT*m_DT))) + static_cast<double>(T > m_Tm + m_DT/2);
}
