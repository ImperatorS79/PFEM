#include "Problem.hpp"
#include "Solver.hpp"
#include "MomEquation.hpp"
#include "ContEquation.hpp"
#include "HeatEquation.hpp"
#include "../../utility/StatesFromToQ.hpp"
#include "../../utility/Clock.hpp"

#define REGISTER_EQ(Eq, dim) \
std::make_unique<Eq<dim>>( \
    m_pProblem, this, m_pMesh, m_solverParams, materialParams, \
    bcFlags, statesIndex \
); \

SolverWCompNewton::SolverWCompNewton(Problem* pProblem, Mesh* pMesh, std::vector<SolTable> problemParams):
Solver(pProblem, pMesh, problemParams)
{
	//Check if the asked problem and solver are supported
    if(m_pProblem->getID() != "WCompNewtonNoT" && m_pProblem->getID() != "BoussinesqWC")
        throw std::runtime_error("this solver cannot be used with problem whose id is " + m_pProblem->getID());

    if(m_id != "CDS_Meduri" && m_id != "CDS_FIC")
        throw std::runtime_error("this solver does not know id " + m_id);

    //Load material params for equations
    std::vector<SolTable> materialParams(problemParams.size());
    for(std::size_t i = 0 ; i < problemParams.size() ; ++i)
    {
        materialParams[i] = SolTable("Material", problemParams[i]);
    }

    unsigned int dim = m_pMesh->getDim();

	//Load equations depending of the problem
    std::vector<unsigned short> bcFlags;
    std::vector<unsigned int> statesIndex;
    if(m_pProblem->getID() == "WCompNewtonNoT")
    {
        //Momentum and continuity equations
        m_pEquations.resize(2);

		bcFlags = {0};
		statesIndex = {dim, dim + 1, 0};
		if(m_pMesh->getDim() == 2)
            m_pEquations[0] = REGISTER_EQ(ContEqWCompNewton, 2)
        else
            m_pEquations[0] = REGISTER_EQ(ContEqWCompNewton, 3)

		bcFlags = {0};
		statesIndex = {0, dim + 2, dim, dim + 1};
		if(m_pMesh->getDim() == 2)
            m_pEquations[1] = REGISTER_EQ(MomEqWCompNewton, 2)
        else
            m_pEquations[1] = REGISTER_EQ(MomEqWCompNewton, 3)

        //Set the right node flag if the boundary condition is present
        SolTable bcParam = m_pEquations[1]->getBCParam(0);
        for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
        {
            const Node& node = m_pMesh->getNode(n);
            if(node.isBound())
            {
                bool res = checkBC(bcParam, n, node, "V", m_pMesh->getDim());

                if(res)
                    m_pMesh->setNodeFlag(n, 0);
            }
        }

        m_solveFunc = std::bind(&SolverWCompNewton::m_solveWCompNewtonNoT, this);
    }
    else if(m_pProblem->getID() == "BoussinesqWC")
    {
        //The momentum, continuity and the heat equation
        m_pEquations.resize(3);

        bcFlags = {0};
		statesIndex = {dim, dim + 1, 0, 2*dim + 2};
		if(m_pMesh->getDim() == 2)
            m_pEquations[0] = REGISTER_EQ(ContEqWCompNewton, 2)
        else
            m_pEquations[0] = REGISTER_EQ(ContEqWCompNewton, 3)

		bcFlags = {0};
		statesIndex = {0, dim + 2, dim, dim + 1, 2*dim + 2};
		if(m_pMesh->getDim() == 2)
            m_pEquations[1] = REGISTER_EQ(MomEqWCompNewton, 2)
        else
            m_pEquations[1] = REGISTER_EQ(MomEqWCompNewton, 3)

        bcFlags = {1, 2}; //Dirichlet and Neumann
        statesIndex = {2*dim + 2, dim + 1};
        if(m_pMesh->getDim() == 2)
            m_pEquations[2] = REGISTER_EQ(HeatEqWCompNewton, 2)
        else
            m_pEquations[2] = REGISTER_EQ(HeatEqWCompNewton, 3)

        //Set the right node flag if the boundary condition is present
        SolTable bcParamMomCont = m_pEquations[1]->getBCParam(0);
        SolTable bcParamHeat = m_pEquations[2]->getBCParam(0);
        for(std::size_t n = 0 ; n < m_pMesh->getNodesCount() ; ++n)
        {
            const Node& node = m_pMesh->getNode(n);
            if(node.isBound())
            {
                bool resV = checkBC(bcParamMomCont, n, node, "V", m_pMesh->getDim());
                bool resT = checkBC(bcParamHeat, n, node, "T", 1);
                bool resQ = checkBC(bcParamHeat, n, node, "Q", 1);

                if(resV)
                    m_pMesh->setNodeFlag(n, 0);

                if(resT)
                    m_pMesh->setNodeFlag(n, 1);

                if(resQ)
                    m_pMesh->setNodeFlag(n, 2);

                if(resT && resQ)
                    throw std::runtime_error("the boundary " + m_pMesh->getNodeType(n) +
                                             "has a BC for both T and Q. This is forbidden!");
            }
        }

        m_solveFunc = std::bind(&SolverWCompNewton::m_solveBoussinesqWC, this);
	}

    m_securityCoeff = m_solverParams[0].checkAndGet<double>("securityCoeff");

    m_timeStep = m_initialDT;

    //Should we compute the normals and the curvature ?
    m_pMesh->setComputeNormalCurvature(false);
    for(auto& pEq : m_pEquations)
    {
        if(pEq->isNormalCurvNeeded())
        {
            m_pMesh->setComputeNormalCurvature(true);
            break;
        }
    }
}

SolverWCompNewton::~SolverWCompNewton()
{

}

void SolverWCompNewton::displayParams() const
{
    std::cout << "Maximum dt: " << m_maxDT << "\n"
              << "Initial dt: " << m_initialDT << "\n"
              << "Security coeff: " << m_securityCoeff << std::endl;

    for(auto& pEquation : m_pEquations)
        pEquation->displayParams();
}

std::size_t SolverWCompNewton::getAdditionalStateCount() const
{
    return (m_id == "CDS_FIC") ? 1 : 0;
}

bool SolverWCompNewton::solveOneTimeStep()
{
    return m_solveFunc();
}

void SolverWCompNewton::computeNextDT()
{
    if(m_adaptDT)
    {
        bool heat = false;
        if(m_pProblem->getID() == "BoussinesqWC")
            heat = true;

        m_timeStep = std::numeric_limits<double>::max();
        m_remeshTimeStep = std::numeric_limits<double>::max();
        #pragma omp parallel for reduction(min:m_timeStep)
        for(std::size_t elm = 0 ; elm < m_pMesh->getElementsCount() ; ++elm)
        {
            const Element& element = m_pMesh->getElement(elm);
            double he = element.getRin();

            double maxSquaredSpeedFluidPressureVN = 0;
            double maxSquaredSpeedFluid = 0;
            for(std::size_t n = 0 ; n < m_pMesh->getNodesPerElm() ; ++n)
            {
                const Node& node = element.getNode(n);

                double c2 = m_pEquations[0]->getSquaredSpeedEquiv(node);
                double u2 = m_pEquations[1]->getSquaredSpeedEquiv(node);
                double alpha = m_pEquations[1]->getDiffusionParam(node);
                if(heat)
                    alpha = std::max(alpha, m_pEquations[2]->getDiffusionParam(node));

                maxSquaredSpeedFluidPressureVN = std::max(std::max(u2, c2), maxSquaredSpeedFluidPressureVN);
                maxSquaredSpeedFluidPressureVN = std::max(maxSquaredSpeedFluidPressureVN, 4*alpha*alpha/(he*he));
                maxSquaredSpeedFluid = std::max(maxSquaredSpeedFluid, u2);
            }
            m_timeStep = std::min(m_timeStep, m_securityCoeff*m_securityCoeff*he*he/maxSquaredSpeedFluidPressureVN);
            m_remeshTimeStep = std::max(m_remeshTimeStep, maxSquaredSpeedFluid);
        }

        m_timeStep = std::min(std::sqrt(m_timeStep), m_maxDT);
        m_remeshTimeStep = std::min(std::sqrt(m_remeshTimeStep), m_maxRemeshDT);

        if(std::isnan(m_timeStep) || std::isnan(m_remeshTimeStep))
            throw std::runtime_error("NaN time step!");
    }
}

bool SolverWCompNewton::m_solveWCompNewtonNoT()
{
    m_clock.start();
    unsigned int dim = m_pMesh->getDim();
    Eigen::VectorXd qVPrev = getQFromNodesStates(m_pMesh, 0, dim - 1);           //The precedent speed.
    Eigen::VectorXd qAccPrev = getQFromNodesStates(m_pMesh, dim + 2, 2*dim + 1);  //The precedent acceleration.
    m_accumalatedTimes["Update solutions"] += m_clock.end();

    m_clock.start();
    m_pEquations[0]->preCompute();
    m_accumalatedTimes["Solving continuity eq"] += m_clock.end();

    m_clock.start();
    Eigen::VectorXd qV1half = qVPrev + 0.5*m_timeStep*qAccPrev;

    setNodesStatesfromQ(m_pMesh, qV1half, 0, dim - 1);
    Eigen::VectorXd deltaPos = qV1half*m_timeStep;
    m_pMesh->updateNodesPosition(deltaPos);
    m_accumalatedTimes["Update solutions"] += m_clock.end();

    m_clock.start();
    m_pEquations[0]->solve();
    m_accumalatedTimes["Solving continuity eq"] += m_clock.end();
    m_clock.start();
    m_pEquations[1]->solve();
    m_accumalatedTimes["Solving momentum eq"] += m_clock.end();

    m_clock.start();
    m_pProblem->updateTime(m_timeStep);
    if(m_pProblem->getCurrentSimTime() > m_nextTimeToRemesh)
    {
        m_pMesh->remesh(m_pProblem->isOutputVerbose());
        m_nextTimeToRemesh += m_maxDT;
    }
    m_accumalatedTimes["Remeshing"] += m_clock.end();

    return true;
}

bool SolverWCompNewton::m_solveBoussinesqWC()
{
    m_clock.start();
    unsigned int dim = m_pMesh->getDim();
    Eigen::VectorXd qVPrev = getQFromNodesStates(m_pMesh, 0, dim - 1);           //The precedent speed.
    Eigen::VectorXd qAccPrev = getQFromNodesStates(m_pMesh, dim + 2, 2*dim + 1);  //The precedent acceleration.
    m_accumalatedTimes["Update solutions"] += m_clock.end();

    m_clock.start();
    m_pEquations[0]->preCompute();
    m_accumalatedTimes["Solving continuity eq"] += m_clock.end();

    m_clock.start();
    Eigen::VectorXd qV1half = qVPrev + 0.5*m_timeStep*qAccPrev;

    setNodesStatesfromQ(m_pMesh, qV1half, 0, dim - 1);
    Eigen::VectorXd deltaPos = qV1half*m_timeStep;
    m_pMesh->updateNodesPosition(std::vector<double> (deltaPos.data(), deltaPos.data() + deltaPos.cols()*deltaPos.rows()));
    m_accumalatedTimes["Update solutions"] += m_clock.end();

    m_clock.start();
    m_pEquations[2]->solve();
    m_accumalatedTimes["Solving heat eq"] += m_clock.end();
    m_clock.start();
    m_pEquations[0]->solve();
    m_accumalatedTimes["Solving continuity eq"] += m_clock.end();
    m_clock.start();
    m_pEquations[1]->solve();
    m_accumalatedTimes["Solving momentum eq"] += m_clock.end();

    m_clock.start();
    m_pProblem->updateTime(m_timeStep);
    if(m_pProblem->getCurrentSimTime() > m_nextTimeToRemesh)
    {
        m_pMesh->remesh(m_pProblem->isOutputVerbose());
        m_nextTimeToRemesh += m_maxDT;
    }
    m_accumalatedTimes["Remeshing"] += m_clock.end();

    return true;
}
