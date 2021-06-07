#include "MomContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"


template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_buildMatFracStep(const std::vector<Eigen::VectorXd>& qPrev)
{
    m_clock.start();
    constexpr unsigned short nodPerEl = dim + 1;
    constexpr unsigned int tripletPerElmM = (dim*nodPerEl*nodPerEl);
    constexpr unsigned int tripletPerElmMKDT = (dim*nodPerEl*nodPerEl + dim*nodPerEl*dim*nodPerEl);
    constexpr unsigned int tripletPerElmDtL = (nodPerEl*nodPerEl);
    constexpr unsigned int doubletPerElm = (3*dim*nodPerEl);

    const std::size_t nElm = m_pMesh->getElementsCount();
    const std::size_t nNodes = m_pMesh->getNodesCount();
    const double dt = m_pSolver->getTimeStep();

    std::vector<Eigen::Triplet<double>> indexM(tripletPerElmM*nElm);
    std::vector<Eigen::Triplet<double>> indexMK_dt(tripletPerElmMKDT*nElm);
    std::vector<Eigen::Triplet<double>> indexL(tripletPerElmDtL*nElm);
    std::vector<std::pair<std::size_t, double>> indexbVappStep(doubletPerElm*nElm); m_bVAppStep.setZero();
    m_DTelm.resize(nElm);
    m_Lelm.resize(nElm);
    m_accumalatedTimes["Prepare matrices assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        GradNmatType<dim> gradNe = m_pMatBuilder->getGradN(element);
        BmatType<dim> Be = m_pMatBuilder->getB(gradNe);
        auto Me_s = m_pMatBuilder->getM(element);
        auto Me_dt_s = static_cast<Eigen::Matrix<double, nodPerEl, nodPerEl>>((1/dt)*Me_s);
        auto Me = MatrixBuilder<dim>::diagBlock(Me_s);
        auto Me_dt = MatrixBuilder<dim>::diagBlock(Me_dt_s);
        auto Ke = m_pMatBuilder->getK(element, Be);
        m_DTelm[elm] = m_pMatBuilder->getD(element, Be).transpose();
        m_Lelm[elm] = m_pMatBuilder->getL(element, Be, gradNe);
        auto Fe = m_pMatBuilder->getF(element, m_bodyForce, Be);

        auto vPrev = getElementVecState<dim>(qPrev[0], element, 0, nNodes);
        auto pPrev = getElementState<dim>(qPrev[1], element, 0, nNodes);

        auto MvPreve_dt = Me_dt*vPrev;
        auto DTpPreve = m_gammaFS*m_DTelm[elm]*pPrev;

        std::size_t countM = 0;
        std::size_t countMK_dt = 0;
        std::size_t countL = 0;
        std::size_t countbVappStep = 0;

        for(unsigned short i = 0 ; i < nodPerEl ; ++i)
        {
            const Node& ni = m_pMesh->getNode(element.getNodeIndex(i));

            for(unsigned short j = 0 ; j < nodPerEl ; ++j)
            {
                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    /********************************************************************
                                                 Build M and M/dt
                    ********************************************************************/
                    if(!(ni.isBound() || ni.isFree()))
                    {
                        indexM[tripletPerElmM*elm + countM] =
                            Eigen::Triplet<double>(element.getNodeIndex(i) + d*nNodes,
                                                   element.getNodeIndex(j) + d*nNodes,
                                                   Me(i + d*nodPerEl, j + d*nodPerEl));

                        indexMK_dt[tripletPerElmMKDT*elm + countMK_dt] =
                            Eigen::Triplet<double>(element.getNodeIndex(i) + d*nNodes,
                                                   element.getNodeIndex(j) + d*nNodes,
                                                   Me_dt(i + d*nodPerEl, j + d*nodPerEl));
                    }

                    countM++;
                    countMK_dt++;

                    /********************************************************************
                                                  Build K
                    ********************************************************************/
                    for(unsigned short d2 = 0 ; d2 < dim ; ++d2)
                    {

                        if(!(ni.isBound() || ni.isFree()))
                        {
                            indexMK_dt[tripletPerElmMKDT*elm + countMK_dt] =
                                Eigen::Triplet<double>(element.getNodeIndex(i) + d*nNodes,
                                                       element.getNodeIndex(j) + d2*nNodes,
                                                       Ke(i + d*nodPerEl, j + d2*nodPerEl));
                        }
                        countMK_dt++;
                    }
                }

                /********************************************************************
                                            Build L
                ********************************************************************/
                if(!ni.isFree() && !ni.isOnFreeSurface())
                {
                    indexL[tripletPerElmDtL*elm + countL] =
                        Eigen::Triplet<double>(element.getNodeIndex(i),
                                               element.getNodeIndex(j),
                                               m_Lelm[elm](i,j));
                }
                countL++;
            }

            /************************************************************************
                                              Build f
            ************************************************************************/
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                if(!(ni.isBound() || ni.isFree()))
                {
                    indexbVappStep[doubletPerElm*elm + countbVappStep] = std::make_pair(element.getNodeIndex(i) + d*nNodes, Fe(i + d*nodPerEl));
                    countbVappStep++;

                    indexbVappStep[doubletPerElm*elm + countbVappStep] = std::make_pair(element.getNodeIndex(i) + d*nNodes, MvPreve_dt(i + d*nodPerEl));
                    countbVappStep++;

                    indexbVappStep[doubletPerElm*elm + countbVappStep] = std::make_pair(element.getNodeIndex(i) + d*nNodes, DTpPreve(i + d*nodPerEl));
                    countbVappStep++;
                }
            }
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    //Best would be to know the number of nodes in which case :/
    //This can still be fasten using OpenMP but will never be as good as using []
    //with preallocated memory
    m_clock.start();
    for(std::size_t n = 0 ; n < nNodes ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree() || node.isOnFreeSurface())
        {
            indexL.push_back(Eigen::Triplet<double>(n, n, 1));
        }

        if(node.isBound() || node.isFree())
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                indexMK_dt.push_back(Eigen::Triplet<double>(n + d*nNodes,
                                                            n + d*nNodes,
                                                            1));

                indexM.push_back(Eigen::Triplet<double>(n + d*nNodes,
                                                        n + d*nNodes,
                                                        1));
            }
        }
    }
    m_accumalatedTimes["Push back (n, n, 1)"] += m_clock.end();

    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    m_clock.start();
    m_M.setFromTriplets(indexM.begin(), indexM.end());
    m_MK_dt.setFromTriplets(indexMK_dt.begin(), indexMK_dt.end());
    m_L.setFromTriplets(indexL.begin(), indexL.end());
    m_accumalatedTimes["Assemble matrices"] += m_clock.end();

    m_clock.start();
    for(const auto& doublet : indexbVappStep)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        m_bVAppStep[doublet.first] += doublet.second;
    }
    m_accumalatedTimes["Assemble b v app step"] += m_clock.end();
}

template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_applyBCVAppStep(const Eigen::VectorXd& qPrev)
{
    const std::size_t nodesCount = m_pMesh->getNodesCount();
    const std::size_t facetsCount = m_pMesh->getFacetsCount();
    constexpr unsigned short noPerFacet = dim;

    if(m_gamma < 1e-15)
        goto applyBC; //Heresy ^^

    for(std::size_t f = 0 ; f < facetsCount ; ++f)
    {
        const Facet& facet = m_pMesh->getFacet(f);

        bool onFS = true;
        for(unsigned short n = 0 ; n < noPerFacet ; ++n)
        {
            if(!facet.getNode(n).isOnFreeSurface())
            {
                onFS = false;
                break;
            }
        }
        if(!onFS)
            continue;

        auto MGamma_s = m_pMatBuilder->getMGamma(facet);
        auto MGamma = MatrixBuilder<dim>::diagBlock(MGamma_s);

        Eigen::Matrix<double, dim*noPerFacet, 1> nVec;
        for(uint8_t n = 0 ; n < noPerFacet; ++n)
        {
            std::array<double, 3> normal = m_pMesh->getBoundFSNormal(facet.getNodeIndex(n));

            for(uint8_t d = 0 ; d < dim ; ++d)
                nVec(n + d*noPerFacet) = normal[d];
        }

        Eigen::Matrix<double, dim*noPerFacet, 1> Ff = MGamma*nVec;

        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
                m_bVAppStep(facet.getNodeIndex(i) + d*nodesCount) += Ff[d*noPerFacet + i];
        }
    }

    //Do not parallelize this (lua)
    applyBC:
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree())
        {
            if(!node.isBound())
            {
                for(uint8_t d = 0 ; d < dim  ; ++d)
                {
                    m_bVAppStep(n + d*nodesCount) = qPrev(n + d*nodesCount) + m_pSolver->getTimeStep()*m_bodyForce[d];
                }
            }
        }

        if(node.isBound())
        {
            if(node.getFlag(m_bcFlags[0]))
            {
                std::array<double, dim> result; //TO do: try to change that
                result = m_bcParams[0].call<std::array<double, dim>>(m_pMesh->getNodeType(n) + "V",
                                                        node.getPosition(),
                                                        m_pMesh->getBoundNodeInitPos(n),
                                                        m_pProblem->getCurrentSimTime() +
                                                        m_pSolver->getTimeStep());

                for(uint8_t d = 0 ; d < dim ; ++d)
                {
                    m_bVAppStep(n + d*nodesCount) = result[d];
                    for(Eigen::SparseMatrix<double>::InnerIterator it(m_MK_dt, n + d*nodesCount); it; ++it)
                    {
                        Eigen::Index row = it.row();
                        if(row == it.col())
                            continue;

                        double value = it.value();
                        m_bVAppStep(row) -= value*result[d];
                        it.valueRef() = 0;
                    }
                }
            }
        }
    }

    m_MK_dt.makeCompressed();
}

template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_buildMatPcorrStep(const Eigen::VectorXd& qVTilde, const Eigen::VectorXd& qPprev)
{
    m_clock.start();
    constexpr unsigned short nodPerEl = dim + 1;
    constexpr unsigned int doubletPerElm = 2*nodPerEl;

    const std::size_t nElm = m_pMesh->getElementsCount();
    const std::size_t nNodes = m_pMesh->getNodesCount();

    const double dt = m_pSolver->getTimeStep();

    std::vector<std::pair<std::size_t, double>> indexbPcorrStep(doubletPerElm*nElm); m_bPcorrStep.setZero();
    m_accumalatedTimes["Prepare matrices assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        auto vTilde = getElementVecState<dim>(qVTilde, element, 0, nNodes);
        auto pPrev = getElementState<dim>(qPprev, element, 0, nNodes);

        auto rho_dt_DvTilde = (m_rho/dt)*m_DTelm[elm].transpose()*vTilde;
        auto LpPrev = m_gammaFS*m_Lelm[elm]*pPrev;

        std::size_t countbPcorrStep = 0;
        for(unsigned short i = 0 ; i < nodPerEl ; ++i)
        {
            const Node& node = element.getNode(i);
            if(!node.isFree() && !node.isOnFreeSurface())
            {
                /************************************************************************
                                                  Build f
                ************************************************************************/
                indexbPcorrStep[doubletPerElm*elm + countbPcorrStep] = std::make_pair(element.getNodeIndex(i), -rho_dt_DvTilde(i));
                countbPcorrStep++;

                indexbPcorrStep[doubletPerElm*elm + countbPcorrStep] = std::make_pair(element.getNodeIndex(i), LpPrev(i));
                countbPcorrStep++;
            }
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    m_clock.start();
    for(const auto& doublet : indexbPcorrStep)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        m_bPcorrStep[doublet.first] += doublet.second;
    }
    m_accumalatedTimes["Assemble p corr step"] += m_clock.end();
}

template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_applyBCPCorrStep()
{
    std::size_t nodesCount = m_pMesh->getNodesCount();

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree() || node.isOnFreeSurface())
        {
            m_bPcorrStep[n] = 0;

            for(Eigen::SparseMatrix<double>::InnerIterator it(m_L, n); it; ++it)
            {
                Eigen::Index row = it.row();
                if(row == it.col())
                {
                    it.valueRef() = 1;
                    continue;
                }

                it.valueRef() = 0;
            }
        }
    }

    m_L.makeCompressed();
}

template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_buildMatVStep(const Eigen::VectorXd& qDeltaP)
{
    m_clock.start();
    constexpr unsigned short nodPerEl = dim + 1;
    const unsigned int doubletPerElm = dim*nodPerEl;
    const double dt = m_pSolver->getTimeStep();

    const std::size_t nElm = m_pMesh->getElementsCount();
    const std::size_t nNodes = m_pMesh->getNodesCount();

    std::vector<std::pair<std::size_t, double>> indexbVStep(doubletPerElm*nElm); m_bVStep.setZero();
    m_accumalatedTimes["Prepare matrices assembly"] += m_clock.end();

    m_clock.start();
    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        auto deltaP = getElementState<dim>(qDeltaP, element, 0, nNodes);
        auto DTdeltaP = m_DTelm[elm]*deltaP;

        std::size_t countbVStep = 0;
        for(unsigned short i = 0 ; i < nodPerEl ; ++i)
        {
            const Node& node = element.getNode(i);

            if(!node.isFree() && !node.isBound())
            {
                /************************************************************************
                                                  Build f
                ************************************************************************/
                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    indexbVStep[doubletPerElm*elm + countbVStep] = std::make_pair(element.getNodeIndex(i) + d*nNodes, dt*DTdeltaP(i + d*nodPerEl));
                    countbVStep++;
                }
            }
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());
    m_accumalatedTimes["Compute triplets"] += m_clock.end();

    m_clock.start();
    for(const auto& doublet : indexbVStep)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        m_bVStep[doublet.first] += doublet.second;
    }
    m_accumalatedTimes["Assemble b v corr step"] += m_clock.end();
}

template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_applyBCVStep()
{
    const std::size_t nodesCount = m_pMesh->getNodesCount();

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree() || node.isBound())
        {
            for(uint8_t d = 0 ; d < dim  ; ++d)
            {
                m_bVStep(n + d*nodesCount) = 0;
                for(Eigen::SparseMatrix<double>::InnerIterator it(m_M, n + d*nodesCount); it; ++it)
                {
                    Eigen::Index row = it.row();
                    if(row == it.col())
                    {
                        it.valueRef() = 1;
                        continue;
                    }

                    it.valueRef() = 0;
                }
            }
        }
    }

    m_M.makeCompressed();
}

template<unsigned short dim>
void MomContEqIncompNewton<dim>::m_setupPicardFracStep(unsigned int maxIter, double minRes)
{
    m_pPicardAlgo = std::make_unique<PicardAlgo>([&](const auto& qPrevVec){
        m_clock.start();
        m_M.resize(qPrevVec[0].rows(), qPrevVec[0].rows());
        m_MK_dt.resize(qPrevVec[0].rows(), qPrevVec[0].rows());
        m_L.resize(qPrevVec[1].rows(), qPrevVec[1].rows());
        m_bVAppStep.resize(qPrevVec[0].rows()); m_bVAppStep.setZero();
        m_bPcorrStep.resize(qPrevVec[1].rows()); m_bPcorrStep.setZero();
        m_bVStep.resize(qPrevVec[0].rows()); m_bVStep.setZero();
        m_accumalatedTimes["Prepare Picard algorithm"] += m_clock.end();

        m_clock.start();
        m_pMesh->saveNodesList();
        m_accumalatedTimes["Save/restore nodelist"] += m_clock.end();
    },
    [&](auto& qIterVec, const auto& qPrevVec){
        m_buildMatFracStep(qPrevVec);
        m_clock.start();
        m_applyBCVAppStep(qPrevVec[0]);
        m_accumalatedTimes["Apply boundary conditions v app step"] += m_clock.end();
        m_clock.start();
        m_solverIt.compute(m_MK_dt);
        m_accumalatedTimes["Compute matrix v app step"] += m_clock.end();
        Eigen::VectorXd qVTilde(qPrevVec[0].rows());

        if(m_solverIt.info() == Eigen::Success)
        {
            m_clock.start();
            qVTilde = m_solverIt.solve(m_bVAppStep);
            m_accumalatedTimes["Solve system v app step"] += m_clock.end();
        }
        else
        {
            if(m_pProblem->isOutputVerbose())
                std::cout << "\t * The Eigen solver failed to factorize the matrix for the velocity guess step!" << std::endl;
            m_clock.start();
            m_pMesh->restoreNodesList();
            m_accumalatedTimes["Save/restore nodelist"] += m_clock.end();
            return false;
        }

        m_buildMatPcorrStep(qVTilde, qPrevVec[1]);
        m_clock.start();
        m_applyBCPCorrStep();
        m_accumalatedTimes["Apply boundary conditions p corr step"] += m_clock.end();
        m_clock.start();
        m_solverIt.compute(m_L);
        m_accumalatedTimes["Compute matrix p corr step"] += m_clock.end();
        Eigen::VectorXd qDeltaP(qPrevVec[1].rows());

        if(m_solverIt.info() == Eigen::Success)
        {
            m_clock.start();
            qIterVec[1] = m_solverIt.solve(m_bPcorrStep);
            m_accumalatedTimes["Solve system p corr step"] += m_clock.end();
        }
        else
        {
            if(m_pProblem->isOutputVerbose())
                std::cout << "\t * The Eigen solver failed to factorize the matrix for the pressure correction step!" << std::endl;
            m_clock.start();
            m_pMesh->restoreNodesList();
            m_accumalatedTimes["Save/restore nodelist"] += m_clock.end();
            return false;
        }

        m_clock.start();
        qDeltaP = qIterVec[1] - m_gammaFS*qPrevVec[1];
        m_accumalatedTimes["Update solutions"] += m_clock.end();

        m_buildMatVStep(qDeltaP);
        m_clock.start();
        m_applyBCVStep();
        m_accumalatedTimes["Apply boundary conditions v corr step"] += m_clock.end();
        m_clock.start();
        m_solverIt.compute(m_M);
        m_accumalatedTimes["Compute matrix v corr step"] += m_clock.end();
        Eigen::VectorXd qDeltaV(qPrevVec[0].rows());

        if(m_solverIt.info() == Eigen::Success)
        {
            m_clock.start();
            qDeltaV = m_solverIt.solve(m_bVStep);
            m_accumalatedTimes["Solve system v corr step"] += m_clock.end();
        }
        else
        {
            if(m_pProblem->isOutputVerbose())
                std::cout << "\t * The Eigen solver failed to factorize the matrix for the pressure correction step!" << std::endl;
            m_clock.start();
            m_pMesh->restoreNodesList();
            m_accumalatedTimes["Save/restore nodelist"] += m_clock.end();
            return false;
        }

        m_clock.start();
        qIterVec[0] = qVTilde + qDeltaV;

        setNodesStatesfromQ(m_pMesh, qIterVec[0], m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim() - 1);
        setNodesStatesfromQ(m_pMesh, qIterVec[1], m_pMesh->getDim(), m_pMesh->getDim());
        Eigen::VectorXd deltaPos = qIterVec[0]*m_pSolver->getTimeStep();
        m_pMesh->updateNodesPositionFromSave(deltaPos);
        m_accumalatedTimes["Update solutions"] += m_clock.end();

        return true;
     },
    [&](const auto& qIterVec, const auto& qIterPrevVec) -> double {
        m_clock.start();
        Mesh* p_Mesh = this->m_pMesh;
        const std::size_t nNodes = p_Mesh->getNodesCount();
        double resV = 0, resP = 0;

        double num = 0, den = 0;
        for(std::size_t n = 0 ; n < nNodes ; ++n)
        {
            const Node& node = p_Mesh->getNode(n);

            if(!node.isFree())
            {
                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    num += (qIterVec[0](n + d*nNodes) - qIterPrevVec[0](n + d*nNodes))*(qIterVec[0](n + d*nNodes) - qIterPrevVec[0](n + d*nNodes));
                    den += qIterPrevVec[0](n + d*nNodes)*qIterPrevVec[0](n + d*nNodes);
                }
            }
        }

        if(den == 0)
            resV = std::numeric_limits<double>::max();
        else
            resV = std::sqrt(num/den);

        if(m_computePres)
        {
            num = 0, den = 0;
            for(std::size_t n = 0 ; n < nNodes ; ++n)
            {
                const Node& node = p_Mesh->getNode(n);

                if(!node.isFree())
                {
                    num += (qIterVec[1](n) - qIterPrevVec[1](n))*(qIterVec[1](n) - qIterPrevVec[1](n));
                    den += qIterPrevVec[1](n)*qIterPrevVec[1](n);
                }
            }

            if(den == 0)
                resP = std::numeric_limits<double>::max();
            else
                resP = std::sqrt(num/den);
        }

        std::cout << "\t * Residual: v -> " << resV << ", p -> " << resP << std::endl;

        m_accumalatedTimes["Compute Picard Algo residual"] += m_clock.end();

        return std::max(resV, resP);
    }, maxIter, minRes);
}
