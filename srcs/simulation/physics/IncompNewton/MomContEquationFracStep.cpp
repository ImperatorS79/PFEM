#include "MomContEquation.hpp"
#include "../../Problem.hpp"
#include "../../Solver.hpp"
#include "../../utility/StatesFromToQ.hpp"


void MomContEqIncompNewton::m_buildMatFracStep(Eigen::SparseMatrix<double>& M,
                                               Eigen::SparseMatrix<double>& MK_dt,
                                               Eigen::SparseMatrix<double>& dtL,
                                               std::vector<Eigen::MatrixXd>& DTelm,
                                               std::vector<Eigen::MatrixXd>& Lelm, Eigen::VectorXd& bVAppStep,
                                               const std::vector<Eigen::VectorXd>& qPrev)
{
    unsigned int dim = m_pMesh->getDim();
    unsigned int noPerEl = m_pMesh->getNodesPerElm();
    unsigned int tripletPerElmM = (dim*noPerEl*noPerEl);
    unsigned int tripletPerElmMKDT = (dim*noPerEl*noPerEl + dim*noPerEl*dim*noPerEl);
    unsigned int tripletPerElmDtL = (noPerEl*noPerEl);
    unsigned int doubletPerElm = (3*dim*noPerEl);

    std::size_t nElm = m_pMesh->getElementsCount();
    std::size_t nNodes = m_pMesh->getNodesCount();
    double dt = m_pSolver->getTimeStep();

    std::vector<Eigen::Triplet<double>> indexM(tripletPerElmM*nElm);
    std::vector<Eigen::Triplet<double>> indexMK_dt(tripletPerElmMKDT*nElm);
    std::vector<Eigen::Triplet<double>> indexDtL(tripletPerElmDtL*nElm);
    std::vector<std::pair<std::size_t, double>> indexbVappStep(doubletPerElm*nElm); bVAppStep.setZero();
    DTelm.resize(nElm);
    Lelm.resize(nElm);
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        DTelm[elm].resize(dim*noPerEl, noPerEl);
        Lelm[elm].resize(noPerEl, noPerEl);
    }

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::MatrixXd gradNe = m_pMatBuilder->getGradN(element);
        Eigen::MatrixXd Be = m_pMatBuilder->getB(gradNe);
        Eigen::MatrixXd Me = m_pMatBuilder->getM(element);
        Eigen::MatrixXd Me_dt = (1/dt)*Me;
        Me = MatrixBuilder::diagBlock(Me, dim);
        Me_dt = MatrixBuilder::diagBlock(Me_dt, dim);
        Eigen::MatrixXd Ke = m_pMatBuilder->getK(element, Be);
        DTelm[elm] = m_pMatBuilder->getD(element, Be).transpose();
        Lelm[elm] = m_pMatBuilder->getL(element, Be, gradNe);
        Eigen::VectorXd Fe = m_pMatBuilder->getF(element, m_bodyForce, Be);

        Eigen::VectorXd vPrev(noPerEl*dim);
        Eigen::VectorXd pPrev(noPerEl);
        if(dim == 2)
        {
            vPrev << qPrev[0][element.getNodeIndex(0)], qPrev[0][element.getNodeIndex(1)], qPrev[0][element.getNodeIndex(2)],
                     qPrev[0][element.getNodeIndex(0) + nNodes], qPrev[0][element.getNodeIndex(1) + nNodes], qPrev[0][element.getNodeIndex(2) + nNodes];

            pPrev << qPrev[1][element.getNodeIndex(0)], qPrev[1][element.getNodeIndex(1)], qPrev[1][element.getNodeIndex(2)];
        }
        else
        {
            vPrev << qPrev[0][element.getNodeIndex(0)], qPrev[0][element.getNodeIndex(1)], qPrev[0][element.getNodeIndex(2)], qPrev[0][element.getNodeIndex(3)],
                     qPrev[0][element.getNodeIndex(0) + nNodes], qPrev[0][element.getNodeIndex(1) + nNodes], qPrev[0][element.getNodeIndex(2) + nNodes], qPrev[0][element.getNodeIndex(3) + nNodes],
                     qPrev[0][element.getNodeIndex(0) + 2*nNodes], qPrev[0][element.getNodeIndex(1) + 2*nNodes], qPrev[0][element.getNodeIndex(2) + 2*nNodes], qPrev[0][element.getNodeIndex(3) + 2*nNodes];

            pPrev << qPrev[1][element.getNodeIndex(0)], qPrev[1][element.getNodeIndex(1)], qPrev[1][element.getNodeIndex(2)], qPrev[1][element.getNodeIndex(3)];
        }

        Eigen::VectorXd MvPreve_dt = Me_dt*vPrev;
        Eigen::VectorXd DTpPreve = m_gammaFS*DTelm[elm]*pPrev;

        std::size_t countM = 0;
        std::size_t countMK_dt = 0;
        std::size_t countDtL = 0;
        std::size_t countbVappStep = 0;

        for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            const Node& ni = m_pMesh->getNode(element.getNodeIndex(i));

            for(unsigned short j = 0 ; j < noPerEl ; ++j)
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
                                                   Me(i + d*noPerEl, j + d*noPerEl));

                        indexMK_dt[tripletPerElmMKDT*elm + countMK_dt] =
                            Eigen::Triplet<double>(element.getNodeIndex(i) + d*nNodes,
                                                   element.getNodeIndex(j) + d*nNodes,
                                                   Me_dt(i + d*noPerEl, j + d*noPerEl));
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
                                                       Ke(i + d*noPerEl, j + d2*noPerEl));
                        }
                        countMK_dt++;
                    }
                }

                /********************************************************************
                                            Build L
                ********************************************************************/
                if(!ni.isFree() && !ni.isOnFreeSurface())
                {
                    indexDtL[tripletPerElmDtL*elm + countDtL] =
                        Eigen::Triplet<double>(element.getNodeIndex(i),
                                               element.getNodeIndex(j),
                                               Lelm[elm](i,j));
                }
                countDtL++;
            }

            /************************************************************************
                                              Build f
            ************************************************************************/
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                if(!(ni.isBound() || ni.isFree()))
                {
                    indexbVappStep[doubletPerElm*elm + countbVappStep] = std::make_pair(element.getNodeIndex(i) + d*nNodes, Fe(i + d*noPerEl));
                    countbVappStep++;

                    indexbVappStep[doubletPerElm*elm + countbVappStep] = std::make_pair(element.getNodeIndex(i) + d*nNodes, MvPreve_dt(i + d*noPerEl));
                    countbVappStep++;

                    indexbVappStep[doubletPerElm*elm + countbVappStep] = std::make_pair(element.getNodeIndex(i) + d*nNodes, DTpPreve(i + d*noPerEl));
                    countbVappStep++;
                }
            }
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());;
//    char c;
//    std::cin >> c;

    //Best would be to know the number of nodes in which case :/
    //This can still be fasten using OpenMP but will never be as good as using []
    //with preallocated memory
    for(std::size_t n = 0 ; n < nNodes ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);

        if(node.isFree() || node.isOnFreeSurface())
        {
            indexDtL.push_back(Eigen::Triplet<double>(n, n, 1));
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


    /********************************************************************************
                                        Compute A and b
    ********************************************************************************/
    M.setFromTriplets(indexM.begin(), indexM.end());
    MK_dt.setFromTriplets(indexMK_dt.begin(), indexMK_dt.end());
    dtL.setFromTriplets(indexDtL.begin(), indexDtL.end());

    for(const auto& doublet : indexbVappStep)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        bVAppStep[doublet.first] += doublet.second;
    }

}

void MomContEqIncompNewton::m_applyBCVAppStep(Eigen::SparseMatrix<double>& MK_dt, Eigen::VectorXd& bVAppStep, const Eigen::VectorXd& qPrev)
{
    const uint8_t dim = m_pMesh->getDim();
    const std::size_t nodesCount = m_pMesh->getNodesCount();

    const std::size_t facetsCount = m_pMesh->getFacetsCount();
    const std::size_t noPerFacet = m_pMesh->getNodesPerFacet();

    for(std::size_t f = 0 ; f < facetsCount ; ++f)
    {
        if(m_gamma < 1e-15)
            continue;

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

        Eigen::MatrixXd MGamma = m_pMatBuilder->getMGamma(facet);
        MGamma = MatrixBuilder::diagBlock(MGamma, m_pMesh->getDim());

        Eigen::VectorXd kappa_n(dim*noPerFacet);
        for(uint8_t n = 0 ; n < noPerFacet; ++n)
        {
            double curvature = m_pMesh->getFreeSurfaceCurvature(facet.getNodeIndex(n));
            std::array<double, 3> normal = m_pMesh->getBoundFSNormal(facet.getNodeIndex(n));

            for(uint8_t d = 0 ; d < dim ; ++d)
                kappa_n(n + d*noPerFacet) = curvature*normal[d];
        }

        Eigen::VectorXd Ff = MGamma*kappa_n;

        for(unsigned short i = 0 ; i < noPerFacet ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
                bVAppStep(facet.getNodeIndex(i) + d*nodesCount) += Ff[d*noPerFacet + i];
        }
    }

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree())
        {
            if(!node.isBound())
            {
                for(uint8_t d = 0 ; d < dim  ; ++d)
                {
                    bVAppStep(n + d*nodesCount) = qPrev(n + d*nodesCount) + m_pSolver->getTimeStep()*m_bodyForce[d];
                }
            }
        }

        if(node.isBound())
        {
            if(node.getFlag(m_bcFlags[0]))
            {
                std::vector<double> result;
                result = m_bcParams[0].call<std::vector<double>>(m_pMesh->getNodeType(n) + "V",
                                                        node.getPosition(),
                                                        m_pMesh->getBoundNodeInitPos(n),
                                                        m_pProblem->getCurrentSimTime() +
                                                        m_pSolver->getTimeStep());

                for(uint8_t d = 0 ; d < dim ; ++d)
                {
                    bVAppStep(n + d*nodesCount) = result[d];
                    for(Eigen::SparseMatrix<double>::InnerIterator it(MK_dt, n + d*nodesCount); it; ++it)
                    {
                        Eigen::Index row = it.row();
                        if(row == it.col())
                            continue;

                        double value = it.value();
                        bVAppStep(row) -= value*result[d];
                        it.valueRef() = 0;
                    }
                }
            }
        }
    }

    MK_dt.makeCompressed();
}

void MomContEqIncompNewton::m_buildMatPcorrStep(Eigen::VectorXd& bPcorrStep, const Eigen::VectorXd& qVTilde, const Eigen::VectorXd& qPprev)
{
    unsigned int dim = m_pMesh->getDim();
    unsigned int noPerEl = m_pMesh->getNodesPerElm();
    unsigned int doubletPerElm = 2*noPerEl;

    std::size_t nElm = m_pMesh->getElementsCount();
    std::size_t nNodes = m_pMesh->getNodesCount();

    double dt = m_pSolver->getTimeStep();

    std::vector<std::pair<std::size_t, double>> indexbPcorrStep(doubletPerElm*nElm); bPcorrStep.setZero();

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::VectorXd vTilde(noPerEl*dim);
        Eigen::VectorXd pPrev(noPerEl);
        if(dim == 2)
        {
            vTilde << qVTilde[element.getNodeIndex(0)], qVTilde[element.getNodeIndex(1)], qVTilde[element.getNodeIndex(2)],
                      qVTilde[element.getNodeIndex(0) + nNodes], qVTilde[element.getNodeIndex(1) + nNodes], qVTilde[element.getNodeIndex(2) + nNodes];
            pPrev << qPprev[element.getNodeIndex(0)], qPprev[element.getNodeIndex(1)], qPprev[element.getNodeIndex(2)];
        }
        else
        {
            vTilde << qVTilde[element.getNodeIndex(0)], qVTilde[element.getNodeIndex(1)], qVTilde[element.getNodeIndex(2)], qVTilde[element.getNodeIndex(3)],
                      qVTilde[element.getNodeIndex(0) + nNodes], qVTilde[element.getNodeIndex(1) + nNodes], qVTilde[element.getNodeIndex(2) + nNodes], qVTilde[element.getNodeIndex(3) + nNodes],
                      qVTilde[element.getNodeIndex(0) + 2*nNodes], qVTilde[element.getNodeIndex(1) + 2*nNodes], qVTilde[element.getNodeIndex(2) + 2*nNodes], qVTilde[element.getNodeIndex(3) + 2*nNodes];
            pPrev << qPprev[element.getNodeIndex(0)], qPprev[element.getNodeIndex(1)], qPprev[element.getNodeIndex(2)], qPprev[element.getNodeIndex(3)];
        }

        Eigen::VectorXd rho_dt_DvTilde = (m_rho/dt)*m_DTelm[elm].transpose()*vTilde;
        Eigen::VectorXd LpPrev = m_gammaFS*m_Lelm[elm]*pPrev;

        std::size_t countbPcorrStep = 0;
        for(unsigned short i = 0 ; i < noPerEl ; ++i)
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

    for(const auto& doublet : indexbPcorrStep)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        bPcorrStep[doublet.first] += doublet.second;
    }
}

void MomContEqIncompNewton::m_applyBCPCorrStep(Eigen::SparseMatrix<double>& L, Eigen::VectorXd& bPcorrStep)
{
    std::size_t nodesCount = m_pMesh->getNodesCount();

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree() || node.isOnFreeSurface())
        {
            bPcorrStep[n] = 0;

            for(Eigen::SparseMatrix<double>::InnerIterator it(L, n); it; ++it)
            {
                Eigen::Index row = it.row();
                if(row == it.col())
                    continue;

                it.valueRef() = 0;
            }
        }
    }

    L.makeCompressed();
}

void MomContEqIncompNewton::m_buildMatVStep(Eigen::VectorXd& bVStep, const Eigen::VectorXd& qDeltaP)
{
    unsigned int dim = m_pMesh->getDim();
    unsigned int noPerEl = m_pMesh->getNodesPerElm();
    unsigned int doubletPerElm = dim*noPerEl;
    double dt = m_pSolver->getTimeStep();

    std::size_t nElm = m_pMesh->getElementsCount();
    std::size_t nNodes = m_pMesh->getNodesCount();

    std::vector<std::pair<std::size_t, double>> indexbVStep(doubletPerElm*nElm); bVStep.setZero();

    Eigen::setNbThreads(1);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_pMesh->getElement(elm);

        Eigen::VectorXd deltaP(noPerEl);
        if(dim == 2)
            deltaP << qDeltaP[element.getNodeIndex(0)], qDeltaP[element.getNodeIndex(1)], qDeltaP[element.getNodeIndex(2)];
        else
            deltaP << qDeltaP[element.getNodeIndex(0)], qDeltaP[element.getNodeIndex(1)], qDeltaP[element.getNodeIndex(2)], qDeltaP[element.getNodeIndex(3)];

        Eigen::VectorXd DTdeltaP = m_DTelm[elm]*deltaP;

        std::size_t countbVStep = 0;
        for(unsigned short i = 0 ; i < noPerEl ; ++i)
        {
            const Node& node = element.getNode(i);

            if(!node.isFree() && !node.isBound())
            {
                /************************************************************************
                                                  Build f
                ************************************************************************/
                for(unsigned short d = 0 ; d < dim ; ++d)
                {
                    indexbVStep[doubletPerElm*elm + countbVStep] = std::make_pair(element.getNodeIndex(i) + d*nNodes, dt*DTdeltaP(i + d*noPerEl));
                    countbVStep++;
                }
            }
        }
    }
    Eigen::setNbThreads(m_pProblem->getThreadCount());

    for(const auto& doublet : indexbVStep)
    {
        //std::cout << doublet.first << ", " << doublet.second << std::endl;
        bVStep[doublet.first] += doublet.second;
    }
}

void MomContEqIncompNewton::m_applyBCVStep(Eigen::SparseMatrix<double>& M, Eigen::VectorXd& bVStep)
{
    std::size_t nodesCount = m_pMesh->getNodesCount();
    std::size_t dim = m_pMesh->getDim();

    //Do not parallelize this (lua)
    for (std::size_t n = 0 ; n < nodesCount ; ++n)
    {
        const Node& node = m_pMesh->getNode(n);
        if(node.isFree() || node.isBound())
        {
            for(uint8_t d = 0 ; d < dim  ; ++d)
            {
                bVStep(n + d*nodesCount) = 0;
                for(Eigen::SparseMatrix<double>::InnerIterator it(M, n + d*nodesCount); it; ++it)
                {
                    Eigen::Index row = it.row();
                    if(row == it.col())
                        continue;

                    it.valueRef() = 0;
                }
            }
        }
    }

    M.makeCompressed();
}

void MomContEqIncompNewton::m_setupPicardFracStep(unsigned int maxIter, double minRes)
{
    m_pPicardAlgo = std::make_unique<PicardAlgo>([&](const auto& qPrevVec){
        m_M.resize(qPrevVec[0].rows(), qPrevVec[0].rows());
        m_MK_dt.resize(qPrevVec[0].rows(), qPrevVec[0].rows());
        m_L.resize(qPrevVec[1].rows(), qPrevVec[1].rows());
        m_bVAppStep.resize(qPrevVec[0].rows()); m_bVAppStep.setZero();
        m_bPcorrStep.resize(qPrevVec[1].rows()); m_bPcorrStep.setZero();
        m_bVStep.resize(qPrevVec[0].rows()); m_bVStep.setZero();

        m_pMesh->saveNodesList();
    },
    [&](auto& qIterVec, const auto& qPrevVec){
        m_buildMatFracStep(m_M, m_MK_dt, m_L, m_DTelm, m_Lelm, m_bVAppStep, qPrevVec);
        m_applyBCVAppStep(m_MK_dt, m_bVAppStep, qPrevVec[0]);
        m_solverIt.compute(m_MK_dt);
        Eigen::VectorXd qVTilde(qPrevVec[0].rows());

        if(m_solverIt.info() == Eigen::Success)
            qVTilde = m_solverIt.solve(m_bVAppStep);
        else
        {
            if(m_pProblem->isOutputVerbose())
                std::cout << "\t * The Eigen solver failed to factorize the matrix for the velocity guess step!" << std::endl;
            m_pMesh->restoreNodesList();
            return false;
        }

        m_buildMatPcorrStep(m_bPcorrStep, qVTilde, qPrevVec[1]);
        m_applyBCPCorrStep(m_L, m_bPcorrStep);
        m_solverIt.compute(m_L);
        Eigen::VectorXd qDeltaP(qPrevVec[1].rows());

        if(m_solverIt.info() == Eigen::Success)
            qIterVec[1] = m_solverIt.solve(m_bPcorrStep);
        else
        {
            if(m_pProblem->isOutputVerbose())
                std::cout << "\t * The Eigen solver failed to factorize the matrix for the pressure correction step!" << std::endl;
            m_pMesh->restoreNodesList();
            return false;
        }

        qDeltaP = qIterVec[1] - m_gammaFS*qPrevVec[1];

        m_buildMatVStep(m_bVStep, qDeltaP);
        m_applyBCVStep(m_M, m_bVStep);
        m_solverIt.compute(m_M);
        Eigen::VectorXd qDeltaV(qPrevVec[0].rows());

        if(m_solverIt.info() == Eigen::Success)
            qDeltaV = m_solverIt.solve(m_bVStep);
        else
        {
            if(m_pProblem->isOutputVerbose())
                std::cout << "\t * The Eigen solver failed to factorize the matrix for the pressure correction step!" << std::endl;
            m_pMesh->restoreNodesList();
            return false;
        }

        qIterVec[0] = qVTilde + qDeltaV;

        setNodesStatesfromQ(m_pMesh, qIterVec[0], m_statesIndex[0], m_statesIndex[0] + m_pMesh->getDim() - 1);
        setNodesStatesfromQ(m_pMesh, qIterVec[1], m_pMesh->getDim(), m_pMesh->getDim());
        Eigen::VectorXd deltaPos = qIterVec[0]*m_pSolver->getTimeStep();
        m_pMesh->updateNodesPositionFromSave(std::vector<double> (deltaPos.data(), deltaPos.data() + m_pMesh->getDim()*m_pMesh->getNodesCount()));

        return true;
     },
    [&](const auto& qIterVec, const auto& qIterPrevVec) -> double {
        double num = 0, den = 0;
        Mesh* p_Mesh = this->m_pMesh;

        for(std::size_t n = 0 ; n < p_Mesh->getNodesCount() ; ++n)
        {
            const Node& node = p_Mesh->getNode(n);

            if(!node.isFree())
            {
                for(unsigned short d = 0 ; d < p_Mesh->getDim() ; ++d)
                {
                    num += (qIterVec[0](n + d*p_Mesh->getNodesCount()) - qIterPrevVec[0](n + d*p_Mesh->getNodesCount()))*(qIterVec[0](n + d*p_Mesh->getNodesCount()) - qIterPrevVec[0](n + d*p_Mesh->getNodesCount()));
                    den += qIterPrevVec[0](n + d*p_Mesh->getNodesCount())*qIterPrevVec[0](n + d*p_Mesh->getNodesCount());
                }
            }
        }

        if(den == 0)
            return std::numeric_limits<double>::max();

        return std::sqrt(num/den);
    }, maxIter, minRes);
}
