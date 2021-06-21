#include "PSmoother.hpp"

#include <set>
#include <vector>

template<unsigned short dim>
PSmoother<dim>::PSmoother(Problem* pProblem, Mesh& mesh, double a, double epsADRtoll, double betaInit):
MeshSmoother(pProblem, mesh),
m_a(a),
m_epsADRtoll(epsADRtoll),
m_beta(betaInit)
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

    m_pMatBuilder = std::make_unique<MatrixBuilder<dim>>(m_mesh, nGPHD, nGPLD);

    m_pMatBuilder->setKcomputeFactor([&](const Element& /** element **/,
                                         const NmatTypeHD<dim>& /** N **/,
                                         const BmatType<dim>& /** B **/,
                                         const DdevMatType<dim>& /** ddev **/) -> double {
        return 1; //G = 1
    });

    m_pMatBuilder->setDcomputeFactor([&](const Element& /** element **/,
                                         const NmatTypeHD<dim>& /** N **/,
                                         const BmatType<dim>& /** B **/) -> double {
        return 1; //K = 1
    });

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
}

template<unsigned short dim>
PSmoother<dim>::~PSmoother()
{

}

template<unsigned short dim>
void PSmoother<dim>::smooth(bool verboseOuput)
{
    bool shouldSmooth = m_init(verboseOuput);

    if(!shouldSmooth)
        return;

    m_mesh.saveNodesList();

    double res = std::numeric_limits<double>::max();
    double initBeta = m_beta;

    unsigned int step = 0;

    while(res > m_epsADRtoll)
    {
        Eigen::VectorXd vn_32 = m_vn_12;
        m_vn_12 = m_vn12;
        Eigen::VectorXd un_1 = m_un1;

        double maxU = 0;
        for(auto i = 0 ; i < m_un1.rows() ; ++i)
        {
            if(std::abs(m_un1[i]) > std::abs(maxU))
                maxU = m_un1[i];
        }
        std::cout << "Max u: " << maxU << std::endl;

        double maxV = 0;
        for(auto i = 0 ; i < m_un1.rows() ; ++i)
        {
            if(std::abs(m_vn12[i]) > std::abs(maxV))
                maxV = m_vn12[i];
        }
        std::cout << "Max v: " << maxV << std::endl;

        m_buildSystem(m_vn_12, m_un1);
        m_applyBC();

        m_vn12 = m_lhs*m_rhs;
        m_un1 = un_1 + m_dt*m_vn12;
        m_mesh.updateNodesPositionFromSave(m_un1, m_nodeProbToMesh);

        bool shouldGoOn = true;
        for(auto i = 0 ; i < m_un1.rows() ; ++i)
        {
            if(m_un1[i] == un_1[i])
                continue;

            double epsStab = (m_dt/4)*std::abs(m_vn12[i] - 2*m_vn_12[i] + vn_32[i])/std::abs(m_un1[i] - un_1[i]);
            if(epsStab > 1)
            {
                m_mesh.updateNodesPositionFromSave(un_1, m_nodeProbToMesh);
                m_un1 = un_1;
                m_vn12 = m_vn_12;
                m_vn_12 = vn_32;
                for(std::size_t j = 0 ; j < m_Fint.size() ; ++j)
                    m_Fint[j] = m_FintPrec[j];
                m_beta *= 2*epsStab;
                std::cout << "Unstability detected, increasing beta: " << epsStab << std::endl;
                shouldGoOn = false;
                break;
            }
        }

        if(!shouldGoOn)
            continue;

        char c;
        std::cin >> c;

        if(c == 'p')
            break;

        if(m_vn_12.squaredNorm() == 0)
            res = std::numeric_limits<double>::max();
        else
            res = (m_vn12 - m_vn_12).squaredNorm()/m_vn_12.squaredNorm();

        std::cout << "Step " << step << ": " << res << " vs " << m_epsADRtoll << std::endl;

        step++;
    }

    m_beta = initBeta;
}

template<unsigned short dim>
bool PSmoother<dim>::m_init(bool verboseOuput)
{
    const std::size_t nNodes = m_mesh.getNodesCount();
    const std::size_t nElm = m_mesh.getElementsCount();
    const double hchar = m_mesh.getHchar();

    constexpr double fact = (dim == 2) ? std::sqrt(12) : std::sqrt(24);


    //Keep only deformed elements and their neighbours
    m_elementsToKeep.resize(nElm, false);
    m_nodesMeshToBC.resize(nNodes, -1);
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);
        double rin = element.getRin();

        if(rin > m_a*hchar/fact)
            continue;

        //The element is too flat, will be smoothed with his neighbours
        m_elementsToKeep[elm] = true;

        //Only nodes from the deformed element that are not bound or on the FS will be allowed to move
        for(std::size_t n = 0 ; n < dim + 1 ; ++n)
        {
            const Node& node = element.getNode(n);
            std::size_t nodeIndexInMesh = element.getNodeIndex(n);
            if(!node.isOnFreeSurface() && !node.isBound())
                m_nodesMeshToBC[nodeIndexInMesh] = 0;
            else
                m_nodesMeshToBC[nodeIndexInMesh] = 1;
        }

        unsigned int nNeighbourElm = element.getNeighbourElementsCount();
        for(unsigned int  i = 0 ; i < nNeighbourElm ; ++i)
        {
            m_elementsToKeep[element.getNeighbourElmIndex(i)] = true;
            const Element& neigbourElement = element.getNeighbourElement(i);
            for(unsigned int  n = 0 ; n < dim + 1 ; ++n)
            {
                std::size_t nodeIndexInMesh = neigbourElement.getNodeIndex(n);
                if(m_nodesMeshToBC[nodeIndexInMesh] != 0)
                    m_nodesMeshToBC[nodeIndexInMesh] = 1;
            }
        }
    }

    //Associated with each kept element his index in the smoothing
    std::size_t nElmInProb = 0;
    m_elmMeshToProb.resize(nElm);
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        if(!m_elementsToKeep[elm])
            continue;

        m_elmMeshToProb[elm] = nElmInProb;
        nElmInProb++;
    }

    if(nElmInProb == 0)
        return false;

    //Count the number of nodes in the smoothing
    std::size_t nNodesInProb = 0;
    std::size_t nNodesInProbFixed = 0;
    for(std::size_t n = 0 ; n < nNodes ; ++n)
    {
        if(m_nodesMeshToBC[n] == -1)
            continue;

        if(m_nodesMeshToBC[n] == 1)
            nNodesInProbFixed++;

        nNodesInProb++;
    }

    m_nodeProbToMesh.resize(nNodesInProb);

    //Associate mesh index and smoothing index for the nodes
    std::size_t counter = 0;
    m_nodeMeshToProb.resize(nNodes);
    for(std::size_t n = 0 ; n < nNodes ; ++n)
    {
        if(m_nodesMeshToBC[n] == -1)
            continue;

        m_nodeMeshToProb[n] = counter;
        m_nodeProbToMesh[counter] = n;
        counter++;
    }

    m_un1.resize(nNodesInProb*dim); m_un1.setZero();
    m_vn12.resize(nNodesInProb*dim); m_vn12.setZero();
    m_vn_12.resize(nNodesInProb*dim); m_vn_12.setZero();
    m_Fint.resize(nElmInProb);
    m_FintPrec.resize(nElmInProb);
    for(auto& F : m_FintPrec)
        F.setZero();
    m_lhs.resize(nNodesInProb*dim);
    m_rhs.resize(nNodesInProb*dim);

    if(verboseOuput)
    {
        std::cout << "# elements to smooth: " << nElmInProb << " vs " << nElm << std::endl;
        std::cout << "# nodes in problem : " << nNodesInProb << " vs " << nNodes << std::endl;
        std::cout << "# nodes with BC : " << nNodesInProbFixed << std::endl;
    }

    if(nNodesInProb == nNodesInProbFixed)
        return false;

    return true;
}

template<unsigned short dim>
void PSmoother<dim>::m_buildSystem(const Eigen::VectorXd& vBar, const Eigen::VectorXd& uBar)
{
    const std::size_t nElm = m_mesh.getElementsCount();
    const std::size_t nNodesInProblem = m_nodeProbToMesh.size();
    const std::size_t nElmInProblem = m_Fint.size();
    const double hchar = m_mesh.getHchar();

    constexpr unsigned short nodPerEl = dim + 1;

    constexpr double fact = (dim == 2) ? std::sqrt(12) : std::sqrt(24);

    m_lhs.setZero();
    std::vector<Eigen::DiagonalMatrix<double, dim*nodPerEl>> lhse(nElmInProblem);
    m_rhs.setZero();
    std::vector<Eigen::Matrix<double, dim*nodPerEl, 1>> rhse(nElmInProblem);

    auto getElementVecStateSpec = [&](const Element& element, const Eigen::VectorXd& q) -> Eigen::Matrix<double, dim*nodPerEl, 1> {
        Eigen::Matrix<double, dim*nodPerEl, 1> vec;
        for(unsigned short n = 0 ; n < nodPerEl ; ++n)
        {
            std::size_t nodeIndex = element.getNodeIndex(n);
            if constexpr (dim == 2)
            {
                vec(n) = q[m_nodeMeshToProb[nodeIndex]];
                vec(n + nodPerEl) = q[m_nodeMeshToProb[nodeIndex] + nNodesInProblem];
            }
            else if constexpr (dim == 3)
            {
                vec(n) = q(m_nodeMeshToProb[nodeIndex]);
                vec(n + nodPerEl) = q(m_nodeMeshToProb[nodeIndex] + nNodesInProblem);
                vec(n + 2*nodPerEl) = q(m_nodeMeshToProb[nodeIndex] + 2*nNodesInProblem);
            }
        }

        return vec;
    };

    for(std::size_t i = 0 ; i < nElmInProblem ; ++i)
    {
        m_FintPrec[i] = m_Fint[i];
    }

    #pragma omp parallel for default(shared) schedule(dynamic)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        if(!m_elementsToKeep[elm])
            continue;


        const Element& element = m_mesh.getElement(elm);
        double rin = element.getRin();
        if(rin < 0)
            std::cerr << "Negative rin !" << std::endl;

        GradNmatType<dim> gradNe = m_pMatBuilder->getGradN(element);
        BmatType<dim> Be = m_pMatBuilder->getB(gradNe);
        Eigen::Matrix<double, dim*nodPerEl, dim*nodPerEl> Ke = m_pMatBuilder->getK(element, Be);
        Eigen::Matrix<double, nodPerEl, dim*nodPerEl> De = m_pMatBuilder->getD(element, Be);
        Eigen::Matrix<double, nodPerEl, 1> p; p.setOnes(); p *= (hchar/(rin*fact));
        Eigen::Matrix<double, dim*nodPerEl, 1> u = getElementVecStateSpec(element, uBar);
        Eigen::Matrix<double, dim*nodPerEl, 1> v = getElementVecStateSpec(element, vBar);
        m_Fint[m_elmMeshToProb[elm]] = Ke*u - De.transpose()*p;
        Eigen::DiagonalMatrix<double, dim*nodPerEl> S;
        auto& Sdiag = S.diagonal();
        for(unsigned short i = 0 ; i < S.rows() ; ++i)
        {
            if(std::abs(v[i]) > 1e-15)
                Sdiag[i] = (m_Fint[m_elmMeshToProb[elm]][i] - m_FintPrec[m_elmMeshToProb[elm]][i])/(m_dt*v[i]);
            else
                Sdiag[i] = 1/hchar;//(m_Fint[m_elmMeshToProb[elm]][i] - m_FintPrec[m_elmMeshToProb[elm]][i])/(hchar);
        }

        Eigen::DiagonalMatrix<double, dim*nodPerEl> M = (m_beta*m_dt*m_dt/4)*S;
        Eigen::DiagonalMatrix<double, dim*nodPerEl> C = std::sqrt(m_beta)*m_dt*S;

        lhse[m_elmMeshToProb[elm]].diagonal() = 2*M.diagonal() + m_dt*C.diagonal();
        rhse[m_elmMeshToProb[elm]] = 2*M*v - m_dt*C*v - 2*m_dt*m_Fint[m_elmMeshToProb[elm]];
    }

    auto& lhsdiag = m_lhs.diagonal();

    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);

        if(!m_elementsToKeep[elm])
            continue;

        for(unsigned short i = 0 ; i < dim + 1 ; ++i)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                /********************************************************************
                                             Build M
                ********************************************************************/
                lhsdiag[m_nodeMeshToProb[element.getNodeIndex(i)] + d*nNodesInProblem] += lhse[m_elmMeshToProb[elm]].diagonal()[i + d*nodPerEl];

                /************************************************************************
                                                Build f
                ************************************************************************/
                m_rhs(m_nodeMeshToProb[element.getNodeIndex(i)] + d*nNodesInProblem) += rhse[m_elmMeshToProb[elm]](i + d*nodPerEl);
            }
        }
    }

    MatrixBuilder<dim>::inverse(m_lhs);
}

template<unsigned short dim>
void PSmoother<dim>::m_applyBC()
{
    auto& lhsdiag = m_lhs.diagonal();

    const std::size_t nNodesInProb = m_nodeProbToMesh.size();

    for (std::size_t n = 0 ; n < nNodesInProb ; ++n)
    {
        if(m_nodesMeshToBC[m_nodeProbToMesh[n]] == 1)
        {
            for(unsigned short d = 0 ; d < dim ; ++d)
            {
                m_rhs[n + d*nNodesInProb] = 0;
                lhsdiag[n + d*nNodesInProb] = 1;
            }
        }
    }
}
