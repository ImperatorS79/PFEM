#include "MatricesBuilder.hpp"

#include "../../mesh/Mesh.hpp"

template<unsigned short dim, unsigned short noPerEl>
MatrixBuilder<dim, noPerEl>::MatrixBuilder(const Mesh& mesh, unsigned int nGPHD, unsigned int nGPLD):
m_mesh(mesh),
m_nGPHD(nGPHD),
m_nGPLD(nGPLD)
{
    static_assert(dim == 2 || dim == 3, "MatrixBuilder can only be used with dimension 2 or 3!");

    m_sfsHD = m_mesh.getShapeFunctions(dim, m_nGPHD);
    m_gradsfsHD = m_mesh.getGradShapeFunctions(dim);
    m_sfsLD = m_mesh.getShapeFunctions(dim - 1, m_nGPLD);
    m_gradsfsLD = m_mesh.getGradShapeFunctions(dim - 1);

    m_NHD.resize(m_nGPHD);
    m_NhdTNhd.resize(m_nGPHD);
    m_NLD.resize(m_nGPLD);
    m_NldTNld.resize(m_nGPLD);
    m_NHDtilde.resize(m_nGPHD);
    m_NLDtilde.resize(m_nGPLD);

    unsigned int counter = 0;
    for(std::vector<double> sf : m_sfsHD)
    {
        m_NHDtilde[counter].setZero();

        if constexpr (dim == 2)
        {
            m_NHDtilde[counter](0,0) = m_NHDtilde[counter](1,3) =  m_NHD[counter](0, 0) = sf[0];
            m_NHDtilde[counter](0,1) = m_NHDtilde[counter](1,4) =  m_NHD[counter](0, 1) = sf[1];
            m_NHDtilde[counter](0,2) = m_NHDtilde[counter](1,5) =  m_NHD[counter](0, 2) = sf[2];
        }
        else if constexpr (dim == 3)
        {
            m_NHDtilde[counter](0,0) = m_NHDtilde[counter](1,4) = m_NHDtilde[counter](2,8)  = m_NHD[counter](0, 0) = sf[0];
            m_NHDtilde[counter](0,1) = m_NHDtilde[counter](1,5) = m_NHDtilde[counter](2,9)  = m_NHD[counter](0, 1) = sf[1];
            m_NHDtilde[counter](0,2) = m_NHDtilde[counter](1,6) = m_NHDtilde[counter](2,10) = m_NHD[counter](0, 2) = sf[2];
            m_NHDtilde[counter](0,3) = m_NHDtilde[counter](1,7) = m_NHDtilde[counter](2,11) = m_NHD[counter](0, 3) = sf[3];
        }

        counter++;
    }

    counter = 0;
    for(std::vector<double> sf : m_sfsLD)
    {
        m_NLDtilde[counter].setZero();
        if constexpr (dim == 2)
        {
            m_NLDtilde[counter](0, 0) = m_NLDtilde[counter](1, 2) = m_NLD[counter](0, 0) = sf[0];
            m_NLDtilde[counter](0, 1) = m_NLDtilde[counter](1, 3) = m_NLD[counter](0, 1) = sf[1];
        }
        else if constexpr (dim == 3)
        {
            m_NLDtilde[counter](0, 0) = m_NLDtilde[counter](1, 3) = m_NLDtilde[counter](2, 6) = m_NLD[counter](0, 0) = sf[0];
            m_NLDtilde[counter](0, 1) = m_NLDtilde[counter](1, 4) = m_NLDtilde[counter](2, 7) = m_NLD[counter](0, 1) = sf[1];
            m_NLDtilde[counter](0, 2) = m_NLDtilde[counter](1, 5) = m_NLDtilde[counter](2, 8) = m_NLD[counter](0, 2) = sf[2];
        }

        counter++;
    }

    m_gaussWeightHD = m_mesh.getGaussWeight(dim, m_nGPHD);
    m_gaussPointsHD = m_mesh.getGaussPoints(dim, m_nGPHD);
    m_gaussPointsLD = m_mesh.getGaussPoints(dim - 1, m_nGPLD);
    m_gaussWeightLD = m_mesh.getGaussWeight(dim - 1, m_nGPLD);

    m_sum_NhdTNhd_w.setZero();
    for(std::size_t i = 0 ; i < nGPHD ; ++i)
    {
        m_NhdTNhd[i] = m_NHD[i].transpose()*m_NHD[i];
        m_sum_NhdTNhd_w += m_NhdTNhd[i]*m_gaussWeightHD[i];
    }

    m_sum_NldTNld_w.setZero();
    for(std::size_t i = 0 ; i < nGPLD ; ++i)
    {
        m_NldTNld[i] = m_NLD[i].transpose()*m_NLD[i];
        m_sum_NldTNld_w += m_NldTNld[i]*m_gaussWeightLD[i];
    }
}

template<unsigned short dim, unsigned short noPerEl>
MatrixBuilder<dim, noPerEl>::~MatrixBuilder()
{

}

template<unsigned short dim, unsigned short noPerEl>
auto MatrixBuilder<dim, noPerEl>::getGradN(const Element& element) -> MatrixBuilder<dim, noPerEl>::GradNmatType
{
    GradNmatType gradN;
    gradN.setZero();
    //I could use m_gradsfHD though...
    if constexpr (dim == 2)
    {
        gradN(0,0) = - element.getInvJ(0, 0) - element.getInvJ(1, 0);
        gradN(0,1) = element.getInvJ(0, 0);
        gradN(0,2) = element.getInvJ(1, 0);

        gradN(1,0) = - element.getInvJ(0, 1) - element.getInvJ(1, 1);
        gradN(1,1) = element.getInvJ(0, 1);
        gradN(1,2) = element.getInvJ(1, 1);
    }
    else if constexpr (dim == 3)
    {
        gradN(0,0) = - element.getInvJ(0, 0) - element.getInvJ(1, 0) - element.getInvJ(2, 0);
        gradN(0,1) = element.getInvJ(0, 0);
        gradN(0,2) = element.getInvJ(1, 0);
        gradN(0,3) = element.getInvJ(2, 0);

        gradN(1,0) = - element.getInvJ(0, 1) - element.getInvJ(1, 1) - element.getInvJ(2, 1);
        gradN(1,1) = element.getInvJ(0, 1);
        gradN(1,2) = element.getInvJ(1, 1);
        gradN(1,3) = element.getInvJ(2, 1);

        gradN(2,0) = - element.getInvJ(0, 2) - element.getInvJ(1, 2) - element.getInvJ(2, 2);
        gradN(2,1) = element.getInvJ(0, 2);
        gradN(2,2) = element.getInvJ(1, 2);
        gradN(2,3) = element.getInvJ(2, 2);
    }

    return gradN;
}

template<unsigned short dim, unsigned short noPerEl>
auto MatrixBuilder<dim, noPerEl>::getB(const GradNmatType& gradN) -> MatrixBuilder<dim, noPerEl>::BmatType
{
    BmatType B; B.setZero();

    //With linear shape functions, Be is constant for all gauss points ^^.
    if constexpr (dim == 2)
    {
        B(0,0) = B(2,3) = gradN(0,0);
        B(0,1) = B(2,4) = gradN(0,1);
        B(0,2) = B(2,5) = gradN(0,2);

        B(1,3) = B(2,0) = gradN(1,0);
        B(1,4) = B(2,1) = gradN(1,1);
        B(1,5) = B(2,2) = gradN(1,2);
    }
    else if constexpr (dim == 3)
    {
        B(0,0) = B(3,4) = B(4,8) = gradN(0,0);
        B(0,1) = B(3,5) = B(4,9) = gradN(0,1);
        B(0,2) = B(3,6) = B(4,10) = gradN(0,2);
        B(0,3) = B(3,7) = B(4,11) = gradN(0,3);

        B(1,4) = B(3,0) = B(5,8) = gradN(1,0);
        B(1,5) = B(3,1) = B(5,9) = gradN(1,1);
        B(1,6) = B(3,2) = B(5,10) = gradN(1,2);
        B(1,7) = B(3,3) = B(5,11) = gradN(1,3);

        B(2,8) = B(4,0) = B(5,4) = gradN(2,0);
        B(2,9) = B(4,1) = B(5,5) = gradN(2,1);
        B(2,10) = B(4,2) = B(5,6) = gradN(2,2);
        B(2,11) = B(4,3) = B(5,7) = gradN(2,3);
    }

    return B;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, dim*dim - 2*dim +3, 1> MatrixBuilder<dim, noPerEl>::getP(const Facet& facet)
{
    Eigen::Matrix<double, dim*dim - 2*dim +3, 1> P;
    std::array<double, 3> normal = facet.getNormal();

    if constexpr (dim == 2)
    {
        P[0] = 1 - normal[0]*normal[0];
        P[1] = 1 - normal[1]*normal[1];
        P[2] = - normal[0]*normal[1];
    }
    else if constexpr (dim == 3)
    {
        P[0] = 1 - normal[0]*normal[0];
        P[1] = 1 - normal[1]*normal[1];
        P[2] = 1 - normal[2]*normal[2];
        P[3] = - normal[0]*normal[1];
        P[4] = - normal[0]*normal[2];
        P[5] = - normal[1]*normal[2];
    }

    return P;
}


template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, dim*dim - 2*dim +3, dim*dim - 2*dim +3> MatrixBuilder<dim, noPerEl>::getT(const Eigen::Matrix<double, dim*dim - 2*dim +3, 1>& P)
{
    Eigen::Matrix<double, dim*dim - 2*dim +3, dim*dim - 2*dim +3> T;

    if constexpr (dim == 2)
    {
        T(0,0) = P[0]*P[0];         T(0,1) = P[2]*P[2];         T(0,2) = 2*P[0]*P[2];
        T(1,0) = T(0,1);            T(1,1) = P[1]*P[1];         T(1,2) = 2*P[2]*P[1];
        T(2,0) = P[0]*P[2];         T(2,1) = P[2]*P[1];         T(2,2) = P[0]*P[1] + P[2]*P[2];
    }
    else if constexpr (dim == 3)
    {
        T(0,0) = P[0]*P[0];         T(0,1) = P[3]*P[3];         T(0,2) = P[4]*P[4];     T(0,3) = 2*P[0]*P[3];           T(0,4) = 2*P[4]*P[3];               T(0,5) = 2*P[0]*P[4];
        T(1,0) = T(0,1);            T(1,1) = P[1]*P[1];         T(1,2) = P[5]*P[5];     T(1,3) = 2*P[3]*P[1];           T(1,4) = 2*P[5]*P[1];               T(1,5) = 2*P[3]*P[5];
        T(2,0) = T(0,2);            T(2,1) = T(1,2);            T(2,2) = P[2]*P[2];     T(2,3) = 2*P[5]*P[4];           T(2,4) = 2*P[5]*P[2];               T(2,5) = 2*P[2]*P[4];
        T(3,0) = P[0]*P[3];         T(3,1) = P[3]*P[1];         T(3,2) = P[4]*P[5];     T(3,3) = P[0]*P[1] + P[3]*P[3]; T(3,4) = P[4]*P[1] + P[5]*P[3];     T(3,5) = P[4]*P[3] + P[0]*P[5];
        T(4,0) = P[3]*P[4];         T(4,1) = P[5]*P[1];         T(4,2) = P[5]*P[2];     T(4,3) = T(3,4);                T(4,4) = P[1]*P[2] + P[5]*P[5];     T(4,5) = P[5]*P[4] + P[3]*P[2];
        T(5,0) = P[0]*P[4];         T(5,1) = P[3]*P[5];         T(5,2) = P[2]*P[4];     T(5,3) = T(3,5);                T(5,4) = T(4,5);                    T(5,5) = P[2]*P[0] + P[4]*P[4];
    }

    return T;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, noPerEl, noPerEl> MatrixBuilder<dim, noPerEl>::getM(const Element& element)
{
    Eigen::Matrix<double, noPerEl, noPerEl>  M; M.setZero();

    for(unsigned int i = 0 ; i < m_NhdTNhd.size() ; ++ i)
    {
        M += m_Mfunc(element, m_NHD[i])*m_NhdTNhd[i]*m_gaussWeightHD[i];
    }

    M *= element.getDetJ()*m_mesh.getRefElementSize(dim);

    return M;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, noPerEl - 1, noPerEl - 1> MatrixBuilder<dim, noPerEl>::getMGamma(const Facet& facet)
{
    Eigen::Matrix<double, noPerEl - 1, noPerEl - 1> MGamma; MGamma.setZero();

    for(unsigned int i = 0 ; i < m_NldTNld.size() ; ++ i)
    {
        MGamma += m_MGammafunc(facet, m_NLD[i])*m_NldTNld[i]*m_gaussWeightLD[i];
    }

    MGamma *= facet.getDetJ()*m_mesh.getRefElementSize(dim - 1);

    return MGamma;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, dim, 1> MatrixBuilder<dim, noPerEl>::getQN(const Facet& facet)
{
    Eigen::Matrix<double, dim, 1> qn; qn.setZero();

    std::array<double, 3> normal = facet.getNormal();
    Eigen::Matrix<double, dim, 1> n;

    if constexpr (dim == 2)
    {
        n[0] = normal[0];
        n[1] = normal[1];
    }
    else if constexpr (dim == 3)
    {
        n[0] = normal[0];
        n[1] = normal[1];
        n[2] = normal[2];
    }

    for(unsigned int i = 0 ; i < m_NLD.size() ; ++i)
    {
        Eigen::Matrix<double, dim, 1> q = m_QFunc(facet, m_gaussPointsLD[i]);
        double fact = (q.transpose()*n).value()*m_gaussWeightLD[i];
        qn += fact*m_NLD[i].transpose();
    }

    return qn*m_mesh.getRefElementSize(dim - 1)*facet.getDetJ();
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, dim*noPerEl, dim*noPerEl> MatrixBuilder<dim, noPerEl>::getK(const Element& element, const BmatType& B)
{
    Eigen::Matrix<double, dim*noPerEl, dim*noPerEl> K;
    double fact = 0;

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++i)
    {
        fact += m_Kfunc(element, m_NHD[i], B, m_ddev)*m_gaussWeightHD[i];
    }

    K = element.getDetJ()*m_mesh.getRefElementSize(dim)*fact*B.transpose()*m_ddev*B;

    return K;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, noPerEl, dim*noPerEl> MatrixBuilder<dim, noPerEl>::getD(const Element& element, const BmatType& B)
{
    Eigen::Matrix<double, noPerEl, dim*noPerEl> D;
    Eigen::Matrix<double, noPerEl, 1> sumWNT; sumWNT.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        sumWNT += m_Dfunc(element, m_NHD[i], B)*m_NHD[i].transpose()*m_gaussWeightHD[i];
    }

    D = element.getDetJ()*m_mesh.getRefElementSize(dim)*sumWNT*m_m.transpose()*B;

    return D;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, noPerEl, noPerEl> MatrixBuilder<dim, noPerEl>::getL(const Element& element, const BmatType& B, const GradNmatType& gradN)
{
    Eigen::Matrix<double, noPerEl, noPerEl> L;
    double fact = 0;

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        fact += m_Lfunc(element, m_NHD[i], B)*m_gaussWeightHD[i];
    }

    L = element.getDetJ()*m_mesh.getRefElementSize(dim)*fact*gradN.transpose()*gradN;

    return L;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, noPerEl, dim*noPerEl> MatrixBuilder<dim, noPerEl>::getC(const Element& element, const BmatType& B, const GradNmatType& gradN)
{
    Eigen::Matrix<double, noPerEl, dim*noPerEl>  C;
    Eigen::Matrix<double, dim, dim*noPerEl>  sumNW; sumNW.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        sumNW += m_Cfunc(element, m_NHD[i], B)*m_NHDtilde[i]*m_gaussWeightHD[i];
    }

    C = element.getDetJ()*m_mesh.getRefElementSize(dim)*gradN.transpose()*sumNW;

    return C;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, dim*noPerEl, 1> MatrixBuilder<dim, noPerEl>::getF(const Element& element,
                                                                 const Eigen::Matrix<double, dim, 1>& vec,
                                                                 const BmatType& B)
{
    Eigen::Matrix<double, dim*noPerEl, 1> F; F.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        F += m_Ffunc(element, m_NHD[i], B)*m_NHDtilde[i].transpose()*vec*m_gaussWeightHD[i];
    }

    F *= element.getDetJ()*m_mesh.getRefElementSize(dim);

    return F;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, noPerEl - 1, 1> MatrixBuilder<dim, noPerEl>::getSGamma(const Facet& facet)
{
    Eigen::Matrix<double, noPerEl - 1, 1> SGamma; SGamma.setZero();

    for(unsigned int i = 0 ; i < m_NLD.size() ; ++ i)
    {
        SGamma += m_SGammafunc(facet, m_NLD[i])*m_NLD[i].transpose()*m_gaussWeightLD[i];
    }

    SGamma *= facet.getDetJ()*m_mesh.getRefElementSize(dim - 1);

    return SGamma;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, noPerEl, 1> MatrixBuilder<dim, noPerEl>::getH(const Element& element, const Eigen::Matrix<double, dim, 1>& vec,
                                    const BmatType& B, const GradNmatType& gradN)
{
    Eigen::Matrix<double, noPerEl, 1> H; H.setZero();

    for(unsigned int i = 0 ; i < m_NHD.size() ; ++ i)
    {
        H += m_Hfunc(element, m_NHD[i], B)*gradN.transpose()*vec*m_gaussWeightHD[i];
    }

    H *= element.getDetJ()*m_mesh.getRefElementSize(dim);

    return H;
}

template<unsigned short dim, unsigned short noPerEl>
Eigen::Matrix<double, dim*noPerEl, 1> MatrixBuilder<dim, noPerEl>::getFST(const Facet& facet, const GradNmatType& gradNe, const BmatType& Be)
{
    Eigen::Matrix<double, dim*noPerEl, 1> FST; FST.setZero();

    for(unsigned int i = 0 ; i < m_NLD.size() ; ++ i)
    {
        Eigen::Matrix<double, dim*dim - 2*dim +3, 1> P = getP(facet);
        Eigen::Matrix<double, dim*dim - 2*dim +3, dim*dim - 2*dim +3> T = getT(P);
        FST -= m_FSTfunc(facet, m_NLD[i], m_NLDtilde[i], gradNe)*Be.transpose()*T*P*m_gaussWeightLD[i];
    }

    FST *= facet.getDetJ()*m_mesh.getRefElementSize(dim - 1);

    return FST;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setddev(DdevMatType ddev)
{
    m_ddev = ddev;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setm(mVecType m)
{
    m_m = m;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setMcomputeFactor(simpleMatFuncElm computeFactor)
{
    m_Mfunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setMGammacomputeFactor(simpleMatFuncFacet computeFactor)
{
    m_MGammafunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setKcomputeFactor(matFuncKElm computeFactor)
{
    m_Kfunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setDcomputeFactor(matFuncElm computeFactor)
{
    m_Dfunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setLcomputeFactor(matFuncElm computeFactor)
{
    m_Lfunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setCcomputeFactor(matFuncElm computeFactor)
{
    m_Cfunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setFcomputeFactor(matFuncElm computeFactor)
{
    m_Ffunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setSGammacomputeFactor(simpleMatFuncFacet computeFactor)
{
    m_SGammafunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setHcomputeFactor(matFuncElm computeFactor)
{
    m_Hfunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setFSTcomputeFactor(matFuncFacet computeFactor)
{
    m_FSTfunc = computeFactor;
}

template<unsigned short dim, unsigned short noPerEl>
void MatrixBuilder<dim, noPerEl>::setQFunc(qFuncFacet func)
{
    m_QFunc = func;
}
