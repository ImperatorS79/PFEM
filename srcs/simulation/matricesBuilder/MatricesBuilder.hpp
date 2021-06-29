#pragma once
#ifndef MATRIXBUILDER_HPP_INCLUDED
#define MATRIXBUILDER_HPP_INCLUDED

#include <functional>
#include <Eigen/Dense>

class Element;
class Facet;
class Mesh;

template<unsigned short dim, unsigned short noPerEl = dim + 1>
using NmatTypeHD = Eigen::Matrix<double, 1, noPerEl>;
template<unsigned short dim, unsigned short noPerEl = dim + 1>
using NmatTildeTypeHD = Eigen::Matrix<double, dim, dim*noPerEl>;
template<unsigned short dim, unsigned short noPerEl = dim + 1>
using NmatTypeLD = Eigen::Matrix<double, 1, noPerEl - 1>;
template<unsigned short dim, unsigned short noPerEl = dim + 1>
using NmatTildeTypeLD = Eigen::Matrix<double, dim, dim*(noPerEl - 1)>;
template<unsigned short dim, unsigned short noPerEl = dim + 1>
using BmatType = Eigen::Matrix<double, dim*dim - 2*dim +3, dim*noPerEl>;
template<unsigned short dim, unsigned short noPerEl = dim + 1>
using GradNmatType = Eigen::Matrix<double, dim, noPerEl> ;
template<unsigned short dim, unsigned short noPerEl = dim + 1>
using DdevMatType = Eigen::Matrix<double, dim*dim - 2*dim +3, dim*dim - 2*dim +3>;
template<unsigned short dim, unsigned short noPerEl = dim + 1>
using mVecType = Eigen::Matrix<double, dim*dim - 2*dim + 3, 1>;

/**
 * \class MatrixBuilder
 * \brief Class responsible to hold code for building matrices.
 */
template<unsigned short dim, unsigned short noPerEl = dim + 1>
class MatrixBuilder
{

    using NmatTypeHD = Eigen::Matrix<double, 1, noPerEl>;
    using NmatTildeTypeHD = Eigen::Matrix<double, dim, dim*noPerEl>;
    using NmatTypeLD = Eigen::Matrix<double, 1, noPerEl - 1>;
    using NmatTildeTypeLD = Eigen::Matrix<double, dim, dim*(noPerEl - 1)>;
    using BmatType = Eigen::Matrix<double, dim*dim - 2*dim + 3, dim*(noPerEl)>;
    using GradNmatType = Eigen::Matrix<double, dim, noPerEl> ;
    using DdevMatType = Eigen::Matrix<double, dim*dim - 2*dim +3, dim*dim - 2*dim +3>;
    using mVecType = Eigen::Matrix<double, dim*dim - 2*dim + 3, 1>;
    using matFuncElm = std::function<double(const Element&, const NmatTypeHD&, const BmatType&)>;
    using matFuncKElm = std::function<double(const Element&, const NmatTypeHD&, const BmatType&, const DdevMatType&)>;
    using simpleMatFuncElm = std::function<double(const Element&, const NmatTypeHD&)>;
    using simpleMatFuncFacet = std::function<double(const Facet&, const NmatTypeLD&)>;
    using matFuncFacet = std::function<double(const Facet&, const NmatTypeLD&,  const NmatTildeTypeLD&, const GradNmatType&)>;
    using matVecFuncFacet = std::function<Eigen::Matrix<double, dim*dim - 2*dim + 3, 1>(const Facet&, const NmatTypeLD&,  const NmatTildeTypeLD&, const BmatType&)>;
    using qFuncFacet = std::function<Eigen::Matrix<double, dim, 1>(const Facet&, const std::array<double, 3>& /** gp **/)>;

	public:
		MatrixBuilder(const Mesh& mesh, unsigned int nGPHD, unsigned int nGPLD);
		~MatrixBuilder();

		BmatType getB(const GradNmatType& gradN);
		GradNmatType getGradN(const Element& element);
		Eigen::Matrix<double, noPerEl, noPerEl>                 getM(const Element& element);
		Eigen::Matrix<double, noPerEl - 1, noPerEl - 1>         getMGamma(const Facet& facet);
		Eigen::Matrix<double, dim*noPerEl, dim*noPerEl>         getK(const Element& element, const BmatType& B);
		Eigen::Matrix<double, noPerEl, dim*noPerEl>             getD(const Element& element, const BmatType& B);
		Eigen::Matrix<double, noPerEl, noPerEl>                 getL(const Element& element, const BmatType& B, const GradNmatType& gradN);
		Eigen::Matrix<double, noPerEl, dim*noPerEl>             getC(const Element& element, const BmatType& B, const GradNmatType& gradN);
		Eigen::Matrix<double, dim*noPerEl, 1>                   getF(const Element& element, const Eigen::Matrix<double, dim, 1>& vec, const BmatType& B);
		Eigen::Matrix<double, noPerEl - 1, 1>                   getSGamma(const Facet& facet);
		Eigen::Matrix<double, noPerEl, 1>                       getH(const Element& element, const Eigen::Matrix<double, dim, 1>& vec, const BmatType& B, const GradNmatType& gradN);
		Eigen::Matrix<double, dim, 1>                           getQN(const Facet& facet);
		Eigen::Matrix<double, dim*dim - 2*dim +3, 1>            getP(const Facet& facet);
        Eigen::Matrix<double, dim*dim - 2*dim +3, dim*dim - 2*dim +3> getT(const Eigen::Matrix<double, dim*dim - 2*dim +3, 1>& P);
        Eigen::Matrix<double, dim*noPerEl, 1>                   getFST(const Facet& facet, const GradNmatType& gradNe, const BmatType& Be);

		void setddev(DdevMatType ddev);
		void setm(mVecType m);
		void setMcomputeFactor(simpleMatFuncElm computeFactor);
		void setMGammacomputeFactor(simpleMatFuncFacet computeFactor);
		void setKcomputeFactor(matFuncKElm computeFactor);
		void setDcomputeFactor(matFuncElm computeFactor);
		void setLcomputeFactor(matFuncElm computeFactor);
		void setCcomputeFactor(matFuncElm computeFactor);
		void setFcomputeFactor(matFuncElm computeFactor);
		void setSGammacomputeFactor(simpleMatFuncFacet computeFactor);
		void setHcomputeFactor(matFuncElm computeFactor);
		void setFSTcomputeFactor(matFuncFacet computeFactor);
		void setQFunc(qFuncFacet func);

		template <unsigned short Size>
		static Eigen::DiagonalMatrix<double, Size> lump2(const Eigen::Matrix<double, Size, Size>& mat)
		{
            Eigen::DiagonalMatrix<double, Size> lumpedMat;
            lumpedMat.setZero();

            auto& lumpedMatDiag = lumpedMat.diagonal();

            for(auto i = 0; i < mat.cols() ; ++i)
            {
                for(auto j = 0; j < mat.cols() ; ++j)
                    lumpedMatDiag(i) += mat(i,j);
            }

            return lumpedMat;
        }

        template <unsigned short Size>
        static void lump(Eigen::Matrix<double, Size, Size>& mat)
		{
            for(auto i = 0; i < mat.cols() ; ++i)
            {
                for(auto j = 0; j < mat.cols() ; ++j)
                {
                    if(i != j)
                    {
                        mat(i, i) += mat(i,j);
                        mat(i, j) = 0;
                    }
                }
            }
        }

		static Eigen::Matrix<double, dim*noPerEl, dim*noPerEl> diagBlock(const Eigen::Matrix<double, noPerEl, noPerEl>& mat)
        {
            Eigen::Matrix<double, dim*noPerEl, dim*noPerEl> finalMat;
            finalMat.setZero();

            for(unsigned int i = 0 ; i < dim ; ++i)
                finalMat.block(i*mat.rows(), i*mat.cols(), mat.rows(), mat.cols()) = mat;

            return finalMat;
        }

        static Eigen::Matrix<double, dim*dim, dim*dim> diagBlock(const Eigen::Matrix<double, dim, dim>& mat)
        {
            Eigen::Matrix<double, dim*dim, dim*dim> finalMat;
            finalMat.setZero();

            for(unsigned int i = 0 ; i < dim ; ++i)
                finalMat.block(i*mat.rows(), i*mat.cols(), mat.rows(), mat.cols()) = mat;

            return finalMat;
        }

        static void inverse(Eigen::DiagonalMatrix<double, Eigen::Dynamic>& mat)
        {
            auto& matDiag = mat.diagonal();

            #pragma omp parallel for default(shared)
            for(auto i = 0 ; i < mat.rows() ; ++i)
                matDiag[i] = 1/matDiag[i];
        }

    private:
        const Mesh& m_mesh;

        unsigned int m_nGPHD; /**< Number of Gauss points use for the elements. **/
        std::vector<double> m_gaussWeightHD; /**< Gauss weights for the elements. **/
        std::vector<std::array<double, 3>> m_gaussPointsHD; /**< Gauss points for the elements. **/
        std::vector<std::vector<double>> m_sfsHD; /**< Element shape functions evaluated at each Gauss Points for elements. **/
        std::vector<std::vector<double>> m_gradsfsHD; /**< Element gradient of shape functions for elements (constant). **/

        unsigned int m_nGPLD; /**< Number of Gauss points use for the facets. **/
        std::vector<double> m_gaussWeightLD; /**< Gauss weights for the facets. **/
        std::vector<std::array<double, 3>> m_gaussPointsLD; /**< Gauss points for the facets. **/
        std::vector<std::vector<double>> m_sfsLD; /**< Element shape functions evaluated at each Gauss Points for facets. **/
        std::vector<std::vector<double>> m_gradsfsLD; /**< Element gradient of shape functions for facets (constant). **/

        std::vector<NmatTypeHD> m_NHD;
        std::vector<Eigen::Matrix<double, noPerEl, noPerEl>> m_NhdTNhd;
        Eigen::Matrix<double, noPerEl, noPerEl> m_sum_NhdTNhd_w;
        std::vector<NmatTypeLD> m_NLD;
        std::vector<Eigen::Matrix<double, noPerEl - 1, noPerEl - 1>> m_NldTNld;
        Eigen::Matrix<double, noPerEl - 1, noPerEl - 1> m_sum_NldTNld_w;
        std::vector<NmatTildeTypeHD> m_NHDtilde;
        std::vector<NmatTildeTypeLD> m_NLDtilde;

        DdevMatType m_ddev;
        mVecType m_m;

        simpleMatFuncElm m_Mfunc;
        simpleMatFuncFacet m_MGammafunc;
        matFuncKElm m_Kfunc;
        matFuncElm m_Dfunc;
        matFuncElm m_Lfunc;
        matFuncElm m_Cfunc;
        matFuncElm m_Ffunc;
        simpleMatFuncFacet m_SGammafunc;
        matFuncElm m_Hfunc;
        matFuncFacet m_FSTfunc;
        qFuncFacet m_QFunc;
};

#include "MatricesBuilder.inl"

#endif // MATRIXBUILDER_HPP_INCLUDED
