#pragma once
#ifndef PSMOOTHER_HPP_INCLUDED
#define PSMOOTHER_HPP_INCLUDED

#include <memory>

#include "../matricesBuilder/MatricesBuilder.hpp"
#include "MeshSmoother.hpp"

template<unsigned short dim>
class PSmoother : public MeshSmoother
{
    public:
        PSmoother(Problem* pProblem, Mesh& mesh, double timeBetweenSmooth, double a, double epsADRtoll, double betaInit);
        ~PSmoother() override;

        void smooth(bool verboseOuput) override;
    private:
        double m_a;
        double m_epsADRtoll;
        double m_beta;
        double m_dt = 1;

        std::unique_ptr<MatrixBuilder<dim>> m_pMatBuilder;

        std::vector<bool> m_elementsToKeep; //which element will be part of the smoothing
        std::vector<std::size_t> m_elmMeshToProb;

        std::vector<int> m_nodesMeshToBC; //-1: the node is not part of the problem, 0: node move, 1: node is fixed
        std::vector<std::size_t> m_nodeMeshToProb;
        std::vector<std::size_t> m_nodeProbToMesh;


        Eigen::DiagonalMatrix<double, Eigen::Dynamic> m_lhs;
        Eigen::VectorXd m_rhs;

        Eigen::VectorXd m_un1;
        Eigen::VectorXd m_vn12;
        Eigen::VectorXd m_vn_12;
        std::vector<Eigen::Matrix<double, dim*(dim + 1), 1>> m_Fint;
        std::vector<Eigen::Matrix<double, dim*(dim + 1), 1>> m_FintPrec;

        bool m_init(bool verboseOuput);
        void m_buildSystem(const Eigen::VectorXd& vBar, const Eigen::VectorXd& uBar);
        void m_applyBC();
};

#include "PSmoother.inl"

#endif // PSMOOTHER_HPP_INCLUDED
