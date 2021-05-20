#pragma once
#ifndef STATESFROMQ_HPP_INCLUDED
#define STATESFROMQ_HPP_INCLUDED

#include <cassert>
#include "../../mesh/Mesh.hpp"
#include <Eigen/Dense>

inline void setNodesStatesfromQ(Mesh* pMesh, const Eigen::VectorXd& q, unsigned int beginState, unsigned int endState)
{
    assert(static_cast<std::size_t>(q.rows()) == (endState - beginState + 1)*pMesh->getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < pMesh->getNodesCount() ; ++n)
    {
        for (unsigned int s = beginState ; s <= endState ; ++s)
            pMesh->setNodeState(n, s, q((s - beginState)*pMesh->getNodesCount() + n));
    }
}

inline Eigen::VectorXd getQFromNodesStates(Mesh* pMesh, unsigned int beginState, unsigned int endState)
{
    Eigen::VectorXd q((endState - beginState + 1)*pMesh->getNodesCount());

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < pMesh->getNodesCount() ; ++n)
    {
        const Node& node = pMesh->getNode(n);
        for (unsigned int s = beginState ; s <= endState ; ++s)
            q((s - beginState)*pMesh->getNodesCount() + n) = node.getState(s);
    }

    return q;
}

template<unsigned short dim>
inline Eigen::Matrix<double, dim + 1, 1> getElementState(const Element& element, unsigned int state)
{
    Eigen::Matrix<double, dim + 1, 1> stateVec;

    if constexpr (dim == 2)
    {
        stateVec[0] = element.getNode(0).getState(state);
        stateVec[1] = element.getNode(1).getState(state);
        stateVec[2] = element.getNode(2).getState(state);
    }
    else if constexpr (dim == 3)
    {
        stateVec[0] = element.getNode(0).getState(state);
        stateVec[1] = element.getNode(1).getState(state);
        stateVec[2] = element.getNode(2).getState(state);
        stateVec[3] = element.getNode(3).getState(state);
    }

    return stateVec;
}

template<unsigned short dim>
inline Eigen::Matrix<double, dim + 1, 1> getElementState(const Eigen::VectorXd& q, const Element& element, unsigned short beginState, std::size_t nNodes)
{
    Eigen::Matrix<double, dim + 1, 1> stateVec;

    if constexpr (dim == 2)
    {
        stateVec[0] = q[element.getNodeIndex(0) + beginState*nNodes];
        stateVec[1] = q[element.getNodeIndex(1) + beginState*nNodes];
        stateVec[2] = q[element.getNodeIndex(2) + beginState*nNodes];
    }
    else if constexpr (dim == 3)
    {
        stateVec[0] = q[element.getNodeIndex(0) + beginState*nNodes];
        stateVec[1] = q[element.getNodeIndex(1) + beginState*nNodes];
        stateVec[2] = q[element.getNodeIndex(2) + beginState*nNodes];
        stateVec[3] = q[element.getNodeIndex(3) + beginState*nNodes];
    }

    return stateVec;
}

template<unsigned short dim>
inline Eigen::Matrix<double, dim*(dim + 1), 1> getElementVecState(const Element& element, unsigned int beginState)
{
    Eigen::Matrix<double, dim*(dim + 1), 1> stateVec;

    if constexpr (dim == 2)
    {
        stateVec[0] = element.getNode(0).getState(beginState);
        stateVec[1] = element.getNode(1).getState(beginState);
        stateVec[2] = element.getNode(2).getState(beginState);

        stateVec[3] = element.getNode(0).getState(beginState + 1);
        stateVec[4] = element.getNode(1).getState(beginState + 1);
        stateVec[5] = element.getNode(2).getState(beginState + 1);
    }
    else if constexpr (dim == 3)
    {
        stateVec[0] = element.getNode(0).getState(beginState);
        stateVec[1] = element.getNode(1).getState(beginState);
        stateVec[2] = element.getNode(2).getState(beginState);
        stateVec[3] = element.getNode(3).getState(beginState);

        stateVec[4] = element.getNode(0).getState(beginState + 1);
        stateVec[5] = element.getNode(1).getState(beginState + 1);
        stateVec[6] = element.getNode(2).getState(beginState + 1);
        stateVec[7] = element.getNode(3).getState(beginState + 1);

        stateVec[8] = element.getNode(0).getState(beginState + 2);
        stateVec[9] = element.getNode(1).getState(beginState + 2);
        stateVec[10] = element.getNode(2).getState(beginState + 2);
        stateVec[11] = element.getNode(3).getState(beginState + 2);
    }

    return stateVec;
}

template<unsigned short dim>
inline Eigen::Matrix<double, dim*(dim + 1), 1> getElementVecState(const Eigen::VectorXd& q, const Element& element, unsigned short beginState, std::size_t nNodes)
{
    Eigen::Matrix<double, dim*(dim + 1), 1> stateVec;

    if constexpr (dim == 2)
    {
        stateVec[0] = q[element.getNodeIndex(0) + beginState*nNodes];
        stateVec[1] = q[element.getNodeIndex(1) + beginState*nNodes];
        stateVec[2] = q[element.getNodeIndex(2) + beginState*nNodes];

        stateVec[3] = q[element.getNodeIndex(0) + (beginState + 1)*nNodes];
        stateVec[4] = q[element.getNodeIndex(1) + (beginState + 1)*nNodes];
        stateVec[5] = q[element.getNodeIndex(2) + (beginState + 1)*nNodes];
    }
    else if constexpr (dim == 3)
    {
        stateVec[0] = q[element.getNodeIndex(0) + beginState*nNodes];
        stateVec[1] = q[element.getNodeIndex(1) + beginState*nNodes];
        stateVec[2] = q[element.getNodeIndex(2) + beginState*nNodes];
        stateVec[3] = q[element.getNodeIndex(3) + beginState*nNodes];

        stateVec[4] = q[element.getNodeIndex(0) + (beginState + 1)*nNodes];
        stateVec[5] = q[element.getNodeIndex(1) + (beginState + 1)*nNodes];
        stateVec[6] = q[element.getNodeIndex(2) + (beginState + 1)*nNodes];
        stateVec[7] = q[element.getNodeIndex(3) + (beginState + 1)*nNodes];

        stateVec[8] = q[element.getNodeIndex(0) + (beginState + 2)*nNodes];
        stateVec[9] = q[element.getNodeIndex(1) + (beginState + 2)*nNodes];
        stateVec[10] = q[element.getNodeIndex(2) + (beginState + 2)*nNodes];
        stateVec[11] = q[element.getNodeIndex(3) + (beginState + 2)*nNodes];
    }

    return stateVec;
}

#endif
