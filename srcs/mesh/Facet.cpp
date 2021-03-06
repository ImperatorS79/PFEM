#include "Facet.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "Node.hpp"
#include "Mesh.hpp"

Facet::Facet(Mesh& mesh):
m_pMesh(&mesh)
{

}

void Facet::computeJ()
{
    m_J = {{{0, 0},
            {0, 0},
            {0, 0}}};

    if(m_nodesIndexes.size() == 2)
    {
        const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
        const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);

        double x0 = n0.getCoordinate(0);
        double x1 = n1.getCoordinate(0);
        double y0 = n0.getCoordinate(1);
        double y1 = n1.getCoordinate(1);

        m_J[0][0] = (x1 - x0)/2; //dx/dxi
        m_J[1][0] = (y1 - y0)/2; //dy/dxi
    }
    else
    {
        const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
        const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);
        const Node& n2 = m_pMesh->getNode(m_nodesIndexes[2]);

        double x0 = n0.getCoordinate(0);
        double x1 = n1.getCoordinate(0);
        double x2 = n2.getCoordinate(0);
        double y0 = n0.getCoordinate(1);
        double y1 = n1.getCoordinate(1);
        double y2 = n2.getCoordinate(1);
        double z0 = n0.getCoordinate(2);
        double z1 = n1.getCoordinate(2);
        double z2 = n2.getCoordinate(2);

        m_J[0][0] = x1 - x0; //dx/dxi
        m_J[1][0] = y1 - y0; //dy/dxi
        m_J[2][0] = z1 - z0; //dz/dxi
        m_J[0][1] = x2 - x0; //dx/deta
        m_J[1][1] = y2 - y0; //dy/deta
        m_J[2][1] = z2 - z0; //dz/deta
    }
}

void Facet::computeDetJ()
{
    if(m_nodesIndexes.size() == 2)
    {
        m_detJ = std::sqrt(m_J[0][0]*m_J[0][0] + m_J[1][0]*m_J[1][0]); // || (dx/dxi; dy/dxi)||
    }
    else
    {
        std::array<double, 3> dGammaVec = {
            m_J[1][0]*m_J[2][1] - m_J[1][1]*m_J[2][0],
            m_J[2][0]*m_J[0][1] - m_J[2][1]*m_J[0][0],
            m_J[0][0]*m_J[1][1] - m_J[1][0]*m_J[0][1]
        };

        // || (dx/dxi; dy/dxi; dz/dxi) x (dx/deta; dy/deta; dz/deta)||
        m_detJ = std::sqrt(dGammaVec[0]*dGammaVec[0] + dGammaVec[1]*dGammaVec[1] + dGammaVec[2]*dGammaVec[2]);
    }
}

void Facet::computeInvJ()
{
    m_invJ = {{{0, 0, 0},
               {0, 0, 0}}};

    if(m_nodesIndexes.size() == 2)
    {
        const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
        const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);

        double x0 = n0.getCoordinate(0);
        double x1 = n1.getCoordinate(0);
        double y0 = n0.getCoordinate(1);
        double y1 = n1.getCoordinate(1);

        m_invJ[0][0] = 1/m_J[0][0]; //dxi/dx
        m_invJ[0][1] = m_invJ[0][0]*((x1-x0)/(y1-y0)); //dxi/dy = dxi/dx*dx/dy
    }
    else
    {
        const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
        const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);
        const Node& n2 = m_pMesh->getNode(m_nodesIndexes[2]);

        double x0 = n0.getCoordinate(0);
        double x1 = n1.getCoordinate(0);
        double x2 = n2.getCoordinate(0);
        double y0 = n0.getCoordinate(1);
        double y1 = n1.getCoordinate(1);
        double y2 = n2.getCoordinate(1);
        double z0 = n0.getCoordinate(2);
        double z1 = n1.getCoordinate(2);
        double z2 = n2.getCoordinate(2);

        double detB = m_J[0][0]*m_J[1][1] - m_J[0][1]*m_J[1][0];

        m_invJ[0][0] = m_J[1][1]/detB;  //dxi/dx
        m_invJ[0][1] = -m_J[0][1]/detB; //dxi/dy
        m_invJ[1][0] = -m_J[1][0]/detB; //deta/dx
        m_invJ[1][1] = m_J[0][0]/detB;  //deta/dy

        double den = x0*(y1 - y2) + x1*(y2 - y0) + x2*(y0 - y1);

        double a_1 = den/(z0*(y1 - y2) + z1*(y1 - y0) + z2*(y0 - y1));
        double b_1 = den/(z0*(x2 - x1) + z1*(x0 - x2) + z2*(x1 - x0));

        m_invJ[0][2] = m_invJ[0][0]*a_1 + m_J[0][1]*b_1; //dxi/dz = dxi/dx*dx/dz + dxi/dy*dy/dz
        m_invJ[1][2] = m_invJ[1][0]*a_1 + m_J[1][1]*b_1; //deta/dz = deta/dx*dx/dz + deta/dy*dy/dz
    }
}

void Facet::computeNormal()
{
    if(m_pMesh->getDim() == 2)
    {
        const Node& node = getNode(0);
        const Node& node1 = getNode(1);
        const Node& outNode = getOutNode();

        double x0 = node.getCoordinate(0);
        double x1 = node1.getCoordinate(0);

        double y0 = node.getCoordinate(1);
        double y1 = node1.getCoordinate(1);

        m_normal = {
            y1 - y0,
            x0 - x1,
            0
        };

        std::array<double, 3> vecToOutNode = {
            outNode.getCoordinate(0) - x0,
            outNode.getCoordinate(1) - y0,
            0
        };

        if(m_normal[0]*vecToOutNode[0] + m_normal[1]*vecToOutNode[1] > 0)
        {
            m_normal[0] *= -1;
            m_normal[1] *= -1;
        }

        double norm = std::sqrt(m_normal[0]*m_normal[0] + m_normal[1]*m_normal[1]);
        m_normal[0] /= norm;
        m_normal[1] /= norm;
    }
    else
    {
        const Node& node = getNode(0);
        const Node& node1 = getNode(1);
        const Node& node2 = getNode(2);
        const Node& outNode = getOutNode();

        std::array<double, 3> a = {
            node1.getCoordinate(0) - node.getCoordinate(0),
            node1.getCoordinate(1) - node.getCoordinate(1),
            node1.getCoordinate(2) - node.getCoordinate(2)
        };

        std::array<double, 3> b = {
            node2.getCoordinate(0) - node.getCoordinate(0),
            node2.getCoordinate(1) - node.getCoordinate(1),
            node2.getCoordinate(2) - node.getCoordinate(2)
        };

        m_normal = {
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]
        };

        std::array<double, 3> vecToOutNode = {
            outNode.getCoordinate(0) - node.getCoordinate(0),
            outNode.getCoordinate(1) - node.getCoordinate(1),
            outNode.getCoordinate(2) - node.getCoordinate(2)
        };

        if(vecToOutNode[0]*m_normal[0] + vecToOutNode[1]*m_normal[1] + vecToOutNode[2]*m_normal[2] > 0)
        {
            m_normal[0] *= -1.0;
            m_normal[1] *= -1.0;
            m_normal[2] *= -1.0;
        }

        double norm = std::sqrt(m_normal[0]*m_normal[0] + m_normal[1]*m_normal[1] + m_normal[2]*m_normal[2]);

        m_normal[0] /= norm;
        m_normal[1] /= norm;
        m_normal[2] /= norm;
    }
}

const Node& Facet::getNode(unsigned int nodeIndex) const noexcept
{
    return m_pMesh->getNode(m_nodesIndexes[nodeIndex]);
}

const Node& Facet::getOutNode() const noexcept
{
    return m_pMesh->getNode(m_outNodeIndex);
}

const Element& Facet::getElement() const noexcept
{
    return m_pMesh->getElement(m_elementIndex);
}

std::array<double, 3> Facet::getPosFromGP(const std::array<double, 3>& gp) const noexcept
{
    const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);

    std::array<double, 3> pos;
    pos[0] = m_J[0][0]*gp[0] + m_J[0][1]*gp[1] + n0.getCoordinate(0);
    pos[1] = m_J[1][0]*gp[0] + m_J[1][1]*gp[1] + n0.getCoordinate(1);
    pos[2] = m_J[2][0]*gp[0] + m_J[2][1]*gp[1] + n0.getCoordinate(2);

    return pos;
}

std::vector<double> Facet::getState(unsigned int stateIndex) const noexcept
{
    std::vector<double> states(m_nodesIndexes.size());

    for(std::size_t i = 0 ; i < states.size() ; ++i)
        states[i] = m_pMesh->getNode(m_nodesIndexes[i]).getState(stateIndex);

    return states;
}

bool Facet::isOnFreeSurface() const noexcept
{
    if(m_pMesh->getDim() == 2)
        return getNode(0).isOnFreeSurface() && getNode(1).isOnFreeSurface();
    else
        return getNode(0).isOnFreeSurface() && getNode(1).isOnFreeSurface() && getNode(2).isOnFreeSurface();
}
