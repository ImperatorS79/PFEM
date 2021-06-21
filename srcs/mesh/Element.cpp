#include "Element.hpp"

#include <cassert>
#include <Eigen/Dense>

#include "Node.hpp"
#include "Mesh.hpp"

Element::Element(Mesh& mesh):
m_pMesh(&mesh)
{

}

void Element::computeJ()
{
    m_J = {{{0, 0, 0},
            {0, 0, 0},
            {0, 0, 0}}};

    if(m_nodesIndexes.size() == 3)
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

        m_J[0][0] = x1 - x0;
        m_J[0][1] = x2 - x0;
        m_J[1][0] = y1 - y0;
        m_J[1][1] = y2 - y0;
    }
    else
    {
        const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);
        const Node& n1 = m_pMesh->getNode(m_nodesIndexes[1]);
        const Node& n2 = m_pMesh->getNode(m_nodesIndexes[2]);
        const Node& n3 = m_pMesh->getNode(m_nodesIndexes[3]);

        double x0 = n0.getCoordinate(0);
        double x1 = n1.getCoordinate(0);
        double x2 = n2.getCoordinate(0);
        double x3 = n3.getCoordinate(0);
        double y0 = n0.getCoordinate(1);
        double y1 = n1.getCoordinate(1);
        double y2 = n2.getCoordinate(1);
        double y3 = n3.getCoordinate(1);
        double z0 = n0.getCoordinate(2);
        double z1 = n1.getCoordinate(2);
        double z2 = n2.getCoordinate(2);
        double z3 = n3.getCoordinate(2);

        m_J[0][0] = x1 - x0;
        m_J[0][1] = x2 - x0;
        m_J[0][2] = x3 - x0;
        m_J[1][0] = y1 - y0;
        m_J[1][1] = y2 - y0;
        m_J[1][2] = y3 - y0;
        m_J[2][0] = z1 - z0;
        m_J[2][1] = z2 - z0;
        m_J[2][2] = z3 - z0;
    }
}

void Element::computeDetJ()
{
    if(m_nodesIndexes.size() == 3)
    {
        m_detJ = m_J[0][0]*m_J[1][1] - m_J[1][0]*m_J[0][1];
    }
    else
    {
        m_detJ = m_J[0][0]*m_J[1][1]*m_J[2][2]
               + m_J[0][1]*m_J[1][2]*m_J[2][0]
               + m_J[0][2]*m_J[1][0]*m_J[2][1]
               - m_J[2][0]*m_J[1][1]*m_J[0][2]
               - m_J[2][1]*m_J[1][2]*m_J[0][0]
               - m_J[2][2]*m_J[1][0]*m_J[0][1];
    }
}

void Element::computeInvJ()
{
    assert(m_detJ != 0);

    m_invJ = {{{0, 0, 0},
               {0, 0, 0},
               {0, 0, 0}}};

    if(m_nodesIndexes.size() == 3)
    {
        m_invJ[0][0] = m_J[1][1]/m_detJ;

        m_invJ[0][1] = - m_J[0][1]/m_detJ;

        m_invJ[1][0] = - m_J[1][0]/m_detJ;

        m_invJ[1][1] = m_J[0][0]/m_detJ;
    }
    else
    {
        m_invJ[0][0] = (m_J[1][1]*m_J[2][2]
                     - m_J[1][2]*m_J[2][1])/m_detJ;

        m_invJ[0][1] = (m_J[2][1]*m_J[0][2]
                     - m_J[2][2]*m_J[0][1])/m_detJ;

        m_invJ[0][2] = (m_J[0][1]*m_J[1][2]
                     - m_J[0][2]*m_J[1][1])/m_detJ;

        m_invJ[1][0] = (m_J[2][0]*m_J[1][2]
                     - m_J[1][0]*m_J[2][2])/m_detJ;

        m_invJ[1][1] = (m_J[0][0]*m_J[2][2]
                     - m_J[2][0]*m_J[0][2])/m_detJ;

        m_invJ[1][2] = (m_J[1][0]*m_J[0][2]
                     - m_J[0][0]*m_J[1][2])/m_detJ;

        m_invJ[2][0] = (m_J[1][0]*m_J[2][1]
                     - m_J[2][0]*m_J[1][1])/m_detJ;

        m_invJ[2][1] = (m_J[2][0]*m_J[0][1]
                     - m_J[0][0]*m_J[2][1])/m_detJ;

        m_invJ[2][2] = (m_J[0][0]*m_J[1][1]
                     - m_J[1][0]*m_J[0][1])/m_detJ;
    }
}

const Element& Element::getNeighbourElement(unsigned int neighbourElmIndex) const noexcept
{
    return m_pMesh->getElement(m_neighbourElements[neighbourElmIndex]);
}

const Node& Element::getNode(unsigned int nodeIndex) const noexcept
{
    return m_pMesh->getNode(m_nodesIndexes[nodeIndex]);
}

std::array<double, 3> Element::getPosFromGP(const std::array<double, 3>& gp) const noexcept
{
    const Node& n0 = m_pMesh->getNode(m_nodesIndexes[0]);

    std::array<double, 3> pos;
    pos[0] = m_J[0][0]*gp[0] + m_J[0][1]*gp[1] + m_J[0][2]*gp[2] + n0.getCoordinate(0);
    pos[1] = m_J[1][0]*gp[0] + m_J[1][1]*gp[1] + m_J[1][2]*gp[2] + n0.getCoordinate(1);
    pos[2] = m_J[2][0]*gp[0] + m_J[2][1]*gp[1] + m_J[2][2]*gp[2] + n0.getCoordinate(2);

    return pos;
}

double Element::getSize() const noexcept
{
    return m_detJ*m_pMesh->getRefElementSize(m_pMesh->getDim());
}

double Element::getMinNodeDist() const noexcept
{
    if(m_pMesh->getDim() == 2)
    {
        const Node& n0 = this->getNode(0);
        const Node& n1 = this->getNode(1);
        const Node& n2 = this->getNode(2);

        std::array<double, 3> dist;
        dist[0] = Node::distance(n0, n1);
        dist[1] = Node::distance(n1, n2);
        dist[2] = Node::distance(n0, n2);

        return *std::min_element(dist.begin(), dist.end());
    }
    else
    {
        const Node& n0 = this->getNode(0);
        const Node& n1 = this->getNode(1);
        const Node& n2 = this->getNode(2);
        const Node& n3 = this->getNode(3);

        std::array<double, 6> dist;
        dist[0] = Node::distance(n0, n1);
        dist[1] = Node::distance(n1, n2);
        dist[2] = Node::distance(n0, n2);
        dist[3] = Node::distance(n0, n3);
        dist[4] = Node::distance(n1, n3);
        dist[5] = Node::distance(n2, n3);

        return *std::min_element(dist.begin(), dist.end());
    }
}

double Element::getRin() const noexcept
{
    if(m_pMesh->getDim() == 2)
    {
        double A = getSize();
        const Node& n0 = this->getNode(0);
        const Node& n1 = this->getNode(1);
        const Node& n2 = this->getNode(2);

        double a = Node::distance(n0, n1);
        double b = Node::distance(n1, n2);
        double c = Node::distance(n0, n2);
        double s = (a + b + c)/2;

        return A/s;
    }
    else
    {
        double V = getSize();
        const Node& n0 = this->getNode(0);
        const Node& n1 = this->getNode(1);
        const Node& n2 = this->getNode(2);
        const Node& n3 = this->getNode(3);

        double x0 = n0.getCoordinate(0);
        double x1 = n1.getCoordinate(0);
        double x2 = n2.getCoordinate(0);
        double x3 = n3.getCoordinate(0);
        double y0 = n0.getCoordinate(1);
        double y1 = n1.getCoordinate(1);
        double y2 = n2.getCoordinate(1);
        double y3 = n3.getCoordinate(1);
        double z0 = n0.getCoordinate(2);
        double z1 = n1.getCoordinate(2);
        double z2 = n2.getCoordinate(2);
        double z3 = n3.getCoordinate(2);

        std::array<double, 3> nf1 = {
            (y1 - y0)*(z2 - z0) - (z1 - z0)*(y2 - y0),
            (z1 - z0)*(x2 - x0) - (x1 - x0)*(z2 - z0),
            (x1 - x0)*(y2 - y0) - (y1 - y0)*(x2 - x0)
        };
        double normnf1 = std::sqrt(nf1[0]*nf1[0] + nf1[1]*nf1[1] + nf1[2]*nf1[2]);

        std::array<double, 3> nf2 = {
            (y3 - y0)*(z2 - z0) - (z3 - z0)*(y2 - y0),
            (z3 - z0)*(x2 - x0) - (x3 - x0)*(z2 - z0),
            (x3 - x0)*(y2 - y0) - (y3 - y0)*(x2 - x0)
        };
        double normnf2 = std::sqrt(nf2[0]*nf2[0] + nf2[1]*nf2[1] + nf2[2]*nf2[2]);

        std::array<double, 3> nf3 = {
            (y1 - y0)*(z3 - z0) - (z1 - z0)*(y3 - y0),
            (z1 - z0)*(x3 - x0) - (x1 - x0)*(z3 - z0),
            (x1 - x0)*(y3 - y0) - (y1 - y0)*(x3 - x0)
        };
        double normnf3 = std::sqrt(nf3[0]*nf3[0] + nf3[1]*nf3[1] + nf3[2]*nf3[2]);

        std::array<double, 3> nf4 = {
            (y1 - y3)*(z2 - z3) - (z1 - z3)*(y2 - y3),
            (z1 - z3)*(x2 - x3) - (x1 - x3)*(z2 - z3),
            (x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3)
        };
        double normnf4 = std::sqrt(nf4[0]*nf4[0] + nf4[1]*nf4[1] + nf4[2]*nf4[2]);


        return 6*V/(normnf1 + normnf2 + normnf3 + normnf4);
    }
}

std::vector<double> Element::getState(unsigned int stateIndex) const noexcept
{
    std::vector<double> states(m_nodesIndexes.size());

    for(std::size_t i = 0 ; i < states.size() ; ++i)
        states[i] = m_pMesh->getNode(m_nodesIndexes[i]).getState(stateIndex);

    return states;
}

bool Element::isContact() const noexcept
{
    for(std::size_t i = 0 ; i < m_nodesIndexes.size() ; ++i)
    {
        const Node& node = m_pMesh->getNode(m_nodesIndexes[i]);
        if(node.isBound())
            return true;
    }

    return false;
}

bool Element::isOnFS() const noexcept
{
    for(std::size_t i = 0 ; i < m_nodesIndexes.size() ; ++i)
    {
        const Node& node = m_pMesh->getNode(m_nodesIndexes[i]);
        if(node.isOnFreeSurface())
            return true;
    }

    return false;
}

