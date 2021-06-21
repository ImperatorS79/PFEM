#include "Mesh.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>

#include <gmsh.h>


Mesh::Mesh(const MeshCreateInfo& meshInfos) :
m_hchar(meshInfos.hchar),
m_alpha(meshInfos.alpha),
m_omega(meshInfos.omega),
m_gamma(meshInfos.gamma),
m_boundingBox(meshInfos.boundingBox),
m_exclusionZones(meshInfos.exclusionZones),
m_addOnFS(meshInfos.addOnFS),
m_deleteFlyingNodes(meshInfos.deleteFlyingNodes),
m_computeNormalCurvature(true)
{
    loadFromFile(meshInfos.mshFile);
}

bool Mesh::addNodes(bool verboseOutput)
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh!");

    bool addedNodes = false;

    std::vector<bool> toBeDeleted(m_elementsList.size(), false);

    std::size_t elementCount = m_elementsList.size(); //avoid problems because the follwing for loop will increase the element number

    double limitSize =  m_omega*std::pow(m_hchar, m_dim);

    for(std::size_t elm = 0 ; elm < elementCount ; ++elm)
    {
        //If an element is too big, we add a node at his centre
        if(m_elementsList[elm].getSize() > limitSize && (m_addOnFS ? true : !m_elementsList[elm].isOnFS()))
        {
            Node newNode(*this);

            for(unsigned short k = 0 ; k < m_dim ; ++k)
            {
                newNode.m_position[k] = 0;
                for(unsigned short d = 0 ; d <= m_dim ; ++d)
                {
                    assert(m_elementsList[elm].m_nodesIndexes[d] < m_nodesList.size()) ;
                    newNode.m_position[k] += m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].m_position[k];
                }
                newNode.m_position[k] /= (m_dim + 1);
            }

            newNode.m_states.resize(m_nodesList[0].m_states.size());
            for(unsigned short k = 0 ; k < m_nodesList[0].m_states.size() ; ++k)
            {
                newNode.m_states[k] = 0;
                for(unsigned short d = 0 ; d <= m_dim ; ++d)
                {
                    newNode.m_states[k] += m_nodesList[m_elementsList[elm].m_nodesIndexes[d]].m_states[k];
                }
                newNode.m_states[k] /= (m_dim + 1);
            }

            m_nodesList.push_back(std::move(newNode));

            for(unsigned short d = 0 ; d <= m_dim ; ++d)
            {
                Element element(*this);
                switch(m_dim)
                {
                    case 2:
                        if(d == m_dim)
                        {
                            element.m_nodesIndexes = {m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[0],
                                                      m_nodesList.size() - 1};
                        }
                        else
                        {
                            element.m_nodesIndexes = {m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[d + 1],
                                                      m_nodesList.size() - 1};
                        }
                        break;

                    default:
                        if(d == m_dim)
                        {
                            element.m_nodesIndexes = {m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[0],
                                                      m_elementsList[elm].m_nodesIndexes[1],
                                                      m_nodesList.size() - 1};
                        }
                        else if (d == m_dim - 1)
                        {
                            element.m_nodesIndexes = {m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[d + 1],
                                                      m_elementsList[elm].m_nodesIndexes[0],
                                                      m_nodesList.size() - 1};
                        }
                        else
                        {
                            element.m_nodesIndexes = {m_elementsList[elm].m_nodesIndexes[d],
                                                      m_elementsList[elm].m_nodesIndexes[d + 1],
                                                      m_elementsList[elm].m_nodesIndexes[d + 2],
                                                      m_nodesList.size() - 1};
                        }
                        break;
                }

                //No recompute of J, detJ and ivJ, everything is recomputed in triangulateAndAlphaShape
                m_elementsList.push_back(std::move(element));
            }

            toBeDeleted[elm] = true;

            if(verboseOutput)
            {
                std::cout << "Adding node (" << "(";
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    std::cout << m_nodesList[m_nodesList.size() - 1].m_position[d];
                    if(d == m_dim - 1)
                        std::cout << ")";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
            }

            addedNodes = true;
        }
    }

    m_elementsList.erase(std::remove_if(m_elementsList.begin(), m_elementsList.end(), [this, &toBeDeleted](const Element& element){
        if(toBeDeleted[&element - &*std::begin(m_elementsList)])
            return true;
        else
            return false;
    }), m_elementsList.end());

    return addedNodes;
}

bool Mesh::checkBoundingBox(bool verboseOutput) noexcept
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh !");

    std::vector<bool> toBeDeletedNodes(m_nodesList.size(), false);   //Should the node be deleted

    bool outofBBNodes = false;

    auto isNodeOutOfBB = [this](const Node& node) -> bool {
        for(unsigned short d = 0 ; d < m_dim ; ++d)
        {
            if(node.m_position[d] < m_boundingBox[d] ||
               node.m_position[d] > m_boundingBox[d + m_dim])
            {
                return true;
            }
        }

        return false;
    };

    //If the whole element is out of the bounding box, we delete it.
    // Bounding box fromat: [xmin, ymin, zmin, xmax, ymax, zmax]
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(isNodeOutOfBB(m_nodesList[n]))
        {
            toBeDeletedNodes[n] = true;
            outofBBNodes = true;
        }
    }

    m_nodesList.erase(
    std::remove_if(m_nodesList.begin(), m_nodesList.end(), [this, &toBeDeletedNodes, verboseOutput](const Node& node)
    {
        if(toBeDeletedNodes[&node - &*std::begin(m_nodesList)])
        {
            if(verboseOutput)
            {
                std::cout << "Removing out of bounding box node (";
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    std::cout << node.m_position[d];
                    if(d == m_dim - 1)
                        std::cout << ")";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
            }

            return true;
        }
        else
           return false;
    }), m_nodesList.end());

    return outofBBNodes;
}

void Mesh::computeMeshDim()
{
    int elementDim = -1;

    // loop over the dimension i to get the maximum element dimension in the mesh
    for(unsigned short i = 2 ; i <= 3 ; ++i)
    {
        std::vector<int> eleTypes;
        gmsh::model::mesh::getElementTypes(eleTypes, i);

        switch(eleTypes.size())
        {
            case 0:
                break;
            case 1:
                elementDim = i;
                break;
            default:
                throw std::runtime_error("do not use hybrid meshes for PFEM simulations!");
        }
    }

    if(elementDim == -1)
        throw std::runtime_error("there is no suitable elements in the .msh file!");
    else
        m_dim = static_cast<unsigned short>(elementDim);
}

void Mesh::computeFSNormalCurvature()
{
    switch(m_dim)
    {
        case 2:
            computeFSNormalCurvature2D();
            break;
        default:

            computeFSNormalCurvature3D();
            break;
    }
}

void Mesh::deleteFlyingNodes(bool verboseOutput) noexcept
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh !");

    std::vector<bool> toBeDeletedNodes(m_nodesList.size(), false);   //Should the node be deleted

	std::size_t nodesToDeleteCount = 0;
    std::vector<std::size_t> nodesIndexesDeleted = {};

    //If the whole element is out of the bounding box, we delete it.
    // Bounding box fromat: [xmin, ymin, zmin, xmax, ymax, zmax]
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(m_nodesList[n].isFree() && !m_nodesList[n].isBound())
        {
            toBeDeletedNodes[n] = true;
			nodesToDeleteCount++;
        }
    }

	nodesIndexesDeleted.resize(nodesToDeleteCount);
	std::size_t count = 0;
    for(std::size_t n = 0 ; n < toBeDeletedNodes.size() ; ++n)
    {
		if(toBeDeletedNodes[n])
		{
			nodesIndexesDeleted[count] = n;
			count++;
		}
	}


    m_nodesList.erase(
    std::remove_if(m_nodesList.begin(), m_nodesList.end(), [this, &toBeDeletedNodes, verboseOutput](const Node& node)
    {
        if(toBeDeletedNodes[&node - &*std::begin(m_nodesList)])
        {
            if(verboseOutput)
            {
                std::cout << "Removing free node (";
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    std::cout << node.m_position[d];
                    if(d == m_dim - 1)
                        std::cout << ")";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
            }

            return true;
        }
        else
           return false;
    }), m_nodesList.end());

    for(std::size_t i = 0 ; i < nodesIndexesDeleted.size() ; ++i)
    {
        for(Element& element : m_elementsList)
        {
            for(std::size_t& n : element.m_nodesIndexes)
            {
                if(n > nodesIndexesDeleted[i])
                    n--;
            }
        }

        for(Node& node : m_nodesList)
        {
            for(std::size_t& n : node.m_neighbourNodes)
            {
                if(n > nodesIndexesDeleted[i])
                    n--;
            }
        }

        for(std::size_t j = i + 1 ; j < nodesIndexesDeleted.size() ; ++j)
            nodesIndexesDeleted[j]--;
    }
}

void Mesh::displayToConsole() const noexcept
{
    std::cout << "Mesh dimension: " << m_dim << "D\n";
    std::cout << "hchar: " << m_hchar << "\n";
    std::cout << "alpha: " << m_alpha << "\n";
    std::cout << "omega: " << m_omega << "\n";
    std::cout << "gamma: " << m_gamma << std::endl;
}

std::vector<std::array<double, 3>> Mesh::getGaussPoints(unsigned int dimension, unsigned int n) const
{
    switch(dimension)
    {
        case 1:
        {
            switch(n)
            {
                case 1:
                    return { {0.0, 0.0, 0.0} };
                case 2:
                    return { {1.0/std::sqrt(3.0), 0, 0},
                             {-1.0/std::sqrt(3.0), 0, 0} };
                case 3:
                    return { {-std::sqrt(3.0/5.0), 0, 0},
                             {0, 0, 0},
                             {std::sqrt(3.0/5.0), 0, 0} };
                case 4:
                    return { {-0.861136311594053, 0, 0},
                             {-0.339981043584856, 0, 0},
                             {0.339981043584856, 0, 0},
                             {0.861136311594053, 0, 0} };
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 1D: " + std::to_string(n));
            }
        }

        case 2:
        {
            switch(n)
            {
                case 1:
                    return { {1.0/3.0, 1.0/3.0, 0.0} };
                case 3:
                    return { {1.0/6.0, 1.0/6.0, 0.0},
                             {1.0/6.0, 2.0/3.0, 0.0},
                             {2.0/3.0, 1.0/6.0, 0.0} };
                case 4:
                    return { {1.0/3.0, 1.0/3.0, 0.0},
                             {0.6, 0.2, 0.0},
                             {0.2, 0.2, 0.0},
                             {0.2, 0.6, 0.0} };
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 2D: " + std::to_string(n));
            }
        }

        case 3:
        {
            switch(n)
            {
                case 1:
                    return { {1.0/4.0, 1.0/4.0, 1.0/4.0} };
                case 4:
                    return{ {0.585410196624968, 0.138196601125011, 0.138196601125011},
                            {0.138196601125011, 0.585410196624968, 0.138196601125011},
                            {0.138196601125011, 0.138196601125011, 0.585410196624968},
                            {0.138196601125011, 0.138196601125011, 0.138196601125011} };
                case 5:
                    return{ {1.0/4.0, 1.0/4.0, 1.0/4.0},
                            {1.0/2.0, 1.0/6.0, 1.0/6.0},
                            {1.0/6.0, 1.0/2.0, 1.0/6.0},
                            {1.0/6.0, 1.0/6.0, 1.0/2.0},
                            {1.0/6.0, 1.0/6.0, 1.0/6.0} };
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 3D: " + std::to_string(n));
            }
        }

        default:
            throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
    }
}

std::vector<double> Mesh::getGaussWeight(unsigned int dimension, unsigned int n) const
{
    switch(dimension)
    {
        case 1:
        {
            switch(n)
            {
                case 1:
                    return {2.0/2.0};
                case 2:
                    return {1.0/2.0, 1.0/2.0};
                case 3:
                    return {5.0/18.0, 8.0/18.0, 5.0/18.0};
                case 4:
                    return {0.347853845137454/2.0,0.652145154862546/2.0 ,0.652145154862546/2.0, 0.347853845137454/2.0};
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 1D: " + std::to_string(n));
            }
        }

        case 2:
        {
            switch(n)
            {
                case 1:
                    return {1.0};
                case 3:
                    return {1.0/3.0, 1.0/3.0, 1.0/3.0};
                case 4:
                    return {-0.5625, 0.520833333333333, 0.520833333333333, 0.520833333333333};
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 2D: " + std::to_string(n));
            }
        }

        case 3:
        {
            switch(n)
            {
                case 1:
                    return {1.0};
                case 4:
                    return {0.25, 0.25, 0.25, 0.25};
                case 5:
                    return {-0.8, 0.45, 0.45, 0.45, 0.45};
                default:
                    throw std::runtime_error("Unexpected number of gauss point in 3D: " + std::to_string(n));
            }
        }

        default:
            throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
    }
}

double Mesh::getRefElementSize(unsigned int dimension) const
{
    switch(dimension)
    {
        case 1:
            return 2.0;

        case 2:
            return 0.5;

        case 3:
            return 0.16666666666666666666666666666667;

        default:
            throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
    }
}

std::vector<std::vector<double>> Mesh::getShapeFunctions(unsigned int dimension, unsigned int n) const
{
    std::vector<std::array<double, 3>> gps = getGaussPoints(dimension, n);

    std::vector<std::vector<double>> sfs(gps.size());

    unsigned int counter = 0;
    for(std::array<double, 3> gp : gps)
    {
        std::vector<double> sf(dimension + 1);

        switch(dimension)
        {
            case 1:
                sf[0] = (1 - gp[0])/2;
                sf[1] = (1 + gp[0])/2;
                break;

            case 2:
                sf[0] = 1 - gp[0] - gp[1];
                sf[1] = gp[0];
                sf[2] = gp[1];
                break;

            case 3:
                sf[0] = 1 - gp[0] - gp[1] - gp[2];
                sf[1] = gp[0];
                sf[2] = gp[1];
                sf[3] = gp[2];
                break;

            default:
                throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
        }

        sfs[counter] = std::move(sf);
        counter++;
    }

    return sfs;
}

std::vector<std::vector<double>> Mesh::getGradShapeFunctions(unsigned int dimension) const
{
    std::vector<std::vector<double>> gradsfs(dimension);

    for(unsigned int i = 0 ; i < gradsfs.size() ; ++i)
        gradsfs[i].resize(dimension + 1);

    switch(dimension)
    {
        case 1:
            gradsfs[0][0] = -1;
            gradsfs[0][1] = 1;
            break;

        case 2:
            gradsfs[0][0] = -1;
            gradsfs[0][1] = 1;
            gradsfs[0][2] = 0;

            gradsfs[1][0] = -1;
            gradsfs[1][1] = 0;
            gradsfs[1][2] = 1;
            break;

        case 3:
            gradsfs[0][0] = -1;
            gradsfs[0][1] = 1;
            gradsfs[0][2] = 0;
            gradsfs[0][3] = 0;

            gradsfs[1][0] = -1;
            gradsfs[1][1] = 0;
            gradsfs[1][2] = 1;
            gradsfs[1][3] = 0;

            gradsfs[1][0] = -1;
            gradsfs[1][1] = 0;
            gradsfs[1][2] = 0;
            gradsfs[1][3] = 1;
            break;

        default:
            throw std::runtime_error("Unexpected dimension: " + std::to_string(dimension));
    }

    return gradsfs;
}

void Mesh::laplacianSmoothingBoundaries()
{
    for(std::size_t elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        const Element& element = m_elementsList[elm];

        if(!element.isContact())
            continue;

        Node* pOutNode = nullptr;
        std::array<Node*, 3> pFacetNodes = {nullptr, nullptr, nullptr};

        unsigned short counter = 0;
        for(std::size_t i = 0 ; i < getNodesPerElm() ; ++i)
        {
            Node& node = m_nodesList[element.getNodeIndex(i)];
            if(!node.isBound())
                pOutNode = &node;
            else
            {
               pFacetNodes[counter] = &node;
               counter++;
            }
        }

        if(counter != m_dim)
            continue;

        double distFromFacet = 0;
        std::array<double, 3> normalF = {0, 0, 0};
        if(m_dim == 2)
        {
            if(pFacetNodes[0]->m_position[0] == pFacetNodes[1]->m_position[0])
            {
                normalF[0] = 1;

                distFromFacet = std::fabs(pFacetNodes[0]->m_position[0] - pOutNode->m_position[0]);
            }
            else if(pFacetNodes[0]->m_position[1] == pFacetNodes[1]->m_position[1])
            {
                normalF[1] = 1;

                distFromFacet = std::fabs(pFacetNodes[0]->m_position[1] - pOutNode->m_position[1]);
            }
            else
            {
                Eigen::Matrix<double, 2, 2> A;
                A << pFacetNodes[0]->m_position[0], pFacetNodes[0]->m_position[1],
                     pFacetNodes[1]->m_position[0], pFacetNodes[1]->m_position[1];

                Eigen::Matrix<double, 2, 1> b;
                b << -1,
                     -1;

                Eigen::Matrix<double, 2, 1> coeff = A.colPivHouseholderQr().solve(b);

                for(unsigned short i = 0 ; i < 2 ; ++i)
                    normalF[i] = coeff[i];

                distFromFacet = std::fabs(normalF[0]*pOutNode->m_position[0] + normalF[1]*pOutNode->m_position[1] + 1)
                          / std::sqrt(normalF[0]*normalF[0] + normalF[1]*normalF[1]);
            }


        }
        else
        {
            if(pFacetNodes[0]->m_position[0] == pFacetNodes[1]->m_position[0] && pFacetNodes[1]->m_position[0] == pFacetNodes[2]->m_position[0])
            {
                normalF[0] = 1;

                distFromFacet = std::fabs(pFacetNodes[0]->m_position[0] - pOutNode->m_position[0]);
            }
            else if(pFacetNodes[0]->m_position[1] == pFacetNodes[1]->m_position[1] && pFacetNodes[1]->m_position[1] == pFacetNodes[2]->m_position[1])
            {
                normalF[1] = 1;

                distFromFacet = std::fabs(pFacetNodes[0]->m_position[1] - pOutNode->m_position[1]);
            }
            else if(pFacetNodes[0]->m_position[2] == pFacetNodes[1]->m_position[2] && pFacetNodes[1]->m_position[1] == pFacetNodes[2]->m_position[2])
            {
                normalF[2] = 1;

                distFromFacet = std::fabs(pFacetNodes[0]->m_position[2] - pOutNode->m_position[2]);
            }
            else
            {
                Eigen::Matrix<double, 3, 3> A;
                A << pFacetNodes[0]->m_position[0], pFacetNodes[0]->m_position[1], pFacetNodes[0]->m_position[2],
                     pFacetNodes[1]->m_position[0], pFacetNodes[1]->m_position[1], pFacetNodes[1]->m_position[2],
                     pFacetNodes[2]->m_position[0], pFacetNodes[2]->m_position[1], pFacetNodes[2]->m_position[2];

                Eigen::Matrix<double, 3, 1> b;
                b << -1,
                     -1,
                     -1;

                Eigen::Matrix<double, 3, 1> coeff = A.colPivHouseholderQr().solve(b);

                for(unsigned short i = 0 ; i < 3 ; ++i)
                    normalF[i] = coeff[i];

                distFromFacet = std::fabs(normalF[0]*pOutNode->m_position[0] + normalF[1]*pOutNode->m_position[1] + normalF[2]*pOutNode->m_position[2] + 1)
                              / std::sqrt(normalF[0]*normalF[0] + normalF[1]*normalF[1] + normalF[2]*normalF[2]);
            }
        }

        if(distFromFacet > m_hchar/3)
            continue;

        if(pOutNode->isOnFreeSurface())
        {
            if(m_dim == 2)
            {
                std::array<double, 3> vecToOutNode = {
                    pOutNode->m_position[0] - pFacetNodes[0]->m_position[0],
                    pOutNode->m_position[1] - pFacetNodes[0]->m_position[1],
                    0
                };

                if(normalF[0]*vecToOutNode[0] + normalF[1]*vecToOutNode[1] < 0)
                {
                    normalF[0] *= -1;
                    normalF[1] *= -1;
                }

                double norm = std::sqrt(normalF[0]*normalF[0] + normalF[1]*normalF[1]);
                normalF[0] /= norm;
                normalF[1] /= norm;

                pOutNode->m_position[0] += m_hchar/3 * normalF[0];
                pOutNode->m_position[1] += m_hchar/3 * normalF[1];
            }
            else
            {
                double norm = std::sqrt(normalF[0]*normalF[0] + normalF[1]*normalF[1] + normalF[2]*normalF[2]);

                normalF[0] /= norm;
                normalF[1] /= norm;
                normalF[2] /= norm;

                std::array<double, 3> vecToOutNode = {
                    pOutNode->getCoordinate(0) - pFacetNodes[0]->getCoordinate(0),
                    pOutNode->getCoordinate(1) - pFacetNodes[0]->getCoordinate(1),
                    pOutNode->getCoordinate(2) - pFacetNodes[0]->getCoordinate(2)
                };

                if(vecToOutNode[0]*normalF[0] + vecToOutNode[1]*normalF[1] + vecToOutNode[2]*normalF[2] < 0)
                {
                    normalF[0] *= -1.0;
                    normalF[1] *= -1.0;
                    normalF[2] *= -1.0;
                }

                pOutNode->m_position[0] += m_hchar/3 * normalF[0];
                pOutNode->m_position[1] += m_hchar/3 * normalF[1];
                pOutNode->m_position[2] += m_hchar/3 * normalF[2];
            }
        }
        else
        {
            std::array<double, 3> newPos = {0, 0, 0};
            for(std::size_t i = 0 ; i < pOutNode->m_neighbourNodes.size() ; ++i)
            {
                const Node& node = m_nodesList[pOutNode->m_neighbourNodes[i]];
                newPos[0] += node.getCoordinate(0);
                newPos[1] += node.getCoordinate(1);
                newPos[2] += node.getCoordinate(2);
            }

            newPos[0] /= static_cast<double>(pOutNode->m_neighbourNodes.size());
            newPos[1] /= static_cast<double>(pOutNode->m_neighbourNodes.size());
            newPos[2] /= static_cast<double>(pOutNode->m_neighbourNodes.size());

            pOutNode->m_position = newPos;
        }
    }
}

void Mesh::loadFromFile(const std::string& fileName)
{
    m_nodesList.clear();

    gmsh::initialize();
#ifndef NDEBUG
    gmsh::option::setNumber("General.Terminal", 1);
#else
    gmsh::option::setNumber("General.Terminal", 0);
#endif

    std::ifstream file(fileName);
    if(file.is_open())
        file.close();
    else
        throw std::runtime_error("the input .msh file does not exist!");

    gmsh::open(fileName);

    // Check that the mesh is not 3D
    computeMeshDim();

    if(m_boundingBox.size() != 2*m_dim)
        throw std::runtime_error("Invalid bounding box size: " + std::to_string(m_boundingBox.size()));

    for(auto& exclusionZone : m_exclusionZones)
    {
        if(exclusionZone.size() != 2*m_dim)
            throw std::runtime_error("Invalid exclusion zone size: " + std::to_string(exclusionZone.size()));
    }

    // We retrieve the tags of the physical groups of dimension m_dim and
    // m_dim-1
    std::vector<std::pair<int, int>> physGroupHandlesHD;
    gmsh::model::getPhysicalGroups(physGroupHandlesHD, m_dim);

    std::vector<std::pair<int, int>> physGroupHandlesLD;
    gmsh::model::getPhysicalGroups(physGroupHandlesLD, m_dim - 1);

    //The file should contain physical group for the boundary, the free surface and the fluid

    for(auto physGroupLD : physGroupHandlesLD)
    {
        std::string name;
        gmsh::model::getPhysicalName(m_dim - 1, physGroupLD.second, name);
        if(name != "FreeSurface" && name != "Fluid")
        {
            std::vector<double> coord;
            std::vector<std::size_t> dummyNodesTagsBoundary;
            gmsh::model::mesh::getNodesForPhysicalGroup(m_dim - 1, physGroupLD.second,
                                                        dummyNodesTagsBoundary, coord);

            for(std::size_t i = 0 ; i < dummyNodesTagsBoundary.size() ; ++i)
            {
                Node node(*this);
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    node.m_position[d] = coord[3*i + d];
                }

                //If the nodes is already on the boundary, we do not add it twice
                if(std::find(m_nodesList.begin(), m_nodesList.end(), node) != m_nodesList.end())
                    continue;

                node.m_isBound = true;

                auto posBCinTagNames = std::find(m_tagNames.begin(), m_tagNames.end(), name);
                if(posBCinTagNames == std::end(m_tagNames))
                {
                    m_tagNames.push_back(name);
                    node.m_tag = static_cast<int>(m_tagNames.size() - 1);
                }
                else
                {
                    node.m_tag = static_cast<int>(std::distance(m_tagNames.begin(), posBCinTagNames));
                }
                m_nodesList.push_back(std::move(node));
                m_boundaryInitialPos[m_nodesList.size() - 1] = m_nodesList.back().m_position;
            }
        }
    }

    for(auto physGroupHD : physGroupHandlesHD)
    {
        std::string name;
        gmsh::model::getPhysicalName(m_dim, physGroupHD.second, name);

        std::vector<std::size_t> dummyNodesTags;
        std::vector<double> coord;
        gmsh::model::mesh::getNodesForPhysicalGroup(m_dim, physGroupHD.second,
                                                    dummyNodesTags, coord);

        for(std::size_t i = 0 ; i < dummyNodesTags.size() ; ++i)
        {
            Node node(*this);
            for(unsigned short d = 0 ; d < m_dim ; ++d)
                node.m_position[d] = coord[3*i + d];

            //If the nodes is already on the boundary, we do not add it twice
            if(std::find(m_nodesList.begin(), m_nodesList.end(), node) != m_nodesList.end())
                continue;

            node.m_isBound = false;

            auto posBCinTagNames = std::find(m_tagNames.begin(), m_tagNames.end(), name);
            if(posBCinTagNames == std::end(m_tagNames))
            {
                m_tagNames.push_back(name);
                node.m_tag = m_tagNames.size() - 1;
            }
            else
            {
                node.m_tag = static_cast<uint16_t>(std::distance(m_tagNames.begin(), posBCinTagNames));
            }

            m_nodesList.push_back(std::move(node));
        }
    }

#ifndef NDEBUG
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        for(std::size_t n2 = 0 ; n2 < m_nodesList.size() ; ++n2)
        {
            if(n != n2)
            {
                if(m_nodesList[n] == m_nodesList[n2])
                {
                    std::cerr << "Duplicate node found: " << "(";
                    for(unsigned short d = 0 ; d < m_dim ; ++d)
                    {
                        std::cerr << m_nodesList[n].m_position[d];
                        if(d == m_dim - 1)
                            std::cout << ")";
                        else
                            std::cout << ", ";
                    }
                    std::cerr << std::endl;
                    throw std::runtime_error("duplicates nodes!");
                }
            }
        }
    }
    std::cout << "List of founded BC names: " << std::endl;
    for (auto& aString : m_tagNames)
    {
        std::cout << aString << std::endl;
    }
#endif

    if(m_nodesList.empty())
        throw std::runtime_error("no nodes loaded! Did you add the physical groups in the .geo file?");

    gmsh::finalize();

    triangulateAlphaShape();
}

void Mesh::remesh(bool verboseOutput)
{
    //laplacianSmoothingBoundaries();
    addNodes(verboseOutput);
    removeNodes(verboseOutput);
    checkBoundingBox(verboseOutput);
    triangulateAlphaShape();
}

bool Mesh::removeNodes(bool verboseOutput) noexcept
{
    assert(!m_elementsList.empty() && !m_nodesList.empty() && "There is no mesh !");

    std::vector<bool> touched(m_nodesList.size(), false);       //Is the node next to a node which should be deleted
    std::vector<bool> toBeDeleted(m_nodesList.size(), false);   //Should the node be deleted

    double limitLength = m_gamma*m_hchar;

    bool removeNodes = false;

    for(std::size_t i = 0 ; i < m_nodesList.size() ; ++i)
    {
        //We the node is free, it is not in an element; if a node is already tagged,
        //we can skip it (do not delete node twice)
        if(toBeDeleted[i] || m_nodesList[i].isFree())
            continue;

        for(unsigned int j = 0 ; j < m_nodesList[i].m_neighbourNodes.size() ; ++j)
        {
            double d = 0;

            for(unsigned short k = 0 ; k < m_dim ; ++k)
            {
                d +=(m_nodesList[i].m_position[k] - m_nodesList[m_nodesList[i].m_neighbourNodes[j]].m_position[k])
                   *(m_nodesList[i].m_position[k] - m_nodesList[m_nodesList[i].m_neighbourNodes[j]].m_position[k]);
            }
            d = std::sqrt(d);

            //Two nodes are too close.
            if(d <= limitLength)
            {
                //If the neighbour nodes is touched, we delete the current nodes
                if(touched[m_nodesList[i].m_neighbourNodes[j]])
                {
                    //Do not delete bounded or free surface nodes
                    if(m_nodesList[i].m_isBound || m_nodesList[i].m_isOnFreeSurface)
                        continue;

                    toBeDeleted[i] = true;
                    removeNodes = true;
                }
                //If the neighbour nodes is not touched, we delete the neoghbour nodes
                else
                {
                    //Do not delete bounded or free surface nodes
                    if(m_nodesList[m_nodesList[i].m_neighbourNodes[j]].m_isBound ||
                       m_nodesList[m_nodesList[i].m_neighbourNodes[j]].m_isOnFreeSurface)
                        continue;

                    touched[i] = true;
                    toBeDeleted[m_nodesList[i].m_neighbourNodes[j]] = true;
                    removeNodes = true;
                }
            }
        }
    }

    m_nodesList.erase(
    std::remove_if(m_nodesList.begin(), m_nodesList.end(), [this, &toBeDeleted, verboseOutput](const Node& node)
    {
       if(toBeDeleted[&node - &*std::begin(m_nodesList)])
       {
           if(verboseOutput)
           {
                std::cout << "Removing node " << "(";
                for(unsigned short d = 0 ; d < m_dim ; ++d)
                {
                    std::cout << node.m_position[d];
                    if(d == m_dim - 1)
                        std::cout << ")";
                    else
                        std::cout << ", ";
                }
                std::cout << std::endl;
           }

           return true;
       }
       else
           return false;
    }), m_nodesList.end());

    return removeNodes;
}

void Mesh::restoreNodesList()
{
    if(m_nodesListSave.empty())
        throw std::runtime_error("the nodes list was not saved before or does not exist!");

    m_nodesList = std::move(m_nodesListSave);

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        m_elementsList[elm].computeJ();
        m_elementsList[elm].computeDetJ();
        m_elementsList[elm].computeInvJ();
    }

    #pragma omp parallel for default(shared)
    for(std::size_t facet = 0 ; facet < m_facetsList.size() ; ++facet)
    {
        m_facetsList[facet].computeJ();
        m_facetsList[facet].computeDetJ();
        m_facetsList[facet].computeInvJ();
    }

    computeFSNormalCurvature();
}

void Mesh::saveNodesList()
{
    if(m_nodesList.empty())
        throw std::runtime_error("the nodes list does not exist!");

    m_nodesListSave = m_nodesList;
}

void Mesh::triangulateAlphaShape()
{
    if(m_dim == 2)
        triangulateAlphaShape2D();
    else
        triangulateAlphaShape3D();
}

void Mesh::updateNodesPosition(const std::vector<double>& deltaPos)
{
    if(deltaPos.size() != m_nodesList.size()*m_dim)
        throw std::runtime_error("invalid size of the deltaPos vector: " + std::to_string(deltaPos.size()) + " VS " + std::to_string(m_nodesListSave.size()*m_dim));

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(!m_nodesList[n].m_isFixed)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].m_position[d] += deltaPos[n + d*m_nodesList.size()];
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        m_elementsList[elm].computeJ();
        m_elementsList[elm].computeDetJ();
        m_elementsList[elm].computeInvJ();
    }

    #pragma omp parallel for default(shared)
    for(std::size_t facet = 0 ; facet < m_facetsList.size() ; ++facet)
    {
        m_facetsList[facet].computeJ();
        m_facetsList[facet].computeDetJ();
        m_facetsList[facet].computeInvJ();
    }

    computeFSNormalCurvature();
}

void Mesh::updateNodesPosition(const Eigen::VectorXd& deltaPos)
{
    if(static_cast<std::size_t>(deltaPos.rows()) < m_nodesList.size()*m_dim)
        throw std::runtime_error("invalid size of the deltaPos vector: " + std::to_string(deltaPos.rows()) + " VS " + std::to_string(m_nodesListSave.size()*m_dim));

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(!m_nodesList[n].m_isFixed)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].m_position[d] += deltaPos[n + d*m_nodesList.size()];
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        m_elementsList[elm].computeJ();
        m_elementsList[elm].computeDetJ();
        m_elementsList[elm].computeInvJ();
    }

    #pragma omp parallel for default(shared)
    for(std::size_t facet = 0 ; facet < m_facetsList.size() ; ++facet)
    {
        m_facetsList[facet].computeJ();
        m_facetsList[facet].computeDetJ();
        m_facetsList[facet].computeInvJ();
    }

    computeFSNormalCurvature();
}

void Mesh::updateNodesPosition(const Eigen::VectorXd& deltaPos,
                               const std::vector<std::size_t> nodesIndexes)
{
    if(static_cast<std::size_t>(deltaPos.rows()/m_dim) != nodesIndexes.size())
        throw std::runtime_error("size of deltaPos and nodesIndexes vector must be dim*x and x!");

    std::set<std::size_t> elementModified;
    std::set<std::size_t> facetModified;
    for(std::size_t n : nodesIndexes)
    {
        const Node& node = m_nodesList[n];
        for(unsigned int e = 0 ; e < node.getElementCount() ; ++e)
        {
            elementModified.insert(node.getElementMeshIndex(e));

        }

        for(unsigned int f = 0 ; f < node.getFacetCount() ; ++f)
        {
            facetModified.insert(node.getFacetMeshIndex(f));

        }
    }

    #pragma omp parallel for default(shared)
    for(auto i = 0 ; i < deltaPos.rows()/m_dim ; ++i)
    {
        std::size_t n = nodesIndexes[i];

        if(!m_nodesList[n].m_isFixed)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].m_position[d] += deltaPos[i + d*nodesIndexes.size()];
            }
        }
    }

    for(auto it = elementModified.begin() ; it != elementModified.end() ; ++it)
    {
        m_elementsList[*it].computeJ();
        m_elementsList[*it].computeDetJ();
        m_elementsList[*it].computeInvJ();
    }

    for(auto it = facetModified.begin() ; it != facetModified.end() ; ++it)
    {
        m_facetsList[*it].computeJ();
        m_facetsList[*it].computeDetJ();
        m_facetsList[*it].computeInvJ();
    }

    //To do: only update normal of facet that have been updated
    computeFSNormalCurvature();
}

void Mesh::updateNodesPositionFromSave(const std::vector<double>& deltaPos)
{
    if(m_nodesListSave.empty())
        throw std::runtime_error("you did not save the nodes list!");

    if(deltaPos.size() != m_nodesListSave.size()*m_dim)
        throw std::runtime_error("invalid size of the deltaPos vector: " + std::to_string(deltaPos.size()) + " VS " + std::to_string(m_nodesListSave.size()*m_dim));

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(!m_nodesList[n].m_isFixed)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].m_position[d] = m_nodesListSave[n].m_position[d] + deltaPos[n + d*m_nodesList.size()];
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        m_elementsList[elm].computeJ();
        m_elementsList[elm].computeDetJ();
        m_elementsList[elm].computeInvJ();
    }

    #pragma omp parallel for default(shared)
    for(std::size_t facet = 0 ; facet < m_facetsList.size() ; ++facet)
    {
        m_facetsList[facet].computeJ();
        m_facetsList[facet].computeDetJ();
        m_facetsList[facet].computeInvJ();
    }

    computeFSNormalCurvature();
}

void Mesh::updateNodesPositionFromSave(const Eigen::VectorXd& deltaPos)
{
    if(m_nodesListSave.empty())
        throw std::runtime_error("you did not save the nodes list!");

    if(static_cast<std::size_t>(deltaPos.rows()) < m_nodesListSave.size()*m_dim)
        throw std::runtime_error("invalid size of the deltaPos vector: " + std::to_string(deltaPos.rows()) + " VS " + std::to_string(m_nodesListSave.size()*m_dim));

    #pragma omp parallel for default(shared)
    for(std::size_t n = 0 ; n < m_nodesList.size() ; ++n)
    {
        if(!m_nodesList[n].m_isFixed)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].m_position[d] = m_nodesListSave[n].m_position[d] + deltaPos[n + d*m_nodesList.size()];
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < m_elementsList.size() ; ++elm)
    {
        m_elementsList[elm].computeJ();
        m_elementsList[elm].computeDetJ();
        m_elementsList[elm].computeInvJ();
    }

    #pragma omp parallel for default(shared)
    for(std::size_t facet = 0 ; facet < m_facetsList.size() ; ++facet)
    {
        m_facetsList[facet].computeJ();
        m_facetsList[facet].computeDetJ();
        m_facetsList[facet].computeInvJ();
    }

    computeFSNormalCurvature();
}

void Mesh::updateNodesPositionFromSave(const Eigen::VectorXd& deltaPos,
                                       const std::vector<std::size_t> nodesIndexes)
{
    if(m_nodesListSave.empty())
        throw std::runtime_error("you did not save the nodes list!");

    if(static_cast<std::size_t>(deltaPos.rows()/m_dim) != nodesIndexes.size())
        throw std::runtime_error("size of deltaPos and nodesIndexes vector must be dim*x and x!");

    std::set<std::size_t> elementModified;
    std::set<std::size_t> facetModified;
    for(std::size_t n : nodesIndexes)
    {
        const Node& node = m_nodesList[n];
        for(unsigned int e = 0 ; e < node.getElementCount() ; ++e)
        {
            elementModified.insert(node.getElementMeshIndex(e));

        }

        for(unsigned int f = 0 ; f < node.getFacetCount() ; ++f)
        {
            facetModified.insert(node.getFacetMeshIndex(f));

        }
    }

    #pragma omp parallel for default(shared)
    for(auto i = 0 ; i < deltaPos.rows()/m_dim ; ++i)
    {
        std::size_t n = nodesIndexes[i];

        if(!m_nodesList[n].m_isFixed)
        {
            for(unsigned short d = 0 ; d < m_dim ; ++d)
            {
                m_nodesList[n].m_position[d] = m_nodesListSave[n].m_position[d] + deltaPos[i + d*nodesIndexes.size()];
            }
        }
    }

    for(auto it = elementModified.begin() ; it != elementModified.end() ; ++it)
    {
        m_elementsList[*it].computeJ();
        m_elementsList[*it].computeDetJ();
        m_elementsList[*it].computeInvJ();
    }

    for(auto it = facetModified.begin() ; it != facetModified.end() ; ++it)
    {
        m_facetsList[*it].computeJ();
        m_facetsList[*it].computeDetJ();
        m_facetsList[*it].computeInvJ();
    }

    //To do: only update normal of facet that have been updated
    computeFSNormalCurvature();
}
