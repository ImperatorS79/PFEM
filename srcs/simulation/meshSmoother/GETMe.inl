#include "GETMe.hpp"

template<unsigned short dim>
GETMe<dim>::GETMe(Problem* pProblem, Mesh& mesh, double timeBetweenSmooth, unsigned int maxIter, double epsToll) :
MeshSmoother(pProblem, mesh, timeBetweenSmooth),
m_maxIter(maxIter),
m_epsToll(epsToll)
{

}

template<unsigned short dim>
GETMe<dim>::~GETMe()
{

}

template<unsigned short dim>
void GETMe<dim>::smooth(bool verboseOutput)
{
    const std::size_t nElm = m_mesh.getElementsCount();

    m_rins.resize(nElm);
    #pragma omp parallel for default(shared)
    for(std::size_t elm = 0 ; elm < nElm ; ++elm)
    {
        const Element& element = m_mesh.getElement(elm);
        m_rins[elm] = element.getRin();
    }

    double res = std::numeric_limits<double>::max();

    //TO DO: store all the normal computation function in mesh to avoid copy paste everywhere
    auto computeNormal2D = [](const Element& element, unsigned short n1, unsigned short n2, unsigned short oppN) -> std::array<double, 3>
    {
        const Node& node1 = element.getNode(n1);
        const Node& node2 = element.getNode(n2);
        const Node& oppNode = element.getNode(oppN);

        double x1 = node1.getCoordinate(0);
        double y1 = node1.getCoordinate(1);
        double x2 = node2.getCoordinate(0);
        double y2 = node2.getCoordinate(1);
        double xOpp = oppNode.getCoordinate(0);
        double yOpp = oppNode.getCoordinate(1);

        std::array<double, 3> normal = {
            y2 - y1,
            -x1 + x2,
            0
        };

        std::array<double, 3> vecToOppNode = {
            xOpp - x1,
            yOpp - y1,
            0
        };

        if(vecToOppNode[0]*normal[0] + vecToOppNode[1]*normal[1] < 0)
        {
            normal[0] *= -1;
            normal[1] *= -1;
        }

        double norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

        normal[0] /= norm;
        normal[1] /= norm;

        return normal;
    };

    auto computeNormal3D = [](const Element& element, unsigned short n1, unsigned short n2, unsigned short n3, unsigned short oppN) -> std::array<double, 3>
    {
        const Node& node1 = element.getNode(n1);
        const Node& node2 = element.getNode(n2);
        const Node& node3 = element.getNode(n3);
        const Node& oppNode = element.getNode(oppN);

        double x1 = node1.getCoordinate(0);
        double y1 = node1.getCoordinate(1);
        double z1 = node1.getCoordinate(2);
        double x2 = node2.getCoordinate(0);
        double y2 = node2.getCoordinate(1);
        double z2 = node2.getCoordinate(2);
        double x3 = node3.getCoordinate(0);
        double y3 = node3.getCoordinate(1);
        double z3 = node3.getCoordinate(2);
        double xOpp = oppNode.getCoordinate(0);
        double yOpp = oppNode.getCoordinate(1);
        double zOpp = oppNode.getCoordinate(2);

        std::array<double, 3> normal = {
            (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1),
            (z2 - z1)*(x3 - x1) - (z3 - z1)*(x2 - x1),
            (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y3)
        };

        std::array<double, 3> vecToOppNode = {
            xOpp - x1,
            yOpp - y1,
            zOpp - z1
        };

        if(vecToOppNode[0]*normal[0] + vecToOppNode[1]*normal[1] + vecToOppNode[2]*normal[2] < 0)
        {
            normal[0] *= -1;
            normal[1] *= -1;
            normal[2] *= -1;
        }

        double norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

        normal[0] /= norm;
        normal[1] /= norm;
        normal[2] /= norm;

        return normal;
    };

    std::size_t counter = 0;
    while(res > m_epsToll)
    {
        std::size_t elmToSmooth = 0;
        double minRin = std::numeric_limits<double>::max();
        bool foundOne = false;
        for(std::size_t elm = 0 ; elm < nElm ; ++elm)
        {
            const Element& element = m_mesh.getElement(elm);

            bool keep = false;
            for(unsigned short n = 0 ; n < dim + 1 ; ++n)
            {
                if(!element.getNode(n).isBound())
                {
                    keep = true;
                    foundOne = true;
                    break;
                }
            }

            if(!keep)
                continue;

            if(m_rins[elm] < minRin)
            {
                minRin = m_rins[elm];
                elmToSmooth = elm;
            }
        }

        if(!foundOne)
        {
            if(verboseOutput)
                std::cout << "No element found to smooth" << std::endl;
            return;
        }

        if(verboseOutput)
            std::cout << "Element with smallest rin: " << elmToSmooth << ", rin: " << minRin << std::endl;

        const Element& elementToSmooth = m_mesh.getElement(elmToSmooth);
        double sigmaE = 0;

        std::array<std::array<double, 3>, dim + 1> opposedNodeNormals; //normal to face opposed to node oriented to element exterior


        if constexpr(dim == 2)
        {
            opposedNodeNormals[0] = computeNormal2D(elementToSmooth, 1, 2, 0);
            opposedNodeNormals[1] = computeNormal2D(elementToSmooth, 0, 2, 1);
            opposedNodeNormals[2] = computeNormal2D(elementToSmooth, 0, 1, 2);
            sigmaE = (Node::distance(elementToSmooth.getNode(0), elementToSmooth.getNode(1))
                    + Node::distance(elementToSmooth.getNode(2), elementToSmooth.getNode(1))
                    + Node::distance(elementToSmooth.getNode(0), elementToSmooth.getNode(2)))/30;
        }
        else if constexpr(dim == 3)
        {
            opposedNodeNormals[0] = computeNormal3D(elementToSmooth, 1, 2, 3, 0);
            opposedNodeNormals[1] = computeNormal3D(elementToSmooth, 0, 2, 3, 1);
            opposedNodeNormals[2] = computeNormal3D(elementToSmooth, 0, 1, 3, 2);
            opposedNodeNormals[3] = computeNormal3D(elementToSmooth, 0, 1, 2, 3);
            sigmaE = (Node::distance(elementToSmooth.getNode(0), elementToSmooth.getNode(1))
                    + Node::distance(elementToSmooth.getNode(2), elementToSmooth.getNode(1))
                    + Node::distance(elementToSmooth.getNode(2), elementToSmooth.getNode(0))
                    + Node::distance(elementToSmooth.getNode(0), elementToSmooth.getNode(3))
                    + Node::distance(elementToSmooth.getNode(2), elementToSmooth.getNode(3))
                    + Node::distance(elementToSmooth.getNode(1), elementToSmooth.getNode(3)))/60;
        }

        Eigen::VectorXd deltaPos(dim*(dim + 1));
        std::vector<std::size_t> nodeIndexProbToMesh(dim + 1);

        bool updateSuccessFull = false;
        double fact = 1;
        double newRin = 0;
        double innerCounter = 0;
        while(!updateSuccessFull)
        {
            std::cout << "sigma: " << fact*sigmaE << std::endl;
            for(unsigned short n = 0 ; n < dim + 1 ; ++n)
            {
                const Node& node = elementToSmooth.getNode(n);

                double actualSigma = fact*sigmaE;
                if(node.isBound())
                    actualSigma = 0;
                else if(node.isOnFreeSurface())
                    actualSigma/=4;

                if constexpr(dim == 2)
                {
                    deltaPos[n] = actualSigma*opposedNodeNormals[n][0];
                    deltaPos[n + (dim + 1)] = actualSigma*opposedNodeNormals[n][1];
                }
                else if constexpr(dim == 3)
                {
                    deltaPos[n] = actualSigma*opposedNodeNormals[n][0];
                    deltaPos[n + (dim + 1)] = actualSigma*opposedNodeNormals[n][1];
                    deltaPos[n + 2*(dim + 1)] = actualSigma*opposedNodeNormals[n][2];
                }
                nodeIndexProbToMesh[n] = elementToSmooth.getNodeIndex(n);
            }

            m_mesh.updateNodesPosition(deltaPos, nodeIndexProbToMesh);

            newRin = elementToSmooth.getRin();

            updateSuccessFull = true;

            std::vector<double> newRins(elementToSmooth.getNeighbourElementsCount());
            for(std::size_t elm = 0 ; elm < elementToSmooth.getNeighbourElementsCount() ; ++elm)
            {
                const Element& element = m_mesh.getElement(elementToSmooth.getNeighbourElmIndex(elm));
                std::size_t elmIndex = elementToSmooth.getNeighbourElmIndex(elm);
                double rin = element.getRin();

                std::cout << elmIndex << ": " << rin << std::endl;

                bool keep = false;
                for(unsigned short n = 0 ; n < dim + 1 ; ++n)
                {
                    if(!element.getNode(n).isBound())
                    {
                        keep = true;
                        break;
                    }
                }

                if(!keep)
                    continue;

                if(rin > minRin)
                    newRins[elm] = rin;
                else
                {
                    updateSuccessFull = false;
                    break;
                }
            }

            if(!updateSuccessFull)
            {
                m_mesh.updateNodesPosition(-deltaPos, nodeIndexProbToMesh);
                fact *= 0.75;

                if(verboseOutput)
                    std::cout << "GETMe produced inverted elements or smaller rin, decreasing sigma" << std::endl;
            }
            else
            {
                m_rins[elmToSmooth] = newRin;
                for(std::size_t elm = 0 ; elm < elementToSmooth.getNeighbourElementsCount() ; ++elm)
                {
                    const Element& element = m_mesh.getElement(elementToSmooth.getNeighbourElmIndex(elm));
                    std::size_t elmIndex = element.getNeighbourElmIndex(elm);

                    m_rins[elmIndex] = newRins[elm];
                }
            }

            innerCounter++;

            if(innerCounter > m_maxIter)
            {
                if(verboseOutput)
                    std::cout << "Stopping GETMe, inner iterations greater than: " << m_maxIter << std::endl;

                return;
            }
        }

        if(newRin < 0)
            std::cerr << "dude, wtf ?" << std::endl;

        if(counter == 0)
            res = std::numeric_limits<double>::max();
        else
            res = (elementToSmooth.getRin() - minRin)/minRin;

        if(verboseOutput)
            std::cout << "GETMe residual: " << res << " vs " << m_epsToll << std::endl;

        counter++;

        if(counter > m_maxIter)
        {
            if(verboseOutput)
                std::cout << "Stopping GETMe, iterations greater than: " << m_maxIter << std::endl;

            return;
        }
    }
}
