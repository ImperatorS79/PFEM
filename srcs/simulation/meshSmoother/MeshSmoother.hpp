#ifndef MESHSMOOTHER_HPP_INCLUDED
#define MESHSMOOTHER_HPP_INCLUDED

#include "../../mesh/Mesh.hpp"
#include "../Problem.hpp"

class MeshSmoother
{
    public:
        MeshSmoother(Problem* pProblem, Mesh& mesh):
        m_pProblem(pProblem),
        m_mesh(mesh)
        {

        }

        virtual ~MeshSmoother()
        {

        }

        virtual void smooth(bool /** verboseOuput **/)
        {

        }

    protected:
        Problem* m_pProblem;
        Mesh& m_mesh;
};

#endif // MESHSMOOTHER_HPP_INCLUDED
