#pragma once
#ifndef GETME_HPP_INCLUDED
#define GETME_HPP_INCLUDED

#include "MeshSmoother.hpp"

template<unsigned short dim>
class GETMe: public MeshSmoother
{
    public:
        GETMe(Problem* pProblem, Mesh& mesh, unsigned int maxIter, double epsToll);
        ~GETMe() override;

        void smooth(bool verboseOuput) override;
    private:
        unsigned int m_maxIter;
        double m_epsToll;

        std::vector<double> m_rins; //1 per elm


};

#include "GETMe.inl"

#endif // GETME_HPP_INCLUDED
