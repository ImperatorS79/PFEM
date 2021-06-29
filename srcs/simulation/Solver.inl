#include "Solver.hpp"

inline bool Solver::getBcTagFlags(int tag, unsigned short flag) const noexcept
{
    auto it = m_bcTagFlags.find(tag);

    if(it == m_bcTagFlags.end())
        return false;
    else
        return it->second[flag];
}

inline double Solver::getTimeStep() const noexcept
{
    return m_timeStep;
}
