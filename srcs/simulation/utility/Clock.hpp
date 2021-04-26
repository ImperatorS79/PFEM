#pragma once
#ifndef CLOCK_HPP_INCLUDED
#define CLOCK_HPP_INCLUDED

#include <chrono>
#include <iostream>

typedef std::conditional<
       std::chrono::high_resolution_clock::is_steady,
       std::chrono::high_resolution_clock,
       std::chrono::steady_clock >::type ClockType;

class Clock
{
    public:
        Clock()
        {
            static_assert( ClockType::is_steady, "Clock is not monotonically-increasing (steady).");
        }
        ~Clock() = default;

        void start()
        {
            m_start = ClockType::now();
        }

        double end()
        {
            m_end = ClockType::now();
            auto ellapsedTimeMeasure = std::chrono::duration_cast<std::chrono::nanoseconds>(m_end - m_start);
            return static_cast<double>(ellapsedTimeMeasure.count())/1000000000.0;
        }

        void displayDT(const std::string& text)
        {
            auto ellapsedTimeMeasure = std::chrono::duration_cast<std::chrono::nanoseconds>(m_end - m_start);
            std::cout << text << static_cast<double>(ellapsedTimeMeasure.count())/1000000000.0 << " s" << std::endl;
        }

    private:
        std::chrono::time_point<ClockType> m_start;
        std::chrono::time_point<ClockType> m_end;
};

#endif // CLOCK_HPP_INCLUDED

