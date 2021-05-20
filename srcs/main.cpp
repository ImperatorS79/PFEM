#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>

#include "simulation/physics/Problems.hpp"
#include "simulation/utility/Utility.hpp"
#include "simulation/utility/SolTable.hpp"
#include "simulation/utility/Clock.hpp"


/**
 * \param  argv[1] .lua file that contains the parameters.
 */
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr   << "Usage: " << argv[0] << " params.lua" <<  std::endl;
        return 1;
    }

    Clock myClock;
    myClock.start();

    std::unique_ptr<Problem> pProblem;

    try
    {
        std::ifstream luaFile(argv[1]);
        if(!luaFile.is_open())
            throw std::runtime_error("cannot open lua file " + std::string(argv[1]) + "!");

        luaFile.close();

        sol::state state;
        auto res = state.safe_script_file(argv[1]);
        if(!res.valid())
        {
            sol::error err = res;
            throw std::runtime_error(std::string("an error occured while loading lua file: ") + err.what());
        }
        state = {};
        state.script_file(argv[1]);

        SolTable table = SolTable("Problem", state);

        std::string problemType = table.checkAndGet<std::string>("id");

        if(problemType == "IncompNewtonNoT" || problemType == "Bingham" || problemType == "Boussinesq" || problemType == "Conduction")
            pProblem = std::make_unique<ProbIncompNewton>(argv[1]);
        else if(problemType == "WCompNewtonNoT" || problemType == "BoussinesqWC")
            pProblem = std::make_unique<ProbWCompNewton>(argv[1]);
        else
            throw std::runtime_error("unknown problem type " + problemType + "!");

        pProblem->simulate();
        std::cout << "======================================" << std::endl;
        std::cout << "======================================" << std::endl;
        pProblem->displayTimeStats();
    }
    catch(const std::exception& e)
    {
        std::cerr << std::endl << "\nSomething went wrong while running the program: " << e.what() << std::endl;
        pProblem->dump();
        return -1;
    }
    catch(...)
    {
        std::cerr << std::endl << "\nAn unknown exception has occurred." << std::endl;
        pProblem->dump();
        return -1;
    }

    myClock.end();
    std::cout << "======================================" << std::endl;
    myClock.displayDT("Ellapsed time for simulation: ");

    return 0;
}
