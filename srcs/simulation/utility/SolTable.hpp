#pragma once
#ifndef SOLTABLE_HPP_INCLUDED
#define SOLTABLE_HPP_INCLUDED

#include <string>

#define SOL_ALL_SAFETIES_ON 1
#include <sol/sol.hpp>

class SolTable
{
    public:
        SolTable() = default;

        SolTable(const std::string& tableName, const sol::state& state):
        m_tableName(tableName)
        {
            auto obj = state[tableName];
            if(!obj.valid())
                throw std::runtime_error("Could not find table" + tableName + "inside sol::state!");

            m_tableInternal = state[tableName];
        }

        SolTable(const std::string& tableName, const SolTable& solTable):
        m_tableName(tableName)
        {
            auto obj = solTable.m_tableInternal[tableName];
            if(!obj.valid())
                throw std::runtime_error("Could not find table " + tableName + " inside table " + solTable.getName());

            m_tableInternal = solTable.m_tableInternal[tableName];
        }

        SolTable(sol::table table):
        m_tableName("anonymous")
        {
            m_tableInternal = table;
        }

        ~SolTable() = default;

        template<typename T, typename... Args>
        T call(const std::string& functionName, Args... args) const noexcept
        {
            sol::unsafe_function_result res = m_tableInternal[functionName](m_tableInternal, args...);

            return res.get<T>();
        }

        template<typename... Args>
        void checkCall(const std::string& functionName, Args... args) const
        {
            auto object = m_tableInternal[functionName](m_tableInternal, args...);
            if(!object.valid())
            {
                sol::error err = object;
                throw std::runtime_error("Cannot find method " + functionName + " inside table " + m_tableName + ": " +
                      std::string(err.what()) + "!");
            }
        }

        template<typename... Args>
        bool checkCallNoThrow(const std::string& functionName, Args... args) const
        {
            auto object = m_tableInternal[functionName](m_tableInternal, args...);
            return object.valid();
        }

        template<typename T>
        T get(const std::string& propertyName) const noexcept
        {
            auto object = m_tableInternal[propertyName];
            return object.get<T>();
        }

        template<typename T>
        T checkAndGet(const std::string& propertyName) const
        {
            auto obj = m_tableInternal[propertyName];
            if(!obj.valid())
            {
                throw std::runtime_error("Property " + propertyName + " was not found in table " + m_tableName + "!");
            }

            return obj.get<T>();
        }

        bool doesVarExist(const std::string& propertyName) const
        {
            auto obj = m_tableInternal[propertyName];
            return obj.valid();
        }

        void for_each(std::function<void(sol::object /*key*/, sol::object /*value*/)> f) const
        {
            return m_tableInternal.for_each(f);
        }

        inline std::string getName() const noexcept
        {
            return m_tableName;
        }

    private:
        sol::table m_tableInternal;
        std::string m_tableName;
};

#endif // SOLTABLE_HPP_INCLUDED
