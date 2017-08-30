
#pragma once

#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL
{

class MediumFactory
{
    public:
    void Register(const std::string& name, Medium* (*)(void));
    Medium* CreateMedium(const std::string&);

    static MediumFactory& Get()
    {
        static MediumFactory instance;
        return instance;
    }

    private:
    MediumFactory();
    ~MediumFactory() { medium_map.clear(); }

    std::map<std::string, Medium* (*)(void)> medium_map;
};

} /* PROPOsAL */

