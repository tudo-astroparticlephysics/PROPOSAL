#include "PROPOSAL/Output.h"


void Output::SetLoggingConfigurationFile(std::string file)
{
    PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT(file));
}

