#include "PROPOSAL/Logging.h"
#include <spdlog/sinks/stdout_color_sinks.h>

using namespace PROPOSAL;

std::unordered_map<std::string, std::shared_ptr<spdlog::logger>> Logging::logger = {};

std::shared_ptr<spdlog::sinks::sink> Logging::sink
    = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

spdlog::level::level_enum Logging::global_loglevel = spdlog::level::level_enum::warn;
