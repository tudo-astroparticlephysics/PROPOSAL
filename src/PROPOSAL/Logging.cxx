#include "PROPOSAL/Logging.h"
#include <spdlog/sinks/stdout_color_sinks.h>

using namespace PROPOSAL;

unordered_map<string, unique_ptr<spdlog::logger>> Logging::logger = {};

shared_ptr<spdlog::sinks::sink> Logging::sink
    = make_shared<spdlog::sinks::stdout_color_sink_mt>();

spdlog::level::level_enum Logging::global_loglevel = spdlog::level::level_enum::warn;