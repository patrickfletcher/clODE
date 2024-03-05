//helper functions for MEX interface to clODE

#ifndef CLODE_MEX_LOGGING_H
#define CLODE_MEX_LOGGING_H

#include "spdlog/spdlog.h"
#include "logging/MatlabSink.hpp"

struct LoggerSingleton {
        std::shared_ptr<MatlabSink_mt> sink;
        std::shared_ptr<spdlog::logger> matlab_logger;
        LoggerSingleton() {
            spdlog::set_level(spdlog::level::info);
            sink = std::make_shared<MatlabSink_mt>();
            matlab_logger = std::make_shared<spdlog::logger>("matlab", sink);
            spdlog::set_default_logger(matlab_logger);
        }

        static LoggerSingleton& instance()
        {
            static LoggerSingleton just_one;
            return just_one;
        }

        void set_log_level(spdlog::level::level_enum level){
            matlab_logger->set_level(level);
        };

        void set_log_pattern(std::string &pattern){
            matlab_logger->set_pattern(pattern);
        };

        spdlog::level::level_enum get_log_level() {
            return matlab_logger->level();
        };
    };

LoggerSingleton logger = LoggerSingleton();

// logger->set_log_level(spdlog::level::warn);
// auto get_logger(){return logger;};

// void set_log_level(spdlog::level::level_enum level){
//     get_logger().set_log_level(level);
// };

#endif // CLODE_MEX_LOGGING_H