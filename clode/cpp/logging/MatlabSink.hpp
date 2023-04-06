//
// Created by Wolf on 27/08/2022.
//

#ifndef CLODE_MATLABSINK_HPP
#define CLODE_MATLABSINK_HPP

#include <iostream>

#include <spdlog/sinks/base_sink.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define mexPrintf printf
#define mexEvalString printf
#else
#define mexPrintf printf
#define mexEvalString printf
#endif


template<typename Mutex>
class MatlabSink : public spdlog::sinks::base_sink <Mutex>
        {
protected:
    void sink_it_(const spdlog::details::log_msg& msg) override
    {
        // log_msg is a struct containing the log entry info like level, timestamp, thread id etc.
        // msg.raw contains pre-formatted log

        // If needed (very likely but not mandatory), the sink formats the message before sending it to its final destination:
        spdlog::memory_buf_t formatted;
        spdlog::sinks::base_sink<Mutex>::formatter_->format(msg, formatted);
        mexPrintf(fmt::to_string(formatted).c_str());
    }

    void flush_() override
    {
        mexEvalString("drawnow;");
    }
};

#include "spdlog/details/null_mutex.h"
#include <mutex>
using MatlabSink_mt = MatlabSink<std::mutex>;
using MatlabSink_st = MatlabSink<spdlog::details::null_mutex>;

#endif //CLODE_MATLABSINK_HPP
