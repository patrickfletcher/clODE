//
// Created by Wolf on 20/09/2022.
//

#ifndef CLODE_PYTHONSINK_HPP
#define CLODE_PYTHONSINK_HPP

#include <spdlog/sinks/base_sink.h>

#include "spdlog/details/null_mutex.h"
#include <mutex>

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>

//TODO: Windows produces extra carriage returns [workaround: use stdout_color_mt?]
// https://github.com/gabime/spdlog/issues/1675

//PYBIND11_MODULE(spdlog, m) {

    template<typename Mutex>
    class PythonSink : public spdlog::sinks::base_sink<Mutex> {
    protected:
        void sink_it_(const spdlog::details::log_msg &msg) override {
            // log_msg is a struct containing the log entry info like level, timestamp, thread id etc.
            // msg.raw contains pre-formatted log

            // If needed (very likely but not mandatory), the sink formats the message before sending it to its final destination:
            spdlog::memory_buf_t formatted;
            spdlog::sinks::base_sink<Mutex>::formatter_->format(msg, formatted);
            pybind11::print(fmt::to_string(formatted));
        }

        void flush_() override {
            //pybind11::print("", "flush"_a = true);
        }
    };

    using PythonSink_mt = PythonSink<std::mutex>;
    using PythonSink_st = PythonSink<spdlog::details::null_mutex>;
//}

#endif //CLODE_PYTHONSINK_HPP
