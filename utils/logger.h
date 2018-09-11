#ifndef LOGGER_H
#define LOGGER_H
#include "log.h"

#include "common/stl_header.h"
//todo: scoped lock for multi-threading debug
#ifdef USE_LOGGING
#define SET_LOG_FILE(file) \
        std::cout<<"\n************************\n"; \
        std::cout<<"\nLogging enabled\n"; \
        std::cout<<"Logging to "<<file; \
        std::cout<<"\n\n************************\n"; \
        Output2FILE::filename = file
#define SET_LOG_LEVEL(level) \
        FILELog::ReportingLevel() = FILELog::FromString(level)
#define MARK_FUN(msg) \
		FunLogger dontmatchmedontdontdontdontdontdontmatchme(msg)
#define LOG_DEBUG5(msg) \
        FILE_LOG(logDEBUG5) << msg
#define LOG_DEBUG4(msg) \
        FILE_LOG(logDEBUG4) << msg
#define LOG_DEBUG3(msg) \
        FILE_LOG(logDEBUG3) << msg
#define LOG_DEBUG2(msg) \
        FILE_LOG(logDEBUG2) << msg
#define LOG_DEBUG1(msg) \
        FILE_LOG(logDEBUG1) << msg
#define LOG_DEBUG(msg) \
        FILE_LOG(logDEBUG) << msg
#define LOG_INFO(msg) \
        FILE_LOG(logINFO) << msg
#define LOG_WARN(msg) \
        FILE_LOG(logWARNING) << msg
#define LOG_ERROR(msg) \
        FILE_LOG(logERROR) << msg
#define LOG_DONE() \
        FILELog::logging_finished()
#else
#define SET_LOG_FILE(file)
#define SET_LOG_LEVEL(level)
#define MARK_FUN(msg)
#define LOG_DEBUG5(msg)
#define LOG_DEBUG4(msg)
#define LOG_DEBUG3(msg)
#define LOG_DEBUG2(msg)
#define LOG_DEBUG1(msg)
#define LOG_DEBUG(msg)
#define LOG_INFO(msg)
#define LOG_WARN(msg)
#define LOG_ERROR(msg)
#define LOG_DONE()
#endif

class FunLogger{
public:
    FunLogger();
    FunLogger(const char* msg)
    : m(msg) {
        LOG_DEBUG1("Entering "<<m);
    }
    FunLogger(std::string msg)
    : m(msg) {
        LOG_DEBUG1("Entering "<<m);
    }
    void setMsg(const char* msg) {
        m = std::string(msg);
    }
    ~FunLogger() {
        LOG_DEBUG1("Leaving "<<m);
    }
private:
    std::string m;
};

#endif /*LOGGER_H*/
