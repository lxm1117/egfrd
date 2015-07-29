#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <cstdarg>
#include <set>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>

#include "Defs.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class LogAppender;
class LoggerManager;
class LoggerManagerRegistry;

// --------------------------------------------------------------------------------------------------------------------------------

class Logger
{
public:
    enum class loglevel
    {
        L_OFF = 0,
        L_DEBUG = 1,
        L_INFO = 2,
        L_WARNING = 3,
        L_ERROR = 4,
        L_FATAL = 5
    };

    void level(loglevel level) { ensure_initialized(); level_ = level; }

    loglevel level() const { const_cast<Logger*>(this)->ensure_initialized(); return level_; }

    char const* name() const { return name_.c_str(); }

    std::shared_ptr<LoggerManager> manager() const { const_cast<Logger*>(this)->ensure_initialized(); return manager_; };

    void debug(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(loglevel::L_DEBUG, format, ap);
        va_end(ap);
    }

    void info(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(loglevel::L_INFO, format, ap);
        va_end(ap);
    }

    void warn(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(loglevel::L_WARNING, format, ap);
        va_end(ap);
    }

    void error(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(loglevel::L_ERROR, format, ap);
        va_end(ap);
    }

    void fatal(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(loglevel::L_FATAL, format, ap);
        va_end(ap);
    }

    void log(loglevel lv, char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(lv, format, ap);
        va_end(ap);
    }

    GF_CLASS void logv(loglevel lv, char const* format, va_list ap);

    GF_CLASS void flush();

    Logger() = delete;
    Logger(const Logger&) = delete;
    Logger(LoggerManagerRegistry const& registry, char const* name) : registry_(registry), name_(name), manager_() {}

    GF_CLASS static Logger& get_logger(char const* name);

    GF_CLASS static char const* stringize_error_level(loglevel lv);

private:
    GF_CLASS void ensure_initialized();

protected:
    LoggerManagerRegistry const& registry_;
    std::string const name_;
    std::shared_ptr<LoggerManager> manager_;
    loglevel level_;
    std::vector<std::shared_ptr<LogAppender> > appenders_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class LoggerManager
{
    friend class Logger;

public:
    void level(Logger::loglevel level)
    {
        level_ = level;
        std::for_each(managed_loggers_.begin(), managed_loggers_.end(), [level](Logger* const log){ log->level(level); });
    }

    Logger::loglevel level() const { return level_; };

    char const* name() const { return name_.c_str(); }

    std::vector<std::shared_ptr<LogAppender> > const& appenders() const { return appenders_; }

    void add_appender(std::shared_ptr<LogAppender> const& appender) { appenders_.push_back(appender); };

    LoggerManager() = delete;
    LoggerManager(const LoggerManager&) = delete;
    LoggerManager(char const* name, Logger::loglevel level = Logger::loglevel::L_INFO) : name_(name), level_(level) {}

    GF_CLASS static std::shared_ptr<LoggerManager> get_logger_manager(char const* logger_name_patern);

protected:
    void manage(Logger* logger) { managed_loggers_.insert(logger); }

    std::string const name_;
    Logger::loglevel level_;
    std::set<Logger*> managed_loggers_;
    std::vector<std::shared_ptr<LogAppender> > appenders_;
};

// --------------------------------------------------------------------------------------------------------------------------------

class LogAppender
{
public:
    virtual ~LogAppender() {};

    virtual void flush() = 0;

    virtual void operator()(Logger::loglevel lv, char const* name, char const** chunks) = 0;
};

// --------------------------------------------------------------------------------------------------------------------------------

#define LOG_DEBUG(args) if (log_.level() == Logger::loglevel::L_DEBUG) log_.debug args

#define LOG_INFO(args) if (log_.level() <= Logger::loglevel::L_INFO) log_.info args

#define LOG_WARNING(args) if (log_.level() <= Logger::loglevel::L_WARNING) log_.warn args

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* LOGGER_HPP */
