#include <string>
#include <map>
#include <utility>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include "Logger.hpp"
#include "ConsoleAppender.hpp"
#include <iostream>

// --------------------------------------------------------------------------------------------------------------------------------

class LoggerManagerRegistry
{
    //typedef std::pair<std::regex, std::shared_ptr<LoggerManager> > entry_type;

public:
    //void register_logger_manager(char const* logger_name_pattern, std::shared_ptr<LoggerManager> const& manager)
    //{
    //    managers_.push_back(entry_type(std::regex(logger_name_pattern), manager));
    //}

    std::shared_ptr<LoggerManager> get_default_logger_manager() const
    {
        return default_manager_;
    }

    std::shared_ptr<LoggerManager> operator()(char const* logger_name) const
    {
        if (!logger_name) return default_manager_;

        //char const* const logger_name_end(logger_name + std::strlen(logger_name));
        //BOOST_FOREACH(entry_type const& i, managers_)
        //{
        //    if (boost::regex_match(logger_name, logger_name_end, i.first))
        //        return i.second;
        //}

        assert(default_manager_.get());
        return default_manager_;
    }

    LoggerManagerRegistry() : default_manager_(new LoggerManager("default"))
    {
        default_manager_->add_appender(std::make_shared<ConsoleAppender>());
    }

private:
    //    std::vector<entry_type> managers_;
    std::shared_ptr<LoggerManager> default_manager_;
};

// --------------------------------------------------------------------------------------------------------------------------------

static LoggerManagerRegistry registry;

// --------------------------------------------------------------------------------------------------------------------------------

//void LoggerManager::register_logger_manager(char const* logger_name_pattern, std::shared_ptr<LoggerManager> const& manager)
//{
//    registry.register_logger_manager(logger_name_pattern, manager);
//}

std::shared_ptr<LoggerManager> LoggerManager::get_logger_manager(char const* logger_name_pattern)
{
    return registry(logger_name_pattern);
}

std::shared_ptr<LoggerManager> Logger::manager() const
{
    const_cast<Logger*>(this)->ensure_initialized();
    return manager_;
}

// --------------------------------------------------------------------------------------------------------------------------------

Logger& Logger::get_logger(char const* name)
{
    typedef std::map<std::string, std::unique_ptr<Logger>> loggers_map;
    static loggers_map loggers;

    std::string _name(name);
    auto i(loggers.insert(loggers_map::value_type(_name, nullptr)));

    if (i.second) (*i.first).second = std::unique_ptr<Logger>(new Logger(registry, name));

    return *(*i.first).second;
}

char const* Logger::stringize_error_level(loglevel lv)
{
    static char const* names[] = { "OFF", "DEBUG", "INFO", "WARN", "ERROR", "FATAL" };
    return static_cast<std::size_t>(lv) >= sizeof(names) / sizeof(*names) ? "???" : names[static_cast<std::size_t>(lv)];
}

// --------------------------------------------------------------------------------------------------------------------------------

void Logger::level(loglevel level)
{
    ensure_initialized();
    level_ = level;
}

Logger::loglevel Logger::level() const
{
    const_cast<Logger*>(this)->ensure_initialized();
    return level_;
}

void Logger::logv(loglevel lv, char const* format, va_list ap)
{
    ensure_initialized();

    if (lv < level_) return;

    char buf[1024];
    std::vsnprintf(buf, sizeof(buf), format, ap);

    const char* chunks[] = { buf, nullptr };
    char const* const name = name_.c_str();
    std::for_each(appenders_.begin(), appenders_.end(), [lv, name, &chunks](std::shared_ptr<LogAppender> const& app){ (*app)(lv, name, chunks); });
}

void Logger::flush()
{
    ensure_initialized();
    std::for_each(appenders_.begin(), appenders_.end(), [](std::shared_ptr<LogAppender> pla){ pla->flush(); });
}

inline void Logger::ensure_initialized()
{
    if (manager_) return;

    std::shared_ptr<LoggerManager> manager(registry_(name_.c_str()));
    std::vector<std::shared_ptr<LogAppender> > appenders(manager->appenders());
    level_ = manager->level();
    appenders_.swap(appenders);
    manager->manage(this);
    manager_ = manager;

}

// --------------------------------------------------------------------------------------------------------------------------------

void LoggerManager::level(Logger::loglevel level)
{
    /* synchronized { */
    level_ = level;
    std::for_each(managed_loggers_.begin(), managed_loggers_.end(), [level](Logger* const log){ log->level(level); });
    /* } */
}

std::vector<std::shared_ptr<LogAppender>> const& LoggerManager::appenders() const
{
    /* synchronized() { */
    return appenders_;
    /* } */
}

// --------------------------------------------------------------------------------------------------------------------------------

void LoggerManager::add_appender(std::shared_ptr<LogAppender> const& appender)
{
    /* synchronized() { */
    appenders_.push_back(appender);
    /* } */
}

//LoggerManager::LoggerManager(char const* name, Logger::loglevel level) : name_(name), level_(level) {}

void LoggerManager::manage(Logger* logger)
{
    /* synchronized { */
    managed_loggers_.insert(logger);
    /* }} */
}

// --------------------------------------------------------------------------------------------------------------------------------
