#ifndef CONSOLEAPPENDER_HPP
#define CONSOLEAPPENDER_HPP

#include "Logger.hpp"

class ConsoleAppender : public LogAppender
{
public:
    typedef LogAppender base_type;              // GNU need this, but it's not used anywhere, weird? (VS doesn't)

public:
    virtual ~ConsoleAppender() {};

    virtual void flush();
    virtual void operator()(Logger::loglevel lv, char const* name, char const** chunks);
};

#endif /* CONSOLEAPPENDER_HPP */
