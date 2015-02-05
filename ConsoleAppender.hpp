#ifndef CONSOLEAPPENDER_HPP
#define CONSOLEAPPENDER_HPP

#include <string>
#include "Logger.hpp"

class ConsoleAppender : public LogAppender
{
public:
    virtual ~ConsoleAppender() {};

    virtual void flush();
    virtual void operator()(enum Logger::level lv, char const* name, char const** chunks);
};

#endif /* CONSOLEAPPENDER_HPP */
