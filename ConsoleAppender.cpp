#include <cstdio>
#include <cstring>
#include "ConsoleAppender.hpp"

void ConsoleAppender::operator()(enum class Logger::loglevel lv, char const* name, char const** chunks)
{
    std::fprintf(stderr, "%s: %-8s ", name, Logger::stringize_error_level(lv));
    for (char const** p = chunks; *p; ++p)
        std::fwrite(*p, sizeof(char), std::strlen(*p), stderr);
    std::fputc('\n', stderr);
}

void ConsoleAppender::flush()
{
    std::fflush(stderr);
}
