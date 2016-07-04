#include <iostream>
#include "console_appender.h"

namespace hasenpfote{ namespace log{

void ConsoleAppender::Write(const std::string& buffer)
{
    std::cout << buffer.c_str() << std::endl;
}

}}