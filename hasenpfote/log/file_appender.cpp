#include "file_appender.h"

namespace hasenpfote{ namespace log{

void FileAppender::Write(const std::string& buffer)
{
    *ofs << buffer.c_str() << std::endl;
}

}}