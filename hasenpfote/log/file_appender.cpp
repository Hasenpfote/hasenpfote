#include "file_appender.h"

namespace hasenpfote{ namespace log{

FileAppender::FileAppender(const std::tr2::sys::path& filepath)
    : ofs(std::make_unique<std::ofstream>(filepath))
{
}

void FileAppender::Write(const std::string& buffer)
{
    *ofs << buffer.c_str() << std::endl;
}

}}