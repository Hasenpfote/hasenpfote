#include "file_appender.h"

namespace hasenpfote{ namespace log{

FileAppender::FileAppender(const std::string& filepath)
    : ofs(std::make_unique<std::ofstream>(filepath))
{
}

#if (__cplusplus > 201402L) || (defined(_MSC_VER) && (_MSVC_LANG > 201402L))
FileAppender::FileAppender(const std::filesystem::path& filepath)
    : ofs(std::make_unique<std::ofstream>(filepath))
{
}
#endif

void FileAppender::Write(const std::string& buffer)
{
    *ofs << buffer.c_str() << std::endl;
}

}}