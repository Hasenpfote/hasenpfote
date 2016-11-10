#include "file_appender.h"

namespace hasenpfote{ namespace log{

FileAppender::FileAppender(const std::string& filepath)
    : ofs(std::make_unique<std::ofstream>(filepath))
{
}

#if (__cplusplus > 201402L) || (defined(_MSC_VER) && (_MSC_VER > 1900))
#error Function not implemented.
#elif defined(_MSC_VER) && (_MSC_VER == 1900)
FileAppender::FileAppender(const std::tr2::sys::path& filepath)
    : ofs(std::make_unique<std::ofstream>(filepath))
{
}
#endif

void FileAppender::Write(const std::string& buffer)
{
    *ofs << buffer.c_str() << std::endl;
}

}}