#include <sstream>
#include <cassert>
#include "rolling_file_appender.h"

namespace hasenpfote{ namespace log{

static std::tr2::sys::path BuildFilePath(const std::tr2::sys::path& filepath, int number)
{
    auto new_filepath(filepath);
    new_filepath.replace_extension("");
    std::ostringstream oss;
    oss << "." << number << ".";
    new_filepath += oss.str();
    new_filepath.replace_extension(filepath.extension());
    return new_filepath;
}

static std::streampos GetNextPosition(std::ofstream& ofs, const std::string& buffer)
{
    auto pos = ofs.tellp();
#ifdef _WIN32
    pos += buffer.length() + 2; // CRLF
#else
    pos += buffer.length() + 1;
#endif
    return pos;
}

RollingFileAppender::RollingFileAppender(const std::tr2::sys::path& filepath, int max_files, std::size_t max_file_size)
    : ofs(std::make_unique<std::ofstream>()), filepath(filepath), max_files(max_files), max_file_size(max_file_size)
{
    assert(max_files > 0);
}

void RollingFileAppender::Write(const std::string& buffer)
{
    if(!ofs->is_open()){
        ofs->open(filepath);
    }
    else
    if(max_files > 1 && (GetNextPosition(*ofs, buffer) >= max_file_size)){
        ofs->close();
        // rolling files.
        std::error_code ret;
        for(auto i = max_files - 2; i >= 0; i--){
            auto current_filepath = (i > 0)? BuildFilePath(filepath, i) : filepath;
            if(std::tr2::sys::exists(current_filepath)){
                auto next_filepath = BuildFilePath(filepath, i+1);
                std::tr2::sys::rename(current_filepath, next_filepath, ret);
            }
        }
        ofs->open(filepath);
    }
    *ofs << buffer.c_str() << std::endl;
}

}}