#include <sstream>
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