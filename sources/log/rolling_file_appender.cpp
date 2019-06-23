#include <sstream>
#include <cassert>
#include "rolling_file_appender.h"

namespace hasenpfote{ namespace log{

#if (__cplusplus > 201402L) || (defined(_MSC_VER) && (_MSVC_LANG > 201402L))
static std::filesystem::path BuildFilePath(const std::filesystem::path& filepath, int number)
{
    auto new_filepath(filepath);
    new_filepath.replace_extension("");
    std::ostringstream oss;
    oss << "." << number << ".";
    new_filepath += oss.str();
    new_filepath.replace_extension(filepath.extension());
    return new_filepath;
}

static bool ExistsFile(const std::filesystem::path& filepath)
{
    return std::filesystem::exists(filepath);
}

static bool RenameFile(const std::filesystem::path& old_filepath, const std::filesystem::path& new_filepath)
{
    std::error_code ret;
    std::filesystem::rename(old_filepath, new_filepath, ret);
    return (ret)? false : true;
}
#else
static std::string BuildFilePath(const std::string& filepath, int number)
{
    auto extension = "." + std::to_string(number);
    const auto n = filepath.rfind(".");
    if(n != std::string::npos){
        extension.append(filepath.substr(n));
    }
    return filepath.substr(0, n) + extension;
}

static bool ExistsFile(const std::string& filepath)
{
    std::ifstream infile(filepath);
    return infile.good();
}

static bool RenameFile(const std::string& old_filepath, const std::string& new_filepath)
{
    auto ret = std::rename(old_filepath.c_str(), new_filepath.c_str());
    return (ret)? false : true;
}
#endif

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

RollingFileAppender::RollingFileAppender(const std::string& filepath, int max_files, std::size_t max_file_size)
    : ofs(std::make_unique<std::ofstream>()), filepath(filepath), max_files(max_files), max_file_size(max_file_size)
{
    assert(max_files > 0);
}

#if (__cplusplus > 201402L) || (defined(_MSC_VER) && (_MSVC_LANG > 201402L))
RollingFileAppender::RollingFileAppender(const std::filesystem::path& filepath, int max_files, std::size_t max_file_size)
    : ofs(std::make_unique<std::ofstream>()), filepath(filepath), max_files(max_files), max_file_size(max_file_size)
{
    assert(max_files > 0);
}
#endif

void RollingFileAppender::Write(const std::string& buffer)
{
    if(!ofs->is_open()){
        ofs->open(filepath);
    }
    else
    if(max_files > 1 && (GetNextPosition(*ofs, buffer) >= max_file_size)){
        ofs->close();
        // rolling files.
        for(auto i = max_files - 2; i >= 0; i--){
            auto current_filepath = (i > 0)? BuildFilePath(filepath, i) : filepath;
            if(ExistsFile(current_filepath)){
                auto next_filepath = BuildFilePath(filepath, i + 1);
                auto ret = RenameFile(current_filepath, next_filepath);
                (void)ret;
            }
        }
        ofs->open(filepath);
    }
    *ofs << buffer.c_str() << std::endl;
}

}}