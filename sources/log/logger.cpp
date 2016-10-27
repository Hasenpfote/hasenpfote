#include <cassert>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <unordered_map>
#include "appender.h"
#include "logger.h"

namespace hasenpfote{ namespace log{

static std::string GetTimestamp(std::string& format)
{
    const auto now = std::chrono::system_clock::now();
    const auto time = std::chrono::system_clock::to_time_t(now);
    std::tm tm;
#ifdef _WIN32
    localtime_s(&tm, &time);
#else
#error Platform not supported.
#endif
    std::stringstream ss;
    ss << std::put_time(&tm, format.c_str());
    return ss.str();
}

static std::string SeverityToString(Logger::Severity severity)
{
    switch(severity){
    case Logger::Severity::Verbose:
        return "V";
    case Logger::Severity::Debug:
        return "D";
    case Logger::Severity::Info:
        return "I";
    case Logger::Severity::Warning:
        return "W";
    case Logger::Severity::Error:
        return "E";
    case Logger::Severity::Fatal:
        return "F";
    default:
        return "U";
    }
}

class Logger::Impl final
{
public:
    Impl()
    {
        //log = null_log;
        format = "%Y/%m/%d %T";
    }
    ~Impl() = default;

    Impl(const Impl&) = delete;
    Impl& operator = (const Impl&) = delete;
    Impl(Impl&&) = delete;
    Impl& operator = (Impl&&) = delete;

    void AddAppender(const std::type_index& index, const std::shared_ptr<IAppender>& appender);
    void RemoveAppender(const std::type_index& index);
    void SetSeverity(Severity severity);
    void Log(Severity severity, const std::string& filename, int line, const std::string& message);
    void SetTimestampFormat(const std::string& format);

private:
    std::unordered_map<std::type_index, std::shared_ptr<IAppender>> appender;
    Severity severity;
    std::string format;
    std::mutex m;
};

void Logger::Impl::AddAppender(const std::type_index& index, const std::shared_ptr<IAppender>& appender)
{
    std::lock_guard<std::mutex> lg(m);
    this->appender[index] = appender;
}

void Logger::Impl::RemoveAppender(const std::type_index& index)
{
    std::lock_guard<std::mutex> lg(m);
    decltype(appender)::const_iterator it = appender.find(index);
    if(it != appender.cend()){
        appender.erase(it);
    }
}

void Logger::Impl::SetSeverity(Severity severity)
{
    std::lock_guard<std::mutex> lg(m);
    this->severity = severity;
}

void Logger::Impl::Log(Severity severity, const std::string& filename, int line, const std::string& message)
{
    std::lock_guard<std::mutex> lg(m);
    if(appender.empty() || (severity < this->severity))
        return;
    auto s = SeverityToString(severity);
    std::ostringstream oss;
    oss << GetTimestamp(format) << " " << s << "/" << filename << "(" << line << ") - " << message;
    const auto buffer = oss.str();
    for(auto& pair : appender){
        pair.second->Write(buffer);
    }
}

void Logger::Impl::SetTimestampFormat(const std::string& format)
{
    std::lock_guard<std::mutex> lg(m);
    this->format = format;
};

Logger::Logger()
    : pimpl(std::make_unique<Impl>())
{
}

Logger::~Logger() = default;

void Logger::AddAppender(const std::type_index& index, const std::shared_ptr<IAppender>& appender)
{
    assert(appender);
    pimpl->AddAppender(index, appender);
}

void Logger::RemoveAppender(const std::type_index& index)
{
    pimpl->RemoveAppender(index);
}

void Logger::SetSeverity(Severity severity)
{
    pimpl->SetSeverity(severity);
}

void Logger::Log(Severity severity, const std::string& filename, int line, const std::string& message)
{
    pimpl->Log(severity, filename, line, message);
}

void Logger::SetTimestampFormat(const std::string& format)
{
    pimpl->SetTimestampFormat(format);
}

}}