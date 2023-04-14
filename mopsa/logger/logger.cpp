#include <mopsa/logger/logger.hpp>
#include <chrono>
#include <ctime>
#include <iomanip>

namespace mopsa 
{

Logger::Logger() 
{
  _file_log_stream.open("mopsa.log");

  std::time_t now = 
    std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  (*this) << "======== MOPSA Start at " 
    << std::put_time(std::localtime(&now), "%Y-%m-%e %H:%M:%S") 
    << " ==========\n\n";
}

Logger::~Logger() 
{
  std::time_t now = 
    std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  (*this) << "======== MOPSA Finish at " 
    << std::put_time(std::localtime(&now), "%Y-%m-%e %H:%M:%S") 
    << " ==========\n";

  (*this) << "\n****** Record ******\n"; 
  for(const auto& [group, submap] : _records) {
    (*this) << group << ":\n";
    for(const auto& [name, val] : submap) {
      (*this) << "  " << name << ": " << val << '\n';
    }
  }
  (*this) << "****** Record ******\n"; 

  _file_log_stream.close();
}

Logger& 
Logger::get_instance() 
{
  static Logger logger;

  return logger;
}

Logger& 
Logger::log(LoggerSensivity sensivity) 
{
  Logger& logger = Logger::get_instance();

  std::time_t now = 
    std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  logger << "[" << _sensivity_names[sensivity] << "] " 
    << std::put_time(std::localtime(&now), "%m:%d:%H:%M:%S") 
    << " ";

  return logger;
}

void 
Logger::add_record(
  const std::string& group, 
  const std::string& name, 
  int val) 
{
  Logger& logger = Logger::get_instance();

  logger._records[group][name] += val;
}

void 
Logger::set_record(
  const std::string& group, 
  const std::string& name, 
  int val) 
{
  Logger& logger = Logger::get_instance();

  logger._records[group][name] = val;
}

void
Logger::set_stdout(const std::filesystem::path &path)
{
  Logger& logger = Logger::get_instance();
  logger._stdout.open(path);
}

void
Logger::reset_stdout()
{
  Logger& logger = Logger::get_instance();
  logger._stdout.close();
}

}
