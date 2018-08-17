#include <Logger.h>

ClassImp(MDImputation_ns::Logger)
ClassImp(MDImputation_ns::ConsoleLogger)
ClassImp(MDImputation_ns::LoggerManager)

namespace MDImputation_ns {


int LoggerManager::m_target;
Logger* LoggerManager::m_logger= 0;

}//close namespace
