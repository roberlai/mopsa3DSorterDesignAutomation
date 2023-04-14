#ifndef MOPSA_LOGGER_DEBUG_H
#define MOPSA_LOGGER_DEBUG_H

#include <iostream>
#include <cstdlib>

#include <mopsa/logger/def.hpp>

namespace mopsa 
{

extern 
void mopsa_exit(int code);

extern
void assert_handler(
  int expr, 
  const char *assertion, 
  const char *fname, 
  const char *func, 
  int line,
  const std::string &message
);

}

#endif
