#include <mopsa/logger/debug.hpp>

#include <iostream>
#include <cstdlib>

namespace mopsa 
{

void 
mopsa_exit(int code) 
{
  std::exit(code);
}

void 
assert_handler(
  int expr, 
  const char *assertion, 
  const char *fname, 
  const char *func, 
  int line,
  const std::string &message = ""
) 
{
  if (!(expr)) {
    std::cerr << "\n\nAssertion \"" << assertion << "\" fail at function \"" 
      << func << "\"\n"
      << "- At file \"" 
      << fname 
      << "\", line number "
      << line << "\n\n";

    if(message.size()) std::cerr << message << '\n';
    mopsa_exit(-1);
  }
}

}
