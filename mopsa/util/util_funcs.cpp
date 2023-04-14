#include <mopsa/util/util_funcs.h>
#include <cmath>
#include <string>

namespace mopsa
{

bool 
float_equal(
  double a, 
  double b
)
{
  return std::fabs(a - b) < 1e-7;
}

std::string 
to_string(double val)
{

  std::string res = std::to_string(val);
  while(res.size()) {
    if(res.back() == '0') res.pop_back();
    else break;
  }

  if(res.size() and res.back() == '.') res.pop_back();

  if(res.size() == 0) res = "0";

  return res;
}

}
