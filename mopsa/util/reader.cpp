#include <mopsa/util/reader.hpp>
#include <mopsa/logger/logger.hpp>

#include <cctype>
#include <fstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

namespace mopsa
{

Reader::Reader(const std::string& specialchar, const std::string& sepchar)
  :_read_cur(0) 
{
  for(const auto& c : specialchar) _special_set.insert(c);
  for(const auto& c : sepchar) _sep_set.insert(c);

  _sep_set.insert('\r');
  _sep_set.insert('\t');
  _sep_set.insert(' ');
}

Reader::~Reader() 
{
  delete _buffer;
}

bool Reader::open_file(const std::filesystem::path& fpath)
{
  std::ifstream in(fpath, std::ios::binary);
  _filename = fpath;

  if(!in.good()) return false;

  in.seekg(0, in.end);
  _length = in.tellg();
  in.seekg(0, in.beg);

  if(_buffer) delete _buffer;

  _buffer = new char[_length+10];
  in.read(_buffer, _length);
  in.close();

  _line_no = 1;
  _column = 1;
  _read_cur = 0;
  //std::cout << "Open " << fpath << ", size = " << _length << '\n';
  return true;
}

bool Reader::_is_end() 
{
  return _read_cur >= _length;
}

bool Reader::_is_comment_prefix() 
{
  if(_read_cur + 1 < _length) {
    return ( _buffer[_read_cur] == '/' and _buffer[_read_cur+1] == '/')
      || ( _buffer[_read_cur] == '/' and _buffer[_read_cur+1] == '*');
  }

  return false;
}

bool Reader::_is_special_char() 
{
  return _special_set.find(_buffer[_read_cur]) != _special_set.end();
}

bool Reader::_is_sep_char() 
{
  return _sep_set.find(_buffer[_read_cur]) != _sep_set.end();
}

std::string Reader::read_token(bool ignore_comment /*= true*/) 
{
  std::string token("");

  while(!_is_end()) {
    while(!_is_end() and _is_sep_char()) _get_char();
    if(_is_end()) break;

    if( _is_comment_prefix() ) {
      if( !ignore_comment ) return read_comment(); 
      else read_comment();
    }
    else if(_is_special_char()) {
      return std::string(1, _get_char());
    }
    //else if(is_number())  return _read_number();
    else {
      while(!_is_end()) {
        if(_is_special_char()) return token;
        if(_is_sep_char()) return token;
        token.push_back(_get_char());
      }
    }
  }

  return token;
}

std::string Reader::read_comment() {

  char c;
  std::string token;
  if(_read_cur + 1 < _length) {

    // //
    if(_buffer[_read_cur] == '/' and _buffer[_read_cur+1] == '/') {
      _get_char(); _get_char(); 
      while(!_is_end()) {
        c = _get_char();
        if(c=='\n') break;
        token.push_back(c);
      }
    } 

    // /*
    else if(_buffer[_read_cur] == '/' and _buffer[_read_cur+1] == '*') {
      _get_char(); _get_char(); 
      while(!_is_end()) {
        c = _get_char();
        if(c == '/' and !token.empty() and token.back() == '*') {
          token.pop_back();
          break;
        }
        token.push_back(c);
      }
    }

    else _fail_to_parse("Fail to parse comment");
  }
  else _fail_to_parse("Fail to parse comment");

  return token;
}

bool Reader::read_string(std::string& val) {

  while(!_is_end() and _is_sep_char()) _get_char();
  val = "";
  char c = _cur_char();
  if(c != '\"') {
    LOG(WARNING) << "Cannot read string at " <<
      _filename << ":" << _line_no << '\n';
    return false;
  }

  bool ok = false;
  _get_char(); // '"'
  while(!_is_end()) {
    c = _get_char();
    if(c == '\"') {
      ok = true;
      break;
    }
    else val.push_back(c);
  }

  if(!ok) {
    LOG(WARNING) << "Fail to read string, expect '\"' at " <<
      _filename << ":" << _line_no << '\n';
  }

  return ok;
}

bool Reader::read_int(int &val) {

  std::string token;
  while(!_is_end() and _is_sep_char()) _get_char();

  char c = _cur_char();
  if(c == '+' or c == '-' or std::isdigit(c)) {
    token = _get_char();
    while(!_is_end()) {
      if(!std::isdigit(_cur_char())) break;
      token += _get_char();
    }
  } else token = _get_char();

  //std::cout << "Token = " << token << std::endl;
  
  try {
    val = std::stoi(token);
  }
  catch(...) {
    LOG(WARNING) << "Cannot convert \'" << token << "\' to int. "
      "Use 0 instead. In " << _filename << ":" << _line_no << '\n';
    val = 0;
    return false;
  }
  return true;
}

bool Reader::read_boolean(bool &val) {

  const auto &token = read_token();

  if(token == "true" or token == "1") {
    val = 1;
  } else if (token == "false" or token == "0") {
    val = 0;
  }
  else {
    LOG(WARNING) << "Cannot convert \'" << token << "\' to boolean. " 
      "Use 0 instread. In " << _filename << ":" << _line_no << '\n';
    val = 0;
    return false;
  }

  return true;
}

bool Reader::read_double(double &val) {

  while(!_is_end() and _is_sep_char()) _get_char();

  std::string token;
  int num_dot = 0;
  char c = _cur_char();
  if(c == '+' or c == '-' or std::isdigit(c) or c == '.') {
    num_dot += (c=='.');
    token = _get_char();
    while(!_is_end()) {
      c = _cur_char();
      if(!std::isdigit(c) and c != '.') break;
      num_dot += (c=='.');
      token += _get_char();
    }
  } else token = _get_char();

  //std::cout << "Token = " << token << std::endl;
  
  try {
    if(num_dot > 1) throw std::invalid_argument("Not Float") ;
    val = std::stod(token);
    //std::cout << "to " << val << std::endl;
  }
  catch(...) {
    LOG(WARNING) << "Cannot convert \'" << token << "\' to double. "
      "Use 0.0 instead. In " << _filename << ":" << _line_no << '\n';
    val = 0.0;
    return false;
  }

  return true;
}

char Reader::_get_char() {

  ASSERT(!_is_end());
  _column += 1;
  char c = _buffer[_read_cur++];
  if(c == '\n') _line_no += 1, _column = 0;
  return c;
}

char Reader::_cur_char() {
  if(_is_end()) return ' ';
  else return _buffer[_read_cur];

}

bool Reader::_expect_token(const std::string& exp, const std::string& expto) {
  if(exp != expto) {
    _fail_to_parse("Expect \'" + exp + "\' to be \'" + expto + "\'");
    return false;
  }
  return true;
}

void Reader::_fail_to_parse(const std::string& message) {
  LOG(ERROR) << "" << message << " at " << _filename
    << ":" << _line_no << ":" << _column << '\n';
  mopsa_exit(-1);
}

}
