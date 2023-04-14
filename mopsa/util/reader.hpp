#ifndef MOPSA_UTIL_READER_H
#define MOPSA_UTIL_READER_H

#include <filesystem>
#include <unordered_set>
#include <vector>

namespace mopsa 
{

class Reader {

public:
  Reader(const std::string& specialchar, const std::string& sepchar);
  ~Reader();

  bool open_file(const std::filesystem::path& fpath);

  virtual std::string read_token(bool ignore_comment=true);
  virtual std::string read_comment();
  virtual bool read_string(std::string&); // "string"
  virtual bool read_int(int &val);
  virtual bool read_double(double &val);
  virtual bool read_boolean(bool &val);

protected:

  bool _is_end();
  char _get_char();
  char _cur_char();

  bool _is_special_char();
  bool _is_sep_char();
  virtual bool _is_comment_prefix();

  virtual bool _expect_token(const std::string& exp, const std::string& expto);
  virtual void _fail_to_parse(const std::string& message);

  char* _buffer{nullptr};
  size_t _length{0};
  size_t _read_cur{0};
  size_t _line_no{1};
  size_t _column{1};
  std::filesystem::path _filename;

  std::unordered_set<char> _special_set, _sep_set;
};

}

#endif
