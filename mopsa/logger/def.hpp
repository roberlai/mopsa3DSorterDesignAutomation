#ifndef MOPSA_LOGGER_DEF_H
#define MOPSA_LOGGER_DEF_H

#define ASSERT(EXPR) \
  assert_handler((EXPR), #EXPR, __FILE__, __func__, __LINE__, "");

#define ASSERT_MESS(EXPR, message) \
  assert_handler((EXPR), #EXPR, __FILE__, __func__, __LINE__, message); 

#define LOG(SENSITIVTY) Logger::get_instance().log(SENSITIVTY)

#endif

