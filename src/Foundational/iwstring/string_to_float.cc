// Those IWString member functions associated with string to float conversions

#ifdef USE_DRAGONBOX
#include "dragonbox/dragonbox_to_chars.h"
#endif

#include "iwstring.h"

static int float_precision = 7;

void
set_default_iwstring_float_concatenation_precision (int s)
{
  float_precision = s;
}

static int double_precision = 10;

void
set_default_iwstring_double_concatenation_precision (int s)
{
  double_precision = s;
}

void
IWString::append_number(float f, int fprecision)
{
  char buffer[32];
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
  gcvt(static_cast<double>(f), fprecision, buffer);
#pragma GCC diagnostic pop
  resizable_array<char>::add(buffer, static_cast<int>(::strlen(buffer)));
  return;
}

void
IWString::append_number(float f) {
  append_number(f, float_precision);
}

void
IWString::append_number(double d)
{
  append_number(d, double_precision);
}

void
IWString::append_number(double d, int dprecision)
{
  char buffer[32];
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
  gcvt(d, dprecision, buffer);
#pragma GCC diagnostic pop
  resizable_array<char>::add(buffer, static_cast<int>(::strlen(buffer)));
  return;
}

void
IWString::append_number(float f, const char * fformat)
{
  char buffer[32];

  int nchars = IW_SPRINTF(buffer, fformat, f);

  resizable_array<char>::add(buffer, nchars);

  return;
}

#ifdef USE_DRAGONBOX
void
IWString::append_number_dragonbox(float f) {
  if (f == 0.0f) {
    resizable_array<char>::add('0');
    return;
  }

  constexpr int kBufferLength = 1 + // for '\0'
    jkj::dragonbox::max_output_string_length<jkj::dragonbox::ieee754_binary64>;

  char buffer[kBufferLength];

  char* end_ptr = jkj::dragonbox::to_chars_n(x, buffer);

}
#endif
