#include <iostream>
#include <limits>

using std::cerr;
using std::endl;

#include "iw_tabular_data.h"

template <typename T>
IW_Tabular_Data<T>::IW_Tabular_Data()
{
  _nrows = 0;
  _ncols = 0;
  _zdata = nullptr;

  _has_header = 0;
  _header = nullptr;

  _missing_value = "N/A";

  _first_column_is_identifier = 0;
  _id = nullptr;

  return;
}

template<typename T>
IW_Tabular_Data<T>::~IW_Tabular_Data()
{
  if (nullptr != _zdata)
    delete [] _zdata;

  if (nullptr != _header)
    delete [] _header;

  if (nullptr != _id)
    delete [] _id;

  _nrows = -1;
  _ncols = -1;

  return;
}

template<typename T>
int
IW_Tabular_Data<T>::resize (const int nr, const int nc)
{
  if (nr <= 0 || nc <= 0 || nr > (std::numeric_limits<int>::max() - 31) / nc) {
    return 0;
  }

  if (nr == _nrows && nc == _ncols)
    return 1;

  if (nullptr != _zdata)
    delete [] _zdata;

  delete [] _header;
  _header = nullptr;
  delete [] _id;
  _id = nullptr;
  _zdata = new T[nr * nc]{};

  _nrows = nr;
  _ncols = nc;

  int nb = _nrows * _ncols;
  if (0 != nb % 32)
    nb = (nb /32 + 1) * 32;

  _missing.allocate_space_for_bits(nb);
  _missing.clear();

  return 1;
}

template <typename T>
int
IW_Tabular_Data<T>::build(const char * fname,
                          const char sep)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "IW_Tabular_Data::build:cannot open '" << fname << "'\n";
    return 0;
  }

  return build(input, sep);
}

template <typename T>
int
IW_Tabular_Data<T>::_determine_number_columns(const const_IWSubstring & buffer,
                                              const char sep)
{
  if (' ' == sep)
    _ncols =  buffer.nwords();
  else
    _ncols = buffer.nwords_single_delimiter(sep);

  if (_first_column_is_identifier)
    _ncols--;

  return _ncols > 0;
}

static int
get_next_token(const const_IWSubstring & buffer,
              const_IWSubstring & token,
              int & i,
              const char sep)
{
  if (' ' == sep)
    return buffer.nextword(token, i);

  return buffer.nextword_single_delimiter(token, i, sep);
}

template <typename T>
int
IW_Tabular_Data<T>::_parse_row(const const_IWSubstring & buffer,
                               const char sep,
                               const int row)
{
  T * d = _zdata + row * _ncols;

  int i = 0;
  const_IWSubstring token;

  if (_first_column_is_identifier)
  {
    if (! get_next_token(buffer, token, i, sep))
      return 0;
    _id[row] = token;
  }

  int c = 0;

  for ( ;c < _ncols && get_next_token(buffer, token, i, sep); ++c)
  {
    _missing.set(row * _ncols + c, 0);
    if (token.numeric_value(d[c]))
      continue;

    if (_missing_value == token)
    {
      d[c] = T{};
      _missing.set(row * _ncols + c, 1);
      continue;
    }

    cerr << "IW_Tabular_Data::_parse_row:invalid numeric '" << token << "'\n";
    return 0;
  }

  if (c != _ncols || get_next_token(buffer, token, i, sep)) {
    cerr << "IW_Tabular_Data::_parse_row:incorrect number of columns\n";
    return 0;
  }

  return 1;
}

template <typename T>
int
IW_Tabular_Data<T>::_parse_header(const const_IWSubstring & buffer, const char sep)
{
  int i = 0;
  const_IWSubstring token;

  if (nullptr != _header)
    delete [] _header;

  _header = new IWString[_ncols];

  if (_first_column_is_identifier && !get_next_token(buffer, token, i, sep)) {
    return 0;
  }
  if (_first_column_is_identifier) {
    _identifier_header = token;
  }
  for (int c = 0; c < _ncols; ++c) {
    if (!get_next_token(buffer, token, i, sep)) {
      return 0;
    }
    _header[c] = token;
  }
  if (get_next_token(buffer, token, i, sep)) {
    return 0;
  }

  return 1;
}

template <typename T>
int
IW_Tabular_Data<T>::_determine_rows_and_columns(iwstring_data_source & input,
                                              const char sep,
                                              const_IWSubstring & buffer)
{
  if (! input.next_record(buffer))
  {
    cerr << "IW_Tabular_Data::_determine_number_columns:empty file\n";
    return 0;
  }

  _nrows = input.records_remaining() + 1;

  if (! _determine_number_columns(buffer, sep))
  {
    cerr << "IW_Tabular_Data::_determine_number_columns:cannot determine column count\n";
    cerr << buffer << endl;
    return 0;
  }

  if (_has_header)
  {
    if (1 == _nrows)
    {
      cerr << "IW_Tabular_Data::_determine_number_columns:only header record\n";
      return 0;
    }

    _nrows--;

    if (! _parse_header(buffer, sep))
    {
      cerr << "IW_Tabular_Data::_determine_number_columns:cannot parse header\n";
      cerr << buffer << endl;
      return 0;
    }

    if (! input.next_record(buffer))
      return 0;
  }

#ifdef DEBUG_DETERMINE_NUMBER_COLUMNS
  cerr << "Based on " << buffer << " data has " << _nrows << " rows and " << _ncols << " cols\n";
  cerr << _has_header << endl;
#endif

  if (_nrows <= 0 || _ncols <= 0 ||
      _nrows > (std::numeric_limits<int>::max() - 31) / _ncols) {
    return 0;
  }
  _zdata = new T[_nrows * _ncols]{};
  _missing.allocate_space_for_bits(((_nrows * _ncols + 31) / 32) * 32);
  _missing.clear();

  if (_first_column_is_identifier)
    _id = new IWString[_nrows];

  return 1;
}

template <typename T>
int
IW_Tabular_Data<T>::build(iwstring_data_source & input,
                          const char sep)
{
  // A new build replaces the previous table, including its dimensions and metadata.
  delete [] _zdata;
  _zdata = nullptr;
  delete [] _header;
  _header = nullptr;
  delete [] _id;
  _id = nullptr;
  _nrows = _ncols = 0;

  const_IWSubstring buffer;
  if (!_determine_rows_and_columns(input, sep, buffer)) {
    return 0;
  }
  for (int r = 0; r < _nrows; ++r) {
    if (!_parse_row(buffer, sep, r)) {
      return 0;
    }
    if (r + 1 < _nrows && !input.next_record(buffer)) {
      return 0;
    }
  }
  if (input.next_record(buffer)) {
    cerr << "IW_Tabular_Data::build:too many rows\n";
    return 0;
  }

  return 1;
}

template <typename T> template <typename O>
int
IW_Tabular_Data<T>::do_write(O & output, const char sep) const
{
  if (nullptr != _header)
  {
    if (_first_column_is_identifier) {
      output << _identifier_header << sep;
    }
    output << _header[0];
    for (int i = 1; i < _ncols; ++i)
    {
      output << sep << _header[i];
    }

    output << '\n';
  }

  for (int r = 0; r < _nrows; ++r)
  {
    const T * row = _zdata + r * _ncols;

    if (_first_column_is_identifier && _id != nullptr) {
      output << _id[r] << sep;
    }
    for (int c = 0; c < _ncols; ++c) {
      if (c > 0) {
        output << sep;
      }
      if (_missing.is_set(r * _ncols + c)) {
        output << _missing_value;
      } else {
        output << row[c];
      }
    }
    output << '\n';

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

template <typename T>
int
IW_Tabular_Data<T>::remove_column(const int rmcol)
{
  if (rmcol < 0 || rmcol >= _ncols || _ncols == 1) {
    return 0;
  }
  const int old_ncols = _ncols;
  int dest = 0;
  for (int r = 0; r < _nrows; ++r) {
    for (int c = 0; c < old_ncols; ++c) {
      if (c == rmcol) {
        continue;
      }
      const int source = r * old_ncols + c;
      _zdata[dest] = _zdata[source];
      _missing.set(dest, _missing.is_set(source));
      ++dest;
    }
  }
  if (_header != nullptr) {
    for (int c = rmcol + 1; c < old_ncols; ++c) {
      _header[c - 1] = _header[c];
    }
  }
  _ncols--;

  return 1;
}

template class IW_Tabular_Data<double>;
template class IW_Tabular_Data<int>;
template int IW_Tabular_Data<double>::do_write(IWString_and_File_Descriptor &, const char) const;
