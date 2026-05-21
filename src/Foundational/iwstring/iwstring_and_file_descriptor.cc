#include <stdlib.h>

#include "Foundational/iwmisc/iwconfig.h"
#ifdef _WIN32
#else
#include <unistd.h>
#endif
#include <fcntl.h>

#include <optional>

#include "iwstring.h"
#include "zlib.h"
#include "zstd.h"

#include <sys/stat.h>
#include <sys/types.h>

#include "iwstring_and_file_descriptor.h"

namespace iwstring {

using std::cerr;

IWString_and_File_Descriptor::IWString_and_File_Descriptor() {
  _fd = -9;

  _gzfile = nullptr;

  return;
}

IWString_and_File_Descriptor::IWString_and_File_Descriptor(int f) : _fd(f) {
  _gzfile = nullptr;

  return;
}

void
IWString_and_File_Descriptor::FreeZstdRelated() {
  ZSTD_freeCCtx(_zstd_cctx);
  _zstd_cctx = nullptr;

  delete [] _zstd_buffer;
  _zstd_buffer = nullptr;
  _zstd_buffer_size = 0;
}

IWString_and_File_Descriptor::~IWString_and_File_Descriptor() {
  // zstd compression uses _fd, so it can be closed below.
  cerr << "Destructor buffer contains " << _number_elements << " characters\n";
  if (_zstd_cctx) {
    _zstd_compress_and_write(ZSTD_e_end);
    FreeZstdRelated();
  }

  cerr << "ADFter _zstd_compress_and_write " << _number_elements << " characters\n";

  if (_gzfile) {
    _compress_and_write();
    gzclose(_gzfile);
    _gzfile = nullptr;
    return;
  }

  if (_fd < 0) {
    return;
  }

  // Special handling for stdout, do not close.
  if (_fd == 1) {
    IWString::write(_fd);
    return;
  } 
  
  if (IWString::length() > 0) {
    IWString::write(_fd);
  }

  IW_FD_CLOSE(_fd);
  _fd = -9;

  return;
}

int
IWString_and_File_Descriptor::open(const char* fname) {
  if (_fd > 0) {
    cerr << "IWString_and_File_Descriptor::open:already open\n";
    return 0;
  }

  // I haven't allowed for '>>foo.gz'. Tough

  int strlen_fname = strlen(fname);

  if (1 == strlen_fname && '-' == fname[0]) {
    _fd = 1;
    return 1;
  }

  if (strlen_fname <= 3)  // cannot be of the form 'xxx.gz'
    ;
  else if (0 == ::strcmp(fname + strlen_fname - 3, ".gz")) {
    return _open_gzipd_stream(fname);
  } else if (strlen_fname > 4 && 0 == ::strcmp(fname + strlen_fname - 4, ".zst")) {
      return _open_zstd_stream(fname);
  } else if (strlen_fname > 5 && 0 == ::strcmp(fname + strlen_fname - 5, ".zstd")) {
      return _open_zstd_stream(fname);
  }

  int mode;

  if (strlen_fname > 2 && 0 == ::strncmp(fname, ">>", 2)) {
    fname += 2;
    mode = O_WRONLY | O_APPEND | O_CREAT;
  } else {
    mode = O_WRONLY | O_TRUNC | O_CREAT;
  }

#ifdef _WIN32
  int flags = 0;
#else
  int flags = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
#endif

  // See if we need to expand any shell variables.
  IWString myfname(fname);
  std::optional<IWString> maybe_expanded = myfname.ExpandEnvironmentVariables();
  if (maybe_expanded) {
    _fd = IW_FD_OPEN(maybe_expanded->null_terminated_chars(), mode, flags);
  } else {
    _fd = IW_FD_OPEN(fname, mode, flags);
  }

  if (_fd < 0) {
    cerr << "IWString_and_File_Descriptor::open:cannot open '" << fname << "'\n";
    return 0;
  }

  return _fd;
}

int
IWString_and_File_Descriptor::_open_gzipd_stream(const char* fname) {
#if (__GNUC__ >= 6)
  _gzfile = gzopen64(fname, "wb");
#else
  _gzfile = gzopen(fname, "wb");
#endif

  if (nullptr != _gzfile) {
    return 1;
  }

  cerr << "IWString_and_File_Descriptor::_open_gzipd_stream:cannot open '" << fname
       << "'\n";
  return 0;
}
int
IWString_and_File_Descriptor::_open_zstd_stream(const char* fname) {
  _fd = IW_FD_OPEN(fname, O_WRONLY | O_TRUNC | O_CREAT,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (_fd < 0) {
    cerr << "IWString_and_File_Descriptor::_open_zstd_stream:cannot open '" << fname << "'\n";
    return 0;
  }

  _zstd_cctx = ZSTD_createCCtx();
  if (_zstd_cctx == nullptr) {
    IW_FD_CLOSE(_fd);
    _fd = -9;
    return 0;
  }

  if (size_t rc = ZSTD_CCtx_setParameter(_zstd_cctx, ZSTD_c_compressionLevel,
                                         _zstd_compression_level); ZSTD_isError(rc)) {
    cerr << "ZSTD_CCtx_setParameter failed " << ZSTD_getErrorName(rc) << '\n';
    return 0;
  }

  _zstd_buffer_size = ZSTD_CStreamOutSize();
  if (_zstd_buffer_size == 0) {
    _zstd_buffer_size = 4096;
  }
  _zstd_buffer = new char[_zstd_buffer_size];

  return 1;
}

int
IWString_and_File_Descriptor::_do_resize(int keep_storage) {
  if (IWSFD_KEEP_STORAGE_ON_FLUSH == keep_storage) {
    IWString::resize_keep_storage(0);
  } else if (0 == keep_storage) {
    IWString::resize(0);
  } else {
    cerr << "IWString_and_File_Descriptor::_do_resize:unrecognised option "
         << keep_storage << '\n';
    return 0;
  }

  return 1;
}

int
IWString_and_File_Descriptor::flush(int keep_storage) {
  if (0 == IWString::length()) {
    return 1;
  }

  if (_gzfile) {
    _compress_and_write();
  } else if (_zstd_cctx) {
    _zstd_compress_and_write(ZSTD_e_continue);
  } else if (_fd < 0) {
    cerr << "IWString_and_File_Descriptor::flush:no file descriptor\n";
    return 0;
  } else {
    IWString::write(_fd);
  }

  return _do_resize(IWSFD_KEEP_STORAGE_ON_FLUSH);
}

int
IWString_and_File_Descriptor::_zstd_compress_and_write(ZSTD_EndDirective mode) {
  ZSTD_inBuffer input = {_things, static_cast<size_t>(_number_elements), 0};
  cerr << "_zstd_compress_and_write writing " << _number_elements << " bytes\n";

  do {
    ZSTD_outBuffer output = {_zstd_buffer, _zstd_buffer_size, 0};

    cerr << "Before ZSTD_compressStream2 pos " << output.pos << " input.pos " << input.pos << '\n';
    size_t rc = ZSTD_compressStream2(_zstd_cctx, &output, &input, mode);
    cerr << "Compression returned " << rc << '\n';
    if (ZSTD_isError(rc)) {
      cerr << "ZSTD_compressStream2 failed " << ZSTD_getErrorName(rc) << '\n';
      return 0;
    }

    cerr << "output.pos " << output.pos << '\n';
    if (output.pos > 0) {
      ssize_t written = IW_FD_WRITE(_fd, _zstd_buffer, output.pos);
      cerr << "WRote " << written << " bytes\n";
      if (written != static_cast<ssize_t>(output.pos)) {
        cerr << "IWString_and_File_Descriptor::_zstd_compress_and_write:write failed\n";
        return 0;
      }
    }

    if (mode == ZSTD_e_end && rc == 0) {
      break;
    }
  } while (input.pos < input.size || mode == ZSTD_e_end);

  cerr << "Before resize " << _number_elements << " characters\n";
  IWString::resize_keep_storage(0);

  return 1;
}

int
IWString_and_File_Descriptor::write_whole_blocks_shift_unwritten() {
  if (0 == IWString::length()) {
    return 0;
  }

  if (_gzfile) {
    if (!_compress_and_write()) {
      return 0;
    }

    return _do_resize(IWSFD_KEEP_STORAGE_ON_FLUSH);
  }

  if (_zstd_buffer != nullptr) {
    if (! _zstd_compress_and_write(ZSTD_e_end)) {
      return 0;
    }

    return 1;
  }

  if (_fd < 0) {
    cerr << "IWString_and_File_Descriptor::write_whole_blocks_shift_unwritten:invalid "
            "file descriptor\n";
    return 0;
  }

  return IWString::write_whole_blocks_shift_unwritten(_fd);
}

int
IWString_and_File_Descriptor::close() {
  if (_gzfile) {
    _compress_and_write();
    gzclose(_gzfile);
    _gzfile = nullptr;
    _do_resize(0);
    return 1;
  }

  if (_zstd_cctx) {
    _zstd_compress_and_write(ZSTD_e_end);
    FreeZstdRelated();
    _do_resize(0);
  }

  if (_fd < 0) {
    cerr << "IWString_and_File_Descriptor::close:not open\n";
    return 0;
  }

  if (_number_elements > 0) {
    flush();
  }

  IW_FD_CLOSE(_fd);

  _fd = -1;

  return 1;
}

int
IWString_and_File_Descriptor::_compress_and_write() {
  if (0 == _number_elements) {
    return 1;
  }

  int bytes_written = gzwrite(_gzfile, _things, static_cast<unsigned>(_number_elements));

  if (bytes_written > 0) {
    return bytes_written;
  }

  cerr << "IWString_and_File_Descriptor::_compress_and_write:cannot write\n";
  return 0;
}

ssize_t
IWString_and_File_Descriptor::write(const char* s, size_t nchars) {
  if (nchars + _number_elements < 32768) {
    IWString::strncat(s, nchars);
    return nchars;
  }

  if (!is_open()) {
    cerr << "IWString_and_File_Descriptor::write:not open\n";
    return 0;
  }

  if (!flush(1)) {  // 1 means keep storage
    return 0;
  }

  ssize_t rc = 0;

  // If everything works, we will return the initial number of bytes requested.
  ssize_t initial_nchars = nchars;

  while (nchars) {
    if (nchars < 32768) {
      IWString::strncat(s, nchars);
      return nchars;
    }

    int chars_to_write;
    if (nchars >= 32768) {
      chars_to_write = 32768;
    } else {
      chars_to_write = nchars;
    }

    ssize_t bytes_written = IW_FD_WRITE(_fd, s, chars_to_write);
    if (bytes_written != static_cast<ssize_t>(chars_to_write)) {
      cerr << "IWString_and_File_Descriptor::write:cannot write " << chars_to_write
           << " to " << _fd << '\n';
      return 0;
    }

    rc += bytes_written;

    s += chars_to_write;
    nchars -= chars_to_write;
  }

  if (rc != initial_nchars) {
    cerr << "IWString_and_File_Descriptor::write:wrote " << rc << " of " << initial_nchars
         << " bytes\n";
    return 0;
  }

  return rc;
}

}  // namespace iwstring
