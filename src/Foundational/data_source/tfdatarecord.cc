#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <unistd.h>

#ifdef __APPLE__
#include <machine/endian.h>
#else
#include <endian.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

// Building crc32c with Bazel is difficult.
// For expediency, we can disable the crc based error checking.

#define LILLYMOL_HAS_CRC32C

#ifdef LILLYMOL_HAS_CRC32C
// if using crc32c directly     #include "crc32c/crc32c.h"
#include "absl/crc/crc32c.h"
#endif

#include "Foundational/iwmisc/misc.h"

#include "tfdatarecord.h"

namespace iw_tf_data_record {

using std::cerr;

unsigned int default_read_buffer_size = 4096; 

constexpr uint64_t sizeof_length = sizeof(uint64_t);
constexpr uint64_t sizeof_crc = sizeof(uint32_t);

#ifdef LILLYMOL_HAS_CRC32C
uint32_t
MaskedCrc(uint32_t crc) {
  return ((crc >> 15) | (crc << 17)) + 0xa282ead8ul;
}
#endif

void
TFDataReader::DefaultValues() {
  _fd = -1;
  _good = true;
  _eof = false;
  _next = 0;

  _read_buffer.resize(default_read_buffer_size);

//_compression_type = kUncompressed;

  _items_read = 0;
}

TFDataReader::TFDataReader() {
  DefaultValues();
}

TFDataReader::TFDataReader(const char * fname) {
  DefaultValues();
  OpenFile(fname);
}

TFDataReader::TFDataReader(IWString& fname) {
  DefaultValues();
  OpenFile(fname.null_terminated_chars());
}

TFDataReader::TFDataReader(const const_IWSubstring& fname) {
  DefaultValues();
  IWString tmp(fname);
  OpenFile(tmp.null_terminated_chars());
}

TFDataReader::~TFDataReader() {
  if (_fd >= 0) {
    close(_fd);
  }
}

int
TFDataReader::Open(const char * fname) {
  return OpenFile(fname);
}

int
TFDataReader::Open(IWString& fname) {
  return OpenFile(fname);
}

int
TFDataReader::Open(const const_IWSubstring& fname) {
  IWString tmp(fname);
  return OpenFile(tmp);
}

bool
TFDataReader::OpenFile(const char * fname) {
  if (_fd >= 0) {
    cerr << "TFDataReader::OpenFile:already open " << _fd << ", no action\n";
    return 0;
  }

  _fd = IW_FD_OPEN(fname, O_RDONLY);
#ifdef DEBUG_IW_TF_DATA
  cerr << "Openfile " << fname << " returned " << _fd << '\n';
#endif
  if (_fd < 0) {
    _good = false;
    return false;
  }

  _good = true;
  _eof = false;
  return true;
}

// _next is pointing at length data. Retrieve that length.
std::optional<uint64_t>
TFDataReader::GetLength() {
  // Records are packed without alignment padding. Copy scalar fields rather
  // than dereferencing potentially unaligned uint64_t/uint32_t pointers.
  uint64_t length;
  uint32_t crc;
  std::memcpy(&length, _read_buffer.rawdata() + _next, sizeof(length));
  std::memcpy(&crc, _read_buffer.rawdata() + _next + sizeof_length, sizeof(crc));
#ifdef LILLYMOL_HAS_CRC32C
  const absl::string_view bytes(_read_buffer.rawdata() + _next, sizeof_length);
  const uint32_t expected = MaskedCrc(static_cast<uint32_t>(absl::ComputeCrc32c(bytes)));
  if (expected != crc) {
    cerr << "TFDataReader::GetLength:invalid length CRC\n";
    _good = false;
    return std::nullopt;
  }
#else
  (void) crc;
#endif
  _next += sizeof_length + sizeof_crc;
  return length;
}

std::optional<const_IWSubstring>
TFDataReader::Next() {
  if (_eof || ! _good) {
    cerr << "TFDataReader::Next:eof or not good\n";
    return std::nullopt;
  }

  // EOF is successful only between complete records. A short header or
  // payload must set good() false, even when the underlying read reached EOF.
  while (_read_buffer.size() - _next < sizeof_length + sizeof_crc) {
    if (!FillReadBuffer(sizeof_length + sizeof_crc)) {
      if (_read_buffer.size() != _next) {
        cerr << "TFDataReader::Next:truncated record header\n";
        _good = false;
      }
      return std::nullopt;
    }
  }

  const std::optional<uint64_t> length = GetLength();
  if (!length) {
    return std::nullopt;
  }
  // The backing resizable_array uses int sizes. Reject invalid lengths before
  // arithmetic, allocation or construction of an out-of-bounds substring.
  if (*length > static_cast<uint64_t>(std::numeric_limits<int>::max()) - sizeof_crc) {
    cerr << "TFDataReader::Next:record exceeds buffer size limit\n";
    _good = false;
    return std::nullopt;
  }
  while (_read_buffer.size() - _next < *length + sizeof_crc) {
    if (!FillReadBuffer(*length + sizeof_crc)) {
      cerr << "TFDataReader::Next:truncated record data\n";
      _good = false;
      return std::nullopt;
    }
  }

  const_IWSubstring result(_read_buffer.rawdata() + _next, *length);

#ifdef LILLYMOL_HAS_CRC32C
  // If using crc32c directly
  // const uint32_t crc_data = crc32c::Crc32c(result.data(), *length);
  // const uint32_t masked_crc = MaskedCrc(crc_data);

  const absl::crc32c_t hash = absl::ComputeCrc32c(absl::string_view(result.data(), *length));
  const uint32_t masked_crc = MaskedCrc(static_cast<uint32_t>(hash));

  uint32_t crc;
  std::memcpy(&crc, _read_buffer.rawdata() + _next + *length, sizeof(crc));
  if (crc != masked_crc) {
    cerr << "TFDataReader::Next:Invalid data crc " << *length << " bytes\n";
    _good = 0;
    return std::nullopt;
  }
#endif

  _next += *length + sizeof_crc;

  ++_items_read;
  return result;
}

// The class needs to be able to read `bytes_needed` into an item.
// Upon exit, the next bytes_needed of data will be in the buffer.
bool
TFDataReader::FillReadBuffer(uint64_t bytes_needed) {
  // If the next bytes_needed bytes are already in _read_buffer we are done.
  if (_next + bytes_needed <= static_cast<uint64_t>(_read_buffer.number_elements())) {
    return true;
  }
  
  // We need to shift the data and maybe resize.
  if (_next > 0) {
    _read_buffer.remove_from_to(0, _next);
    _next = 0;
  }

  // At this stage, the next item will start at zero.

  if (static_cast<uint64_t>(_read_buffer.elements_allocated()) < bytes_needed) {
    _read_buffer.resize(bytes_needed);
  }
  const int bytes_already_present = _read_buffer.number_elements();
  const int to_read = _read_buffer.elements_allocated() - bytes_already_present;

  int bytes_read = IW_FD_READ(_fd, _read_buffer.rawdata() + bytes_already_present, to_read);
#ifdef DEBUG_IW_TF_DATA
  cerr << "Reading " << to_read << " got " << bytes_read << " bytes from file " << _fd << '\n';
#endif
  if (bytes_read < 0) {
    _good = false;
    return false;
  }
  if (bytes_read == 0) {
    _eof = true;
    return false;
  }

  _read_buffer.set_number_elements(bytes_already_present + bytes_read);
  return true;
}

int
TFDataReader::seek_zero() {
  if (_fd < 0) {
    return 0;
  }

  off_t rc = IW_FD_LSEEK(_fd, 0, SEEK_SET);

  if (rc < 0) {
    cerr << "TFDataReader::seek_zero:cannot seek back to start of file\n";
    _good = false;
    return 0;
  }

  _good = true;
  _eof = false;
  _next = 0;
  _read_buffer.resize_keep_storage(0);
  return 1;
}


TFDataWriter::TFDataWriter() {
}

TFDataWriter::TFDataWriter(int fd) : _output(fd) {
}

int
TFDataWriter::Open(IWString& fname) {
  if (_output.is_open()) {
    std::cerr << "TFDataWriter::Open:already open\n";
    return 0;
  }

  if (! _output.open(fname.null_terminated_chars())) {
    return 0;
  }

  return 1;
}

int
TFDataWriter::Open(const char* fname) {
  if (_output.is_open()) {
    std::cerr << "TFDataWriter::Open:already open\n";
    return 0;
  }

  if (! _output.open(fname)) {
    return 0;
  }

  return 1;
}

int
TFDataWriter::Open(const const_IWSubstring& fname) {
  if (_output.is_open()) {
    std::cerr << "TFDataWriter::Open:already open\n";
    return 0;
  }
  IWString tmp(fname);
  if (! _output.open(tmp.null_terminated_chars())) {
    return 0;
  }
  return 1;
}

int
TFDataWriter::Close() {
  return _output.close();
}

ssize_t
TFDataWriter::Write(const void * data, uint64_t nbytes) {
  if (! WriteLength(nbytes)) {
    cerr << "TFDataWriter::Write:cannot write " << nbytes << " length\n";
    return 0;
  }
  return WriteData(data, nbytes);
}

// Write `nbytes` from `value` as well as the crc of `data`.
int
TFDataWriter::CommonWrite(const void * data, const uint64_t nbytes) {
  const char * p = reinterpret_cast<const char *>(data);
  if (nbytes > 0 && ! _output.write(p, nbytes)) {
    cerr << "TFDataWriter::CommonWrite:cannot write " << nbytes << " bytes\n";
    return 0;
  }

#ifdef LILLYMOL_HAS_CRC32C
  // If using crc32c directly.
  // uint32_t crc = crc32c::Crc32c(p, nbytes);
  // crc = MaskedCrc(crc);

  const absl::crc32c_t hash = absl::ComputeCrc32c(absl::string_view(p, nbytes));
  const uint32_t crc = MaskedCrc(static_cast<uint32_t>(hash));
#else
  uint32_t crc = 0;
#endif
  p = reinterpret_cast<const char *>(&crc);
  if (! _output.write(p, sizeof_crc)) {
    cerr << "TFDataWriter::CommonWrite:cannot write crc, nbytes " << nbytes << '\n';
    return 0;
  }

  return _output.write_if_buffer_holds_more_than(8192);
}

ssize_t
TFDataWriter::WriteLength(const uint64_t nbytes) {
  const char * p = reinterpret_cast<const char *>(&nbytes);
  return CommonWrite(p, sizeof_length);
}

ssize_t
TFDataWriter::WriteData(const void * data, uint64_t nbytes) {
  // Even an empty payload has a data CRC in TFRecord framing.
  return CommonWrite(data, nbytes);
}

}  // namespace iw_tf_data_record
