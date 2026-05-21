#ifndef FOUNDATIONAL_IWSTRING_IWSTRING_AND_FILE_DESCRIPTOR_H_
#define FOUNDATIONAL_IWSTRING_IWSTRING_AND_FILE_DESCRIPTOR_H_

#include "zlib.h"
#include "zstd.h"

#include "Foundational/iwstring/iwstring.h"

namespace iwstring {

// A File descriptor based output

//  By default, the flush operator resizes the IWString object down to 0.
//  If we want it to keep its allocated storage, pass flush() this arg

#define IWSFD_KEEP_STORAGE_ON_FLUSH -7232

class IWString_and_File_Descriptor : public IWString
{
  private:
    int _fd;

//  If the file name ends in .gz, we gzip our output

    gzFile _gzfile;

    // If the file name ends with .zstd
    ZSTD_CCtx* _zstd_cctx = nullptr;
    char* _zstd_buffer = nullptr;
    size_t _zstd_buffer_size = 0;
    int _zstd_compression_level = 3;

//  private functions

    int _compress_and_write();
    int _do_resize(int);
    int _open_gzipd_stream(const char * fname);
    int _open_zstd_stream(const char* fname);
    int _zstd_compress_and_write(ZSTD_EndDirective mode);
    void FreeZstdRelated();

  public:
    IWString_and_File_Descriptor();
    IWString_and_File_Descriptor(int);
    ~IWString_and_File_Descriptor();

    int write_whole_blocks_shift_unwritten();
    int write_if_buffer_holds_more_than (int s) { if (_number_elements > s) write_whole_blocks_shift_unwritten(); return 1;}
    int flush(int keep_storage = 0);   // by default, we size our buffer back to 0

    int good () const { return (-1 != _fd && IWString::ok());}  // failed open sets _fd to -1

    ssize_t write(const char * s, size_t nchars);

    int fd() const { return _fd;}

    int active () const { return _fd > 0 || nullptr != _gzfile;}
    int is_open() const { return _fd > 0 || nullptr != _gzfile;}

    // If fname starts with '>>' then the file is opened for append.
    int open(const char* fname);
    int close();
};

//extern IWString_and_File_Descriptor & operator<<(int i, IWString_and_File_Descriptor & os) { return os.append_int(i);}
//extern IWString_and_File_Descriptor & operator<<(unsigned int i, IWString_and_File_Descriptor & os) { return os.append_unsigned_int(i);}
}  // name iwstring
#endif // FOUNDATIONAL_IWSTRING_IWSTRING_AND_FILE_DESCRIPTOR_H_
