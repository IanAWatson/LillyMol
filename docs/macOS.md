# MacOS

The Mac installation is relatively recent and not completely tested. If you encounter
problems, please raise an issue.

Several software packages are needed.

- install [homebrew](https://brew.sh)
- `brew install bazelisk libomp wget gnutls xz bash gcc`
- `brew install ruby protobuf`
- `gem install google-protobuf`
- clone repo
- `cd /path/to/LillyMol`
- `export LILLYMOL_HOME=/path/to/LillyMol`
- `make -j 10 2>&1 | tee make.log`

bash is needed because the default system bash is too old, but the
version installed by brew is new enough.

Currently the macOS build depends on gcc rather than clang. Work is
underway to address that.

I wish to express my thanks to Harry Stern for providing invaluable help with the
LillyMol MacOS port.

All Bazel unit tests should pass. If you run the tests in the test directory
you may find some "failing" simply due to different orderings of C++
hashes between compiler versions and machines. Eventually we will
address these discrepancies - they are not failures, but poorly
designed tests.
