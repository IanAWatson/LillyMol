#include <cstdint>
#include <string>
#include <type_traits>

namespace iwstring {
template <typename Int>
std::string with_commas(Int value) {
    static_assert(std::is_integral_v<Int>, "with_commas requires an integer type");

    std::string s = std::to_string(value);
    const std::size_t start = s.starts_with("-") ? 1 : 0;

    int insert_position = static_cast<int>(s.size()) - 3;

    while (insert_position > static_cast<int>(start)) {
        s.insert(static_cast<std::size_t>(insert_position), ",");
        insert_position -= 3;
    }

    return s;
}

template std::string with_commas<int>(int);
template std::string with_commas<uint32_t>(uint32_t);
template std::string with_commas<int64_t>(int64_t);
template std::string with_commas<long long>(long long);
template std::string with_commas<uint64_t>(uint64_t);

}  // namespace iwstring
