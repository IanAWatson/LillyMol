#ifndef FOUNDATIONAL_IWMISC_COMBINATIONS_H
#define FOUNDATIONAL_IWMISC_COMBINATIONS_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>

namespace combinations {

// Class to facilitate exploring combinations of things.
template <typename T>
class Combinations {
  private:
    // The number if items in each slot.
    std::vector<uint32_t> _count;

    bool _first_call;
    bool _exhausted;

  public:
    Combinations(const std::vector<T>& count);

    // Sizes state on the first call. An empty choice set yields no combinations.
    // Once exhausted, returns zero until Reset is called.
    int Next(std::vector<T>& state);

    void Reset();
};

template <typename T>
Combinations<T>::Combinations(const std::vector<T>& count) : _count(count) {
  _first_call = true;
  _exhausted = false;
}

template <typename T>
int
Combinations<T>::Next(std::vector<T>& state) {
  if (_exhausted) {
    return 0;
  }
  if (_first_call) {
    if (std::find(_count.begin(), _count.end(), 0) != _count.end()) {
      _exhausted = true;
      return 0;
    }
    state.assign(_count.size(), 0);
    _first_call = false;
    return 1;
  }

  //std::cerr << " sizeof _count " << _count.size() << " state " << state.size() << '\n';
  for (int i = static_cast<int>(_count.size()) - 1; i >= 0; --i) {
    //std::cerr << " i = " << i << '\n';
    if (state[i] < _count[i] - 1) {
      ++state[i];
      return 1;
    }
    state[i] = 0;
  }

  _exhausted = true;
  return 0;
}

template <typename T>
void
Combinations<T>::Reset() {
  _first_call = true;
  _exhausted = false;
}

};  // namespace combinations

#endif  // FOUNDATIONAL_IWMISC_COMBINATIONS_H
