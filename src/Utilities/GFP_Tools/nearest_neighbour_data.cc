#include "nearest_neighbour_data.h"

#include <cmath>
#include <iostream>
#include <type_traits>
#include <utility>

#include "Foundational/data_source/tfdatarecord.h"

namespace nearest_neighbour_data {

bool
NearestNeighbourData::SetNeighboursToRetain(int count) {
  if (_loaded || count < 1) {
    std::cerr << "NearestNeighbourData: retention count must be positive and set before "
                 "reading data\n";
    return false;
  }
  _neighbours_to_retain = count;
  return true;
}

template <typename Proto>
bool
NearestNeighbourData::ReadFile(const char* fname) {
  iw_tf_data_record::TFDataReader input(fname);
  if (!input.good()) {
    std::cerr << "NearestNeighbourData: cannot open '" << fname << "'\n";
    return false;
  }
  std::vector<Proto> records;
  std::unordered_map<std::string, size_t> id_to_index;
  while (true) {
    const auto record = input.Next();
    if (!record) {
      if (!input.good() || !input.eof()) {
        std::cerr << "NearestNeighbourData: error reading '" << fname << "'\n";
        return false;
      }
      break;
    }
    Proto source;
    if (!source.ParseFromArray(record->data(), record->length())) {
      std::cerr << "NearestNeighbourData: invalid proto in '" << fname << "'\n";
      return false;
    }
    if (!source.has_name() || source.name().empty()) {
      std::cerr << "NearestNeighbourData: missing query identifier\n";
      return false;
    }
    if (!id_to_index.emplace(source.name(), records.size()).second) {
      std::cerr << "NearestNeighbourData: duplicate query identifier '" << source.name()
                << "'\n";
      return false;
    }

    // Truncating a repeated proto field can retain the discarded allocations.
    // Build a fresh proto so stored memory reflects only retained neighbours.
    Proto retained;
    retained.set_name(source.name());
    if (source.has_smiles()) {
      retained.set_smiles(source.smiles());
    }
    if constexpr (std::is_same_v<Proto, nnbr::NearNeighbours>) {
      if (source.has_cluster()) {
        retained.set_cluster(source.cluster());
      }
      if (source.has_csize()) {
        retained.set_csize(source.csize());
      }
    }
    for (const auto& nbr : source.nbr()) {
      if constexpr (std::is_same_v<Proto, nnbr::NearNeighbours>) {
        if (nbr.id() == source.name()) {
          continue;
        }
      } else if (nbr.id() == records.size()) {
        continue;
      }
      if (!std::isfinite(nbr.dist()) || nbr.dist() < 0.0f) {
        std::cerr << "NearestNeighbourData: invalid distance for '" << source.name()
                  << "'\n";
        return false;
      }
      *retained.add_nbr() = nbr;
      if (retained.nbr_size() == _neighbours_to_retain) {
        break;
      }
    }
    records.push_back(std::move(retained));
  }

  if constexpr (std::is_same_v<Proto, nnbr::NearNeighbours>) {
    _neighbours.swap(records);
  } else {
    // Only the retained, compact records survive the input pass. Resolve
    // forward references now, without rereading potentially large lists.
    std::vector<nnbr::NearNeighbours> neighbours;
    neighbours.reserve(records.size());
    for (const Proto& record : records) {
      nnbr::NearNeighbours retained;
      retained.set_name(record.name());
      if (record.has_smiles()) {
        retained.set_smiles(record.smiles());
      }
      for (const auto& nbr : record.nbr()) {
        if (nbr.id() >= records.size()) {
          std::cerr << "NearestNeighbourData: neighbour index " << nbr.id()
                    << " out of range for '" << record.name() << "' in '"
                    << fname << "'\n";
          return false;
        }
        const Proto& neighbour = records[nbr.id()];
        auto* dest = retained.add_nbr();
        dest->set_id(neighbour.name());
        dest->set_dist(nbr.dist());
        if (neighbour.has_smiles()) {
          dest->set_smi(neighbour.smiles());
        }
      }
      neighbours.push_back(std::move(retained));
    }
    _neighbours.swap(neighbours);
  }
  _id_to_index.swap(id_to_index);
  _loaded = true;
  return true;
}

bool
NearestNeighbourData::Read(const char* fname) {
  return ReadFile<nnbr::NearNeighbours>(fname);
}

bool
NearestNeighbourData::ReadIndexed(const char* fname) {
  return ReadFile<nnbr::NearNeighboursIndices>(fname);
}

const nnbr::NearNeighbours*
NearestNeighbourData::Lookup(const std::string& id) const {
  const auto iter = _id_to_index.find(id);
  if (iter == _id_to_index.end()) {
    return nullptr;
  }
  return &_neighbours[iter->second];
}

std::optional<float>
NearestNeighbourData::NearestDistance(const std::string& id) const {
  const nnbr::NearNeighbours* proto = Lookup(id);
  if (!proto || proto->nbr().empty()) {
    return std::nullopt;
  }
  return proto->nbr(0).dist();
}
}  // namespace nearest_neighbour_data
