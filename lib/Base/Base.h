#pragma once
#include <any>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <vector>

#include "Logger.h"
using std::deque;
using std::string;

FILE *log_to = stderr;
// FILE *log_to = fopen("/dev/null","w");
/**
 * @brief String parsing with specified delimeter
 * @Source Folklore
 */

deque<string> parse(string line, const string &del) {
  deque<string> ret;

  size_t pos = 0;
  string token;
  while ((pos = line.find(del)) != string::npos) {
    token = line.substr(0, pos);
    ret.push_back(token);
    line.erase(0, pos + del.length());
  }
  ret.push_back(line);
  return ret;
}

static std::streampos fileSize(const char *filePath) {
  std::streampos fsize = 0;
  std::ifstream file(filePath, std::ios::binary);

  fsize = file.tellg();
  file.seekg(0, std::ios::end);
  fsize = file.tellg() - fsize;
  file.close();

  return fsize;
}

bool CreateDirectory(const std::string &dirName) {
  std::error_code err;
  if (!std::filesystem::create_directories(dirName, err)) {
    if (std::filesystem::exists(dirName)) {
      return true;
    }
    printf("CREATEDIR: FAILED to create [%s], err:%s\n", dirName.c_str(),
           err.message().c_str());
    return false;
  }
  return true;
}

namespace std {
// from boost (functional/hash):
// see http://www.boost.org/doc/libs/1_35_0/doc/html/hash/combine.html template
template <class T> inline void combine(std::size_t &seed, T const &v) {
  seed ^= hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <> struct hash<std::pair<int, int>> {
  auto operator()(const std::pair<int, int> &x) const -> size_t {
    std::size_t seed = 17;
    combine(seed, x.first);
    combine(seed, x.second);
    return seed;
  }
};

template <> struct hash<std::vector<int>> {
  auto operator()(const std::vector<int> &x) const -> size_t {
    std::size_t seed = 17;
    for (auto it : x)
      combine(seed, it);
    return seed;
  }
};

} // namespace std

template <typename T> void EraseIndex(std::vector<T> &vec, int &idx) {
  vec[idx] = vec.back();
  vec.pop_back();
  --idx;
}

__int128_t atoi128_t(const char *str) {
  __int128_t result = 0;
  __int128_t sign = 1;

  if (*str == '-') {
    sign = -1;
    ++str;
  }

  for (; *str != '\0'; ++str) {
    if (*str < '0' || *str > '9') {
      std::cerr << "Invalid input: Non-digit character encountered."
                << std::endl;
      return 0;
    }
    result = result * 10 + *str - '0';
  }

  return result * sign;
}
