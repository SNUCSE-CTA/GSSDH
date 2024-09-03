#pragma once
#include <any>
#include <map>

#include "Metrics.h"

std::function<double(double)> identity = [](double x) { return x; };
std::function<double(double)> square = [](double x) { return x * x; };
std::function<double(double)> absolute_value = [](double x) {
  return std::abs(x);
};

std::string to_string(__int128_t num) {
  std::string result;
  bool isNegative = num < 0;
  if (isNegative) {
    num = -num;
  }
  do {
    int digit = num % 10;
    result.push_back(digit + '0');
    num /= 10;
  } while (num != 0);
  if (isNegative) {
    result.push_back('-');
  }
  std::reverse(result.begin(), result.end());
  return result;
}

enum ResultType {
  RESULT_INT,
  RESULT_INT64,
  RESULT_INT128,
  RESULT_DOUBLE_FIXED,
  RESULT_DOUBLE_SCIENTIFIC,
  RESULT_STRING
};

class ResultLogger {
public:
  std::map<std::string, std::any> result;
  std::map<std::string, ResultType> result_types;

  double GetNumericResult(const std::string &key) {
    if (!(result.find(key) != result.end()))
      return 0.0;
    std::any value = result[key];
    switch (result_types[key]) {
    case RESULT_INT:
      return 1.0 * std::any_cast<int>(value);
    case RESULT_INT64:
      if (result_types[key] == RESULT_INT64) {
        if (value.type() == typeid(int64_t)) {
          return 1.0 * std::any_cast<int64_t>(value);
        } else if (value.type() == typeid(int)) {
          return 1.0 * std::any_cast<int>(value);
        } else {
          std::cerr << "Unexpected type for key '" << key
                    << "'. Actual type held: " << value.type().name()
                    << std::endl;
        }
      }
    // case RESULT_INT128:
    //     return 1.0 * std::any_cast<__int128_t>(value);
    case RESULT_DOUBLE_FIXED:
      if (value.type() == typeid(int)) {
        return static_cast<double>(std::any_cast<int>(value));
      } else if (value.type() == typeid(double)) {
        return std::any_cast<double>(value);
      } else {
        std::cerr << "Unexpected type for key '" << key
                  << "'. Actual type held: " << value.type().name()
                  << std::endl;
        return 0.0;
      }

    case RESULT_DOUBLE_SCIENTIFIC:
      if (value.type() == typeid(int)) {
        return static_cast<double>(std::any_cast<int>(value));
      } else if (value.type() == typeid(double)) {
        return std::any_cast<double>(value);
      } else {
        std::cerr << "Unexpected type for key '" << key
                  << "'. Actual type held: " << value.type().name()
                  << std::endl;
        return 0.0;
      }

    default:
      return 0.0;
    }
  }

  void PrintResults() {
    for (auto &[key, value] : result) {
      if (!(result.find(key) != result.end()) or
          !(result_types.find(key) != result_types.end()))
        continue;
      switch (result_types[key]) {
      case RESULT_INT:
        fprintf(stdout, "  [Result] %-20s: %d\n", key.c_str(),
                std::any_cast<int>(value));
        break;
      case RESULT_INT64:
        fprintf(stdout, "  [Result] %-20s: %ld\n", key.c_str(),
                std::any_cast<int64_t>(value));
        break;
      // case RESULT_INT128:
      //   fprintf(stdout, "  [Result] %-20s: %s\n", key.c_str(),
      //           to_string(std::any_cast<__int128_t>(value)).c_str());
      //   break;
      case RESULT_DOUBLE_FIXED:
        fprintf(stdout, "  [Result] %-20s: %.04lf\n", key.c_str(),
                std::any_cast<double>(value));
        break;
      case RESULT_DOUBLE_SCIENTIFIC:
        fprintf(stdout, "  [Result] %-20s: %.04e\n", key.c_str(),
                std::any_cast<double>(value));
        break;
      case RESULT_STRING:
        fprintf(stdout, "  [Result] %-20s: %s\n", key.c_str(),
                std::any_cast<std::string>(value).c_str());
        break;
      default:
        break;
      }
    }
  }

  void PrintResults(const std::vector<std::string> &keys) {
    for (auto &key : keys) {
      if (!(result.find(key) != result.end()) or
          !(result_types.find(key) != result_types.end()))
        continue;
      std::any value = result[key];
      switch (result_types[key]) {
      case RESULT_INT:
        fprintf(stdout, "  [Result] %-20s: %d\n", key.c_str(),
                std::any_cast<int>(value));
        break;
      case RESULT_INT64:
        fprintf(stdout, "  [Result] %-20s: %ld\n", key.c_str(),
                std::any_cast<int64_t>(value));
        break;
      // case RESULT_INT128:
      //   fprintf(stdout, "  [Result] %-20s: %s\n", key.c_str(),
      //           to_string(std::any_cast<__int128_t>(value)).c_str());
      //   break;
      case RESULT_DOUBLE_FIXED:
        fprintf(stdout, "  [Result] %-20s: %.04lf\n", key.c_str(),
                std::any_cast<double>(value));
        break;
      case RESULT_DOUBLE_SCIENTIFIC:
        fprintf(stdout, "  [Result] %-20s: %.04e\n", key.c_str(),
                std::any_cast<double>(value));
        break;
      case RESULT_STRING:
        fprintf(stdout, "  [Result] %-20s: %s\n", key.c_str(),
                std::any_cast<std::string>(value).c_str());
        break;
      default:
        break;
      }
    }
  }

  void AddResult(const std::string &key, std::any value,
                 ResultType result_type) {
    result[key] = value;
    result_types[key] = result_type;
  }

  void clear() {
    result.clear();
    result_types.clear();
  }

  void CombineResult(ResultLogger from) {
    for (auto &[key, value] : from.result) {
      result[key] = value;
      result_types[key] = from.result_types[key];
    }
  }
};

void CombineResult(ResultLogger &to, ResultLogger &from) {
  for (auto &[key, value] : from.result) {
    to.result[key] = value;
    to.result_types[key] = from.result_types[key];
  }
}

// Aggregating multiple ResultLoggers: Use for combining per-query results.
double Total(std::vector<ResultLogger> &results, const std::string &which,
             const std::function<double(double)> &func = identity) {
  double total = 0.0;
  for (auto it : results) {
    total += func(it.GetNumericResult(which));
  }
  return total;
}

double Average(std::vector<ResultLogger> &results, const std::string &which,
               const std::function<double(double)> &func = identity) {
  double total = Total(results, which, func);
  return (total / (results.size() * 1.0));
}

double Std(std::vector<ResultLogger> &results, const std::string &which,
           const std::function<double(double)> &func = identity) {
  double total = Total(results, which, identity);
  double sqtotal = Total(results, which, square);
  return sqrt((sqtotal - total * total) / (results.size() * 1.0));
}
