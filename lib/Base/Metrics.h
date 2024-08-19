#pragma once
#include <algorithm>
#include <cmath>
#include <functional>

#include "Base/Base.h"
#include "Logger.h"

double QError(double y_true, double y_measured) {
  return std::max(y_measured, 1.0) / std::max(y_true, 1.0);
}

double logQError(double y_true, double y_measured) {
  return log10(QError(y_true, y_measured));
}
