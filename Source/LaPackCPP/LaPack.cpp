#pragma once


#include "LaPack.hpp"

#include <stdarg.h>

#include <cmath> 
#include <limits>      // uses in lamch to inquire numeric parameters
#include <algorithm>   // for min/max

#include "LibF2C.cpp"
#include "Blas.cpp"
#include "XBlas.cpp"
#include "LapackBanded.cpp"

#include "TemplateInstantiations.cpp"
#include "TypedFunctions.cpp"