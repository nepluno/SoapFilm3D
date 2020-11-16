// ---------------------------------------------------------
//
//  options.h
//  Tyson Brochu 2008
//  Christopher Batty, Fang Da 2014
//
//  Constants and macro defines
//
// ---------------------------------------------------------

#ifndef COMMONOPTIONS_H
#define COMMONOPTIONS_H

#include <iostream>

// ---------------------------------------------------------
// Global constants
// ---------------------------------------------------------

const double UNINITIALIZED_DOUBLE = 0x0F;
const double BIG_DOUBLE = 1e30;

#ifdef _MSC_VER
#include "BaseTsd.h"
typedef SSIZE_T ssize_t;
#endif

const size_t UNINITIALIZED_SIZE_T = static_cast<size_t>(~0);

#endif
