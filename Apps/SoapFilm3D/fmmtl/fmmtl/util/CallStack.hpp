#pragma once

#include "fmmtl/config.hpp"

#ifndef FMMTL_RELEASE

#include <string>
#include <stack>

namespace fmmtl {

std::stack<std::string> callStack;

void PushCallStack(const std::string& s) {
#ifdef FMMTL_WITH_OPENMP
  if (omp_get_thread_num() != 0)
    return;
#endif
  callStack.push(s);
}

void PopCallStack() {
#ifdef FMMTL_WITH_OPENMP
  if (omp_get_thread_num() != 0)
    return;
#endif
  callStack.pop();
}

void DumpCallStack(std::ostream& os) {
#ifdef FMMTL_WITH_OPENMP
  if( omp_get_thread_num() != 0 )
    return;
#endif
  std::ostringstream msg;
  while (!callStack.empty()) {
    msg << "[" << callStack.size() << "]: " << callStack.top()
        << "\n";
    callStack.pop();
  }
  os << msg.str();
  os.flush();
}

struct CallStackEntry {
  CallStackEntry(const std::string& s) {
    if (!std::uncaught_exception())
      PushCallStack(s);
  }
  ~CallStackEntry() {
    if (!std::uncaught_exception())
      PopCallStack();
  }
};

} // end namespace fmmtl

#endif
