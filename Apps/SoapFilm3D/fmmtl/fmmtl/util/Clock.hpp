#pragma once

#include <chrono>
#include <iostream>

/** An RAII class
 * Updates a listener with the amount of time the Ticker was alive.
 */
template <typename Listener>
class TickerNotifier {
 public:
  typedef std::chrono::high_resolution_clock clock;
  typedef typename clock::time_point         time_point;
  typedef typename clock::duration           duration_type;
  typedef typename duration_type::rep        tick_type;

  // Default constructor
  TickerNotifier()
      : starttime_(clock::now()) {}
  // Listener object constructor
  explicit TickerNotifier(const Listener& owner)
      : owner_(owner), starttime_(clock::now()) {}
  // Forwarding listener object constructor
  template <typename... Args>
  explicit TickerNotifier(Args&&... args)
      : owner_(args...), starttime_(clock::now()) {}
  // Disable copying
  //TickerNotifier(const TickerNotifier&) = delete;
  //TickerNotifier& operator=(const TickerNotifier&) = delete;
  // Destructor
  ~TickerNotifier() {
    owner_(*this);
  }
  // Restart the timer on this Ticker
  void start() {
    starttime_ = clock::now();
  }
  // Get the duration on this Ticker
  duration_type duration() const {
    return clock::now() - starttime_;
  }
  // Get the seconds on this Ticker
  double seconds() const {
    typedef std::chrono::duration<double> units;
    return std::chrono::duration_cast<units>(duration()).count();
  }
 private:
  Listener owner_;
  time_point starttime_;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const TickerNotifier<T>& t) {
  return os << t.seconds();
}

/** A quick class for timing code:
 * Usage:
 *
 * Clock clock;
 * ...
 * clock.start();
 * // code to time
 * double time = clock.seconds();
 */
struct noop {
  template <typename T>
  void operator()(const T&) const {}
};
typedef TickerNotifier<noop> Clock;


/** A class for timing a scope of code:
 * Usage:
 *
 * { ScopeClock sc("Section Y: ");
 * // code to time
 * }
 * // Prints out "Section Y: XX secs"
 */
struct TickerPrinter {
  TickerPrinter(const std::string& _msg = "", std::ostream& _os = std::cout)
      : msg(_msg), os(_os) {}
  void operator()(const TickerNotifier<TickerPrinter>& t) {
    os << msg << t << " secs" << std::endl;
  }
  std::string msg;
  std::ostream& os;
};
typedef TickerNotifier<TickerPrinter> ScopeClock;
