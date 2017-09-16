#pragma once

#include <iostream>
#include <iomanip>

#include <map>
#include <string>

#include <atomic>

#include "fmmtl/util/Clock.hpp"

#include "fmmtl/config.hpp"

/** A class for accumulating the time of segments of code:
 * Usage:
 *
 * Timer timer;
 * { auto ts = timer.time_scope();
 *   // code to time
 * }
 * std::cout << timer << std::endl;
 */
class Timer {
 public:
  typedef TickerNotifier<Timer&>       ticker_type;
  typedef ticker_type::clock           clock;
  typedef typename clock::time_point   time_point;
  typedef typename clock::duration     duration_type;
  typedef typename duration_type::rep  tick_type;

  // Start by returning an RAII ticker
  ticker_type time_scope() {
    return ticker_type(*this);
  }
  // Add a ticker to the duration
  void operator()(const ticker_type& t) {
    ticks_ += t.duration().count();
  }
  void operator+=(const duration_type& t) {
    ticks_ += t.count();
  }
  // Reset this timer
  void reset() {
    ticks_ = tick_type(0);
  }
  // Get the duration on this Timer
  duration_type duration() const {
    return duration_type(ticks_);
  }
  // Get the seconds on this Timer
  double seconds() const {
    typedef std::chrono::duration<double> units;
    return std::chrono::duration_cast<units>(duration()).count();
  }
  // Print this Timer
  friend std::ostream& operator<<(std::ostream& s, const Timer& t) {
    return s << t.seconds() << "secs";
  }
 private:
  tick_type ticks_;
};


/** @class Logger
 * @brief Logging class to keep the hit count and total time of sections of code
 *
 * Usage:
 * Logger logger;
 *
 * { Logger::timer timer = logger.log("identifier");
 *  // code to track
 * }
 *
 * // Print all string-identified events' hit counts and timings
 * std::cout << log << std::endl;
 */
class Logger {
  struct EventData;

 public:
  typedef TickerNotifier<EventData&>  ticker_type;
  typedef ticker_type::clock          clock;
  typedef typename clock::time_point  time_point;
  typedef typename clock::duration    duration;
  typedef typename duration::rep      tick_type;

  /** Start a ticker for an event. */
  inline ticker_type log(const std::string& event) {
    auto range = data_.equal_range(event);
    if (range.first == range.second) {
#pragma omp critical
      {
      range.first = data_.insert(range.first,
                                 std::make_pair(event, EventData()));
      }
    }
    return ticker_type((*range.first).second);
  }

  /** Erase all events in timer */
  inline void clear() {
#pragma omp critical
    {
      data_.clear();
    }
  }

  //! Print all events and timing to an ostream
  friend std::ostream& operator<<(std::ostream& s, const Logger& log) {
    for (auto& event : log.data_)
      s << std::setw(20) << std::left << event.first
        << ": " << event.second << std::endl;
    return s;
  }

 private:
  struct EventData {
    void operator()(const ticker_type& ticker) {
      total_ += ticker.duration();
      ++call_;
    }
    double seconds() const {
      return total_.seconds();
    }
    unsigned calls() const {
      return call_;
    }
    friend std::ostream& operator<<(std::ostream& s, const EventData& e) {
      return s << e.calls() << " (calls) * "
               << e.seconds() / e.calls() << " (sec/call) = "
               << e.seconds() << " (sec)";
    }
   private:
    Timer total_;
    unsigned call_;
  };

  // A map of string identifiers to EventData (total_time, #calls)
  std::map<std::string, EventData> data_;
};


#if defined(FMMTL_LOGGING)

//! Global static logger rather than a singleton for efficiency/consistency
static Logger fmmtl_logger;

#  define FMMTL_LOG(STRING) auto t##__LINE__ = fmmtl_logger.log(std::string(STRING) + " [" + std::to_string(omp_get_thread_num()) + ']')
#  define FMMTL_PRINT_LOG(OUT) OUT << fmmtl_logger << std::endl
#  define FMMTL_LOG_CLEAR fmmtl_logger.clear()

#else

#  define FMMTL_LOG(STRING)
#  define FMMTL_PRINT_LOG(OUT)
#  define FMMTL_LOG_CLEAR

#endif
