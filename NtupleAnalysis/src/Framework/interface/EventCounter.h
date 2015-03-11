// -*- c++ -*-
#ifndef Framework_EventCounter_h
#define Framework_EventCounter_h

#include <string>
#include <vector>

class TH1;
class TDirectory;

class EventWeight;

// Represents a single count, possible to construct only via
// EventCounter
class EventCounter;
class Count {
public:
  ~Count();

  void increment();

  friend class EventCounter;

private:
  Count() = delete;
  Count(EventCounter *ec, size_t counterIndex, size_t countIndex);

  EventCounter *fEventCounter;
  size_t fCounterIndex;
  size_t fCountIndex;
};


/// Event counter
class EventCounter {
private:
  class Counter {
  public:
    explicit Counter(const std::string& n);

    // not copyable
    Counter(const Counter&) = delete;
    Counter& operator=(const Counter&) = delete;
    // but movable
    Counter(Counter&&) = default;
    Counter& operator=(Counter&&) = default;

    bool contains(const std::string& l) const;

    size_t insert(const std::string& label);

    void incrementCount(size_t countIndex, double weight);

    void book(TDirectory *dir);
    void bookWeighted(TDirectory *dir);
    void serialize();

    const std::string& getName() const { return name; }

  private:
    std::string name;
    std::vector<std::string> labels;
    std::vector<long int> values;
    std::vector<double> weights;
    std::vector<double> weightsSquared;

    TH1 *counter;
    TH1 *weightedCounter;
  };

public:
  explicit EventCounter(const EventWeight& weight);
  ~EventCounter();

  // non-copyable
  EventCounter(const EventCounter&) = delete;
  EventCounter& operator=(const EventCounter&) = delete;

  Count addCounter(const std::string& name);
  Count addSubCounter(const std::string& subcounterName, const std::string& countName);

  void setOutput(TDirectory *dir);
  void serialize();

private:
  friend class Count;
  void incrementCount(size_t counterIndex, size_t countIndex);

  size_t findOrInsertCounter(const std::string& name);

  const EventWeight& fWeight;
  std::vector<Counter> fCounters; // main counter has index 0
};

inline
void Count::increment() {
  fEventCounter->incrementCount(fCounterIndex, fCountIndex);
}

#endif
