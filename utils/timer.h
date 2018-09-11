#ifndef TIMER_H_
#define TIMER_H_

#include <ctime>
#include <iostream>
#include <iomanip>
#include <sstream>

/*
 * Adpoted from Bowtie's source codes.
 */
namespace utils {
/**
 * Use time() call to keep track of elapsed time between creation and
 * destruction.  If verbose is true, Timer will print a message showing
 * elapsed time to the given output stream upon destruction.
 */
class Timer {
public:
    Timer(std::ostream& out = std::cout, const char *msg = "", bool verbose =
            true) :
            _t(time(0)), _out(out), _msg(msg), _verbose(verbose) {
    }

    /// Optionally print message
    ~Timer() {
        if (_verbose)
            write(_out);
    }

    /// Return elapsed time since Timer object was created
    time_t elapsed() const {
        return time(0) - _t;
    }

    void write(std::ostream& out) {
        time_t passed = elapsed();
        // Print the message supplied at construction time followed
        // by time elapsed formatted HH:MM:SS

        unsigned int hours = (passed / 60) / 60;
        unsigned int minutes = (passed / 60) % 60;
        unsigned int seconds = (passed % 60);
        out << _msg << std::setfill('0') << std::setw(2) << hours << ":"
                << std::setfill('0') << std::setw(2) << minutes << ":"
                << std::setfill('0') << std::setw(2) << seconds << std::endl;
    }

private:
    time_t _t;
    std::ostream& _out;
    const char *_msg;
    bool _verbose;
};

void logDate(std::ostream& os, bool nl = true);
void getDate(std::string& date);
void logTime(std::ostream& os, bool nl = true);

#ifndef NDEBUG
#define TIMER_START(msg) \
		{Timer timer(cout, msg, verbose);
#define TIMER_END() }
#else
#define TIMER_START(msg)
#define TIMER_END()
#endif
}
#endif /*TIMER_H_*/
