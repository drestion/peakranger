#ifndef ASSERT_HELPERS_H_
#define ASSERT_HELPERS_H_

#include <stdexcept>
#include <string>
#include <cassert>
#include <iostream>
/**
 * Assertion for release-enabled assertions
 */
class ReleaseAssertException:public std::runtime_error{
public:
    ReleaseAssertException(const std::string& msg =
    "Runtime error. Run the debug version may potentially figure out the reason") :
    std::runtime_error(msg) {
    }
};

class BadChrException:public std::runtime_error{
public:
    BadChrException(const std::string& msg = "") :
    std::runtime_error(msg) {
    }
};

class GeneralException:public std::runtime_error{
public:
    GeneralException(const std::string& msg = "") :
    std::runtime_error(msg) {
    }
};

class BadFileException:public std::runtime_error{
public:
    BadFileException(const std::string& msg = "") :
    std::runtime_error(msg) {
    }
};

class BadTextLineException:public std::runtime_error{
public:
    BadTextLineException(const std::string& msg = "") :
    std::runtime_error(msg) {
    }
};

class ReadOutOfRegionException:public std::runtime_error{
public:
    ReadOutOfRegionException(const std::string& msg = "") :
    std::runtime_error(msg) {
    }
};

/**
 * Macros for release-enabled assertions, and helper macros to make
 * all assertion error messages more helpful.
 */

#define BADFILE_EXIT(filename) \
		{string error("Fatal:  Can not open "); \
		error += filename; \
		perror(error.c_str());\
		abort();}

#define BADFILE_EXCEPTION(filename) \
		{string error("Fatal:  Can not open "); \
		error += filename; \
		perror(error.c_str());\
		throw BadFileException();}

#define VERB(msg) \
		{if(verbose) cout << msg << endl;}

#define VERB_PRINT_MAP() \
		{if(verbose) for(map<string, int>::iterator it = chr_ind.begin(); \
		 it != chr_ind.end(); it++) cout << it->first << " -> " << it->second  \
		 << endl;}

#define PRINT_RESULT(r) \
		{cout << r.chr <<"\t"<< r.start <<"\t" << r.orientation \
          << "\t" << r.seq << endl;}

#ifndef NDEBUG
#define DEBUG_PRINT(x,msg) \
		std::cout << "DEBUG_PRINT: " << __FILE__ << ":" << __LINE__ << std::endl; \
		std::cout <<  msg << " : " << x << std::endl;
#else
#define DEBUG_PRINT(x,msg)
#endif

#define rt_assert(b)  \
	if(!(b)) { \
		std::cout << "rt_assert at " << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_msg(b,msg)  \
	if(!(b)) { \
		std::cout << msg <<  " at " << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#define rt_assert_eq(ex,ac)  \
	if(!((ex) == (ac))) { \
		std::cout << "rt_assert_eq: expected (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_eq_msg(ex,ac,msg)  \
	if(!((ex) == (ac))) { \
		std::cout << "rt_assert_eq: " << msg <<  ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_eq(ex,ac)  \
	if(!((ex) == (ac))) { \
		std::cout << "assert_eq: expected (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_eq_msg(ex,ac,msg)  \
	if(!((ex) == (ac))) { \
		std::cout << "assert_eq: " << msg <<  ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_eq(ex,ac)
#define assert_eq_msg(ex,ac,msg)
#endif

#define rt_assert_neq(ex,ac)  \
	if(!((ex) != (ac))) { \
	    std::cout << "rt_assert_neq: expected " << ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
	            std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_neq_msg(ex,ac,msg)  \
	if(!((ex) != (ac))) { \
		std::cout << "rt_assert_neq: " << msg << ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_neq(ex,ac)  \
	if(!((ex) != (ac))) { \
		std::cout << "assert_neq: expected not (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_neq_msg(ex,ac,msg)  \
	if(!((ex) != (ac))) { \
		std::cout << "assert_neq: " << msg << ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_neq(ex,ac)
#define assert_neq_msg(ex,ac,msg)
#endif

#define rt_assert_gt(a,b) \
	if(!((a) > (b))) { \
		std::cout << "rt_assert_gt: expected (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_gt_msg(a,b,msg) \
	if(!((a) > (b))) { \
		std::cout << "rt_assert_gt: " << msg << ": (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_gt(a,b) \
	if(!((a) > (b))) { \
		std::cout << "assert_gt: expected (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_gt_msg(a,b,msg) \
	if(!((a) > (b))) { \
		std::cout << "assert_gt: " << msg << ": (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_gt(a,b)
#define assert_gt_msg(a,b,msg)
#endif

#define rt_assert_geq(a,b) \
	if(!((a) >= (b))) { \
		std::cout << "rt_assert_geq: expected (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_geq_msg(a,b,msg) \
	if(!((a) >= (b))) { \
		std::cout << "rt_assert_geq: " << msg << ": (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_geq(a,b) \
	if(!((a) >= (b))) { \
		std::cout << "assert_geq: expected (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_geq_msg(a,b,msg) \
	if(!((a) >= (b))) { \
		std::cout << "assert_geq: " << msg << ": (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_geq(a,b)
#define assert_geq_msg(a,b,msg)
#endif

#define rt_assert_lt(a,b) \
	if(!(a < b)) { \
		std::cout << "rt_assert_lt: expected (" << a << ") < (" << b << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_lt_msg(a,b,msg) \
	if(!(a < b)) { \
		std::cout << "rt_assert_lt: " << msg << ": (" << a << ") < (" << b << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_lt(a,b) \
	if(!(a < b)) { \
		std::cout << "assert_lt: expected (" << a << ") < (" << b << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_lt_msg(a,b,msg) \
	if(!(a < b)) { \
		std::cout << "assert_lt: " << msg << ": (" << a << ") < (" << b << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_lt(a,b)
#define assert_lt_msg(a,b,msg)
#endif

#define rt_assert_leq(a,b) \
	if(!((a) <= (b))) { \
		std::cout << "rt_assert_leq: expected (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_leq_msg(a,b,msg) \
	if(!((a) <= (b))) { \
		std::cout << "rt_assert_leq: " << msg << ": (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_leq(a,b) \
	if(!((a) <= (b))) { \
		std::cout << "assert_leq: expected (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_leq_msg(a,b,msg) \
	if(!((a) <= (b))) { \
		std::cout << "assert_leq: " << msg << ": (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_leq(a,b)
#define assert_leq_msg(a,b,msg)
#endif

#endif /*ASSERT_HELPERS_H_*/
