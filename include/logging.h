///////////////////////////////////////////////////////////////////////////////
// Name:        logging.h
// Purpose:     HSFlow logging with MPI synchronization
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>

//-----------------------------------------------------------------------------

//
// String manipulation
//
std::string hs_string_format(const char* fmt, ...)
#if defined(__GNUC__)
__attribute__((format(printf, 1, 2)))  // ask GCC to type-check against a format string
#endif
;

namespace hsflow
{

	enum ELogLevel
	{
		hsLOG_ERROR,
		hsLOG_WARNING,
		hsLOG_MESSAGE,
		hsLOG_MESSAGE_ALL_RANKS,
		hsLOG_DEBUG,
	};

	/**
	 *  TLog singleton class with pImpl idiom
	 *  for logging to both stdout and file-stream
	 */
	class TLog
	{
	private:
		class Impl;
		static Impl* _instance;

		TLog() = delete;
		TLog(TLog&) = delete;

		static Impl* instance();

	public:
		static bool set_file(const std::string& fn);

		static void log(ELogLevel level, const std::string& msg);
		static void log_raw_from_root(const std::string& msg, bool to_file_only);

		static void flush();
		static void destroy();
	};

} // namespace
//-----------------------------------------------------------------------------

// Print wall-clock time elapsed from the program start
void hsLogWTime(bool isInit = false);
//-----------------------------------------------------------------------------

// Log error with source file name
void hs_log_error_src(const std::string& msg, const char* src, int line);
//-----------------------------------------------------------------------------


/**
 *  Helper logging macros
 */
#define hsLogError(...)  \
	hsflow::TLog::log( hsflow::hsLOG_ERROR, hs_string_format(__VA_ARGS__) )

#define hsLogErrorSrc(...) \
	hs_log_error_src( hs_string_format(__VA_ARGS__), __FILE__, __LINE__ )

#define hsLogWarning(...)  \
	hsflow::TLog::log( hsflow::hsLOG_WARNING, hs_string_format(__VA_ARGS__) )

#define hsLogMessage(...)  \
	hsflow::TLog::log( hsflow::hsLOG_MESSAGE, hs_string_format(__VA_ARGS__) )

#define hsLogMsgAllRanks(...)  \
	hsflow::TLog::log( hsflow::hsLOG_MESSAGE_ALL_RANKS, hs_string_format(__VA_ARGS__) )

#if ! defined(NDEBUG)
#define hsLogDebug(...)  \
		hsflow::TLog::log( hsflow::hsLOG_DEBUG, hs_string_format(__VA_ARGS__) )
#else
#define hsLogDebug(...)  (void(0))
#endif
 //-----------------------------------------------------------------------------

#define hsTHROW(...)  \
	throw TError( hs_string_format(__VA_ARGS__), __FILE__, __LINE__ )

// Use for compile-time checking of format arguments
// #define hsTHROW(...)  printf(__VA_ARGS__)
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

	/**
	 * The TError class
	 *
	 * A class for all HSFlow exceptions
	 */
class TError
{
protected:
	std::string _what;  // description of the error
	std::string _file;  // source file path
	int         _line;  // source file line number

public:
	TError(const std::string& what, const char* src, const int line)
		: _what(what), _file(src), _line(line) {   }

	std::string what() const { return _what; }
	std::string file() const { return _file; }
	int line() const { return _line; }
	std::string what_detailed() const
	{
		return  hs_string_format("%s\n\t@ %s(%d)",
			_what.c_str(), _file.c_str(), _line
		);
	}
};
//-----------------------------------------------------------------------------

 /**
  * Guard for showing messages on exiting a scope (e.g. function)
  *
  * NB: Make sure it's used on each MPI rank simultaneously
  * (e.g. don't use in loops with different limits)
  */
class TLogSyncGuard
{
public:
	void Flush() {
		hsflow::TLog::flush();
	}

	~TLogSyncGuard() {
		Flush();
	}
};
//-----------------------------------------------------------------------------
