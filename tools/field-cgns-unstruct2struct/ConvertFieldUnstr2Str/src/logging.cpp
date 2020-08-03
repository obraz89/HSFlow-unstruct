#define _SCL_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS

///////////////////////////////////////////////////////////////////////////////
// Name:        logging.cpp
// Purpose:     Centralized logging with MPI synchronization
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#include <queue>

#include "logging.h"

#include <stdarg.h>

/**
 * Make std::string with printf-style formatting
 *
 * @param fmt format string like "%s"
 * @return
 */
std::string hs_string_format(const char* fmt, ...)
{
	va_list args;

	static int len_max = 0;
	static char* buf = NULL;

	// Try to format into existing buffer
	va_start(args, fmt);
	int len = vsnprintf(buf, len_max, fmt, args);
	va_end(args);

	if (len >= len_max)
	{
		// Increase buffer
		len_max = len + 1;
		delete[] buf;   buf = new char[len_max];

		// Format again
		va_start(args, fmt);
		vsnprintf(buf, len_max, fmt, args);
		va_end(args);
	}

	return buf;
}

using hsflow::TLog;
using hsflow::ELogLevel;

#if defined(NDEBUG)
static const int MSG_MAX_LEN = 512;
#else
static const int MSG_MAX_LEN = 1024;
#endif

TLog::Impl* TLog::_instance = nullptr;
//-----------------------------------------------------------------------------

/**
 *  TLog::Impl singleton class implementation
 *  for logging to both stdout and file-stream
 */
class TLog::Impl
{
private:
	// log file descriptor (may be NULL)
	FILE* _file = nullptr;

	// Messages queue local to the current MPI rank
	std::queue<std::string> _q_msgs;

	// Put message to the local queue or print
	void log_string(const char* szString);

public:
	Impl() = default;
	~Impl();

	bool set_file(const std::string& fn);
	void log(ELogLevel level, const std::string& msg);
	void log_raw_from_root(const std::string& msg, bool to_file_only);

	// Collects messages received from several MPI ranks
	void sync();

	// Flushes messages to stdout and/or logfile
	void flush();
};
//-----------------------------------------------------------------------------


TLog::Impl* TLog::instance()
{
	if (!_instance)
		_instance = new TLog::Impl();
	return _instance;
}

bool TLog::set_file(const std::string& fn)
{
	return instance()->set_file(fn);
}

void TLog::log(ELogLevel level, const std::string& msg)
{
	instance()->log(level, msg);
}

void TLog::log_raw_from_root(const std::string& msg, bool to_file_only)
{
	instance()->log_raw_from_root(msg, to_file_only);
}

void TLog::flush()
{
	instance()->sync();
}

void TLog::destroy()
{
	delete _instance;
	_instance = nullptr;
}
//-----------------------------------------------------------------------------


TLog::Impl::~Impl()
{
	sync();
	if (_file)
		fclose(_file);
}
//-----------------------------------------------------------------------------

/**
 * Opens text file to log into in addition to stdout
 *
 * @param[in] fn file name
 *
 * @return true if succeded, false otherwise
 */
bool TLog::Impl::set_file(const std::string& fn)
{
	if (_file) {
		sync();
		fclose(_file);
	}

	_file = fopen(fn.c_str(), "at");
	return _file;
}
//-----------------------------------------------------------------------------

/**
 * Decorates a message and passes it for logging
 *
 * @param[in] level logging type
 * @param[in] msg   original message text
 */
void TLog::Impl::log(ELogLevel level, const std::string& msg)
{
	static std::string s;

	// Decorate the message based on the log level
	s.clear();
	bool from_all_ranks = true;
	switch (level)
	{
	case hsLOG_ERROR:   s = "Error: " + msg;  break;
	case hsLOG_MESSAGE: s = msg;  from_all_ranks = false;  break;
	case hsLOG_WARNING: s = "Warning: " + msg;  break;
	case hsLOG_DEBUG:
#if defined(NDEBUG)
		return;
#endif
		s = "Debug: " + msg;  break;
	default:            s = msg;
	}

	log_string(s.c_str());
}
//-----------------------------------------------------------------------------

/**
*  Put string to the message queue of the current MPI rank
*  NB: Flush() must be called ASAP on *each* MPI rank (i.e. collectively)
*
*  @param[in] szString - message to log
*/
void TLog::Impl::log_string(const char* szString)
{
		// Log immediately from root rank and ignore others
			FILE* streams[] = { stdout, _file, nullptr };

			for (FILE** f = streams; *f; ++f) {
				fprintf(*f, "%s\n", szString);
				fflush(*f);
			}
}
//-----------------------------------------------------------------------------


void TLog::Impl::log_raw_from_root(const std::string& msg, bool to_file_only)
{
	{
		if (!to_file_only)
			printf("%s\n", msg.c_str());

		if (_file)
			fprintf(_file, "%s\n", msg.c_str());
	}
}
//-----------------------------------------------------------------------------


/**
 *  Flushes buffered messages to stdout and log-file
 */
void TLog::Impl::flush()
{
	FILE* streams[] = { stdout, _file, nullptr };

	for (FILE** f = streams; *f; ++f) {
		fflush(*f);
	}
}
//-----------------------------------------------------------------------------

void TLog::Impl::sync()
{
	FILE* streams[] = { stdout, _file, nullptr };

	if (true)
	{
		// MPI is not active -> print messages by itself
		while (!_q_msgs.empty())
		{
			const char* msg = _q_msgs.front().c_str();
			for (FILE** f = streams; *f; ++f) {
				fprintf(*f, "[%d] %s\n", 0, msg);
			}
			_q_msgs.pop();
		}
		flush();
		return;
	}
}
//-----------------------------------------------------------------------------


/**
 *  Log error with source file name
 *
 * @param[in] msg  text message to log
 * @param[in] src  source code file name the message comes from
 * @param[in] line line number in the source file where the message is generated
 */
void hs_log_error_src(const std::string& msg, const char* src, int line)
{
	const std::string& s = hs_string_format("%s\n\t@ %s(%d)", msg.c_str(), src, line);
	hsflow::TLog::log(hsflow::hsLOG_ERROR, s);
}
//-----------------------------------------------------------------------------
