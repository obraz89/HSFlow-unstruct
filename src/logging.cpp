#define _SCL_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS

///////////////////////////////////////////////////////////////////////////////
// Name:        logging.cpp
// Purpose:     Centralized logging with MPI synchronization
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#include <queue>
#include <mpi.h>

#include "common_data.h"  // G_State
#include "common_procs.h"
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
	void log_string(const char* szString, bool fromAllRanks);

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
	case hsLOG_MESSAGE_ALL_RANKS: s = "Msg:" + msg;  break;
	case hsLOG_WARNING: s = "Warning: " + msg;  break;
	case hsLOG_DEBUG:
#if defined(NDEBUG)
		return;
#endif
		s = "Debug: " + msg;  break;
	default:            s = msg;
	}

	log_string(s.c_str(), from_all_ranks);
}
//-----------------------------------------------------------------------------

/**
*  Put string to the message queue of the current MPI rank
*  NB: Flush() must be called ASAP on *each* MPI rank (i.e. collectively)
*
*  @param[in] szString - message to log
*  @param[in] fromAllRanks - log from every MPI rank or just from root
*/
void TLog::Impl::log_string(const char* szString, bool fromAllRanks)
{
	if (!fromAllRanks)
	{
		// Log immediately from root rank and ignore others
		if (G_State.mpiRank == 0)
		{
			FILE* streams[] = { stdout, _file, nullptr };

			for (FILE** f = streams; *f; ++f) {
				fprintf(*f, "%s\n", szString);
				fflush(*f);
			}
		}
	}
	else
		_q_msgs.push(szString);
}
//-----------------------------------------------------------------------------


void TLog::Impl::log_raw_from_root(const std::string& msg, bool to_file_only)
{
	if (G_State.mpiRank == 0)
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


/**
 * Collects messages from each MPI rank on the root rank,
 * merges equivalent ones and flushes to the m_stream
 *
 * NB: Must be called on *each* MPI rank (collective operation)
*/

// List of MPI-ranks as bitset
class Tranks
{
	std::vector<bool> _ranks;  // should be terminated by off-bit

public:
	Tranks(int rank) : _ranks(G_State.mpiNProcs + 1, false) {
		_ranks[rank] = true;
	}

	void operator<<(int rank) {
		_ranks[rank] = true;
	}

	/// Convert bitset to intervals string like "[1-3,5,8-9]"
	const char* range() const
	{
		static std::string s;
		s.clear();

		int rs = -1, re = -1;
		for (int i = 0; i < _ranks.size(); ++i)
		{
			if (_ranks[i]) {
				if (rs < 0) rs = i;
				re = i;
			}
			else {
				if (re >= 0) {
					if (!s.empty())  s += ",";
					char num[16];  sprintf(num, (rs == re ? "%d" : "%d-%d"), rs, re);
					s += num;
				}
				rs = -1;  re = -1;
			}
		}

		return s.c_str();
	}
};

void TLog::Impl::sync()
{
	FILE* streams[] = { stdout, _file, nullptr };

	int noMPI = 0;  MPI_Finalized(&noMPI);
	if (noMPI)
	{
		// MPI is not active -> print messages by itself
		while (!_q_msgs.empty())
		{
			const char* msg = _q_msgs.front().c_str();
			for (FILE** f = streams; *f; ++f) {
				fprintf(*f, "[%d] %s\n", G_State.mpiRank, msg);
			}
			_q_msgs.pop();
		}
		flush();
		return;
	}

	static char szMsg[MSG_MAX_LEN];

	const int tag = 'L' + 'O' + 'G';
	const MPI_Comm& comm = MPI_COMM_WORLD;

	if (G_State.mpiRank == 0)
	{
		// Messages with the originating ranks
		using TmapMsgRanks = std::map<std::string, Tranks>;
		static TmapMsgRanks map_msg_ranks;

		// Vector to keep messages order
		static std::vector<TmapMsgRanks::iterator> msgs_order;

		map_msg_ranks.clear();
		msgs_order.clear();

		// Process own messages
		while (!_q_msgs.empty())
		{
			std::pair<TmapMsgRanks::iterator, bool> res =
				map_msg_ranks.emplace(_q_msgs.front(), Tranks(0));

			msgs_order.push_back(res.first);
			_q_msgs.pop();
		}

		// Receive messages from other processors
		for (int p = 1; p < G_State.mpiNProcs; ++p)
		{
			int nMsg = 0;  MPI_Recv(&nMsg, 1, MPI_INT, p, tag, comm, MPI_STATUS_IGNORE);
			for (int m = 0; m < nMsg; ++m)
			{
				MPI_Recv(szMsg, MSG_MAX_LEN, MPI_CHAR, p, tag, comm, MPI_STATUS_IGNORE);
				copy_n("...", 4, szMsg + MSG_MAX_LEN - 4);  // guard for truncated message not fitting the buffer

				std::pair<TmapMsgRanks::iterator, bool> res =
					map_msg_ranks.emplace(szMsg, Tranks(p));

				if (res.second) // new message was inserted
					msgs_order.push_back(res.first);
				else // the same message exists, add new originating rank to it
					res.first->second << p;
			}
		}

		// Print collected messages prefixing each with originating ranks
		for (const auto& mranks : msgs_order)
		{
			for (FILE** f = streams; *f; ++f) {
				fprintf(*f, "[%s] %s\n", mranks->second.range(), mranks->first.c_str());
			}
		}
		flush();
	}
	else  // G_State.mpiRank != 0
	{
		// Send local message queue to root rank
		int nMsg = static_cast<int>(_q_msgs.size());
		MPI_Ssend(&nMsg, 1, MPI_INT, 0, tag, comm);  // synchronous send prevents root flooding
		while (!_q_msgs.empty())
		{
			const std::string& s = _q_msgs.front();
			MPI_Send((void*)s.c_str(), int(s.size() + 1), MPI_CHAR, 0, tag, comm);
			_q_msgs.pop();
		}
	}
}
//-----------------------------------------------------------------------------


/**
 *  Print wall-clock time elapsed from the program start
 *
 *  @param[in] isInit - reset timer to 0
 */
void hsLogWTime(bool isInit)
{
	static double t0 = 0;

	if (isInit)
	{
		t0 = MPI_Wtime();
		return;
	}

	double t = MPI_Wtime() - t0;
	hsLogMessage(" WTime= %.1f", t);
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
