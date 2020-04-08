#pragma once

///////////////////////////////////////////////////////////////////////////////
// Name:        PluginBase.h
// Purpose:     Base classes for plugins
// Author:      Andrey V. Novikov
// Modified:    A. Obraz
///////////////////////////////////////////////////////////////////////////////

#pragma once

#pragma warning(push)
// Disable Visual Studio warnings
#pragma warning(disable:4251)  // class 'std::string' needs to have dll-interface
//-----------------------------------------------------------------------------

#include <memory>
#include <vector>

#include "logging.h"
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


	// ----------------------------------------------------------------------------
	// class TPluginParam
	//
	// Plugin parameter
	// ----------------------------------------------------------------------------
	class TPluginParam
	{
	private:
		std::string _value;

	public:
		enum { ptInt = 0, ptDouble, ptString } type;
		std::string description;

	public:
		bool set_raw_value(const std::string& aVal);

		bool set_value(int val);
		bool set_value(double val);
		bool set_value(const std::string& aVal);

		//--------

		const std::string& get_raw_value() const {
			return _value;
		}

		bool get_value(int& val) const;
		bool get_value(double& val) const;
		bool get_value(const std::string*& val) const;

		//--------

		TPluginParam() : type(ptString) { ; }
		TPluginParam(int val, const std::string& aDescr)
			: type(ptInt), description(aDescr)
		{
			set_value(val);
		}
		TPluginParam(double val, const std::string& aDescr)
			: type(ptDouble), description(aDescr)
		{
			set_value(val);
		}
		TPluginParam(const std::string& val, const std::string& aDescr)
			: type(ptString), description(aDescr)
		{
			set_value(val);
		}
	};
	//-----------------------------------------------------------------------------


	/**
	 * OrderedMap container
	 * A map which preserves the order of insertion.
	 * Minimal subset of methods sutable for the project
	 *
	 * pImpl idiom is employed to hide lookup table implementation
	 *
	 * @tparam K - type of the keys
	 * @tparam T - type of the values
	 */
	template<typename K, typename T>
	class OrderedMap
	{
		using Value = typename std::pair<K, T>;
		std::vector<Value> data_;    // data in insertion order

		class Lookup;
		std::unique_ptr<Lookup> lookup_;  // lookup table

	public:
		OrderedMap();
		OrderedMap(const OrderedMap&);
		OrderedMap& operator=(const OrderedMap&) = delete;
		~OrderedMap();

		using iterator = typename std::vector<Value>::iterator;
		iterator begin() { return data_.begin(); }
		iterator end() { return data_.end(); }

		bool insert(const Value& key_val);
		bool emplace(const K& key, const T& val);

		// throws an std::out_of_range if key not exists
		T& at(const K& key);
		const T& at(const K& key) const;

		bool empty() const {
			return data_.empty();
		}

		bool has(const K& k) const;

		void clear();

		std::vector<K> get_list_of_keys() const;
	};
	//----------------------------------------------------------------------------

	/**
	 * class TPluginParamsGroup
	 *
	 * Logical group of plugin parameters
	 */
	class TPluginParamsGroup
	{
		friend class TPlugin;

	private:
		std::string m_name;
		std::string m_description;
		OrderedMap<std::string, TPluginParam> mapParams;

		const TPluginParam& get_raw_param(const std::string& parName) const;
		void warn_if_exist(const char* szParName)
		{
#if ! defined(NDEBUG)
			if (has_param(szParName))
			{
				hsLogWarning("Parameter '%s' already exists in the group '%s'.",
					szParName, m_name.c_str()
				);
			}
#endif
		}

	public:
		TPluginParamsGroup();
		TPluginParamsGroup(const char* aName, const char aDescr[]);
		~TPluginParamsGroup();

		const std::string& get_name() const { return m_name; }

		bool has_param(const char* szParName) const {
			return mapParams.has(szParName);
		}

		void add(const char* szParName, double value, const char aDescr[]) {
			warn_if_exist(szParName);
			mapParams.emplace(szParName, TPluginParam(value, aDescr));
		}
		void add(const char* szParName, int value, const char aDescr[]) {
			warn_if_exist(szParName);
			mapParams.emplace(szParName, TPluginParam(value, aDescr));
		}
		void add(const char* szParName, const std::string& value, const char aDescr[]) {
			warn_if_exist(szParName);
			mapParams.emplace(szParName, TPluginParam(value, aDescr));
		}

		double get_real_param(const char* pszName) const;
		int get_int_param(const char* pszName) const;
		const std::string& get_string_param(const char* pszName) const;

		std::vector<std::string> get_list_of_names() const;
	};
	//-----------------------------------------------------------------------------


	/**
	 * class TPlugin
	 *
	 * A base class for all plugin classes.
	 */
	class TPlugin
	{
	protected:
		std::string _spec;  // specialization
		OrderedMap<std::string, TPluginParamsGroup> _mapParamsGrps;

		virtual void default_settings() {
			_mapParamsGrps.clear();
		}

	public:
		static bool load_settings(const std::string& file, std::string& ini_data);
		static void save_settings(const std::string& file, const std::string& ini_data);

		TPlugin();
		virtual ~TPlugin();

		virtual std::string get_name() const = 0;
		virtual std::string get_description() const = 0;

		virtual void init(std::string& ini_data, const std::string& spec);

		TPluginParamsGroup& get_settings_grp(const char* szGrpName);

		void parse_settings(std::string& ini_data);

	};
	//-----------------------------------------------------------------------------


#pragma warning(pop)

