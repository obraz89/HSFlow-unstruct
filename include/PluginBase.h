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

#include <map>

#include "logging.h"
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

		class Lookup : public std::map<K, size_t> {
			// nothing to add to the standard map
		};
		std::unique_ptr<Lookup> lookup_;  // lookup table

	public:
		OrderedMap() : lookup_{ new Lookup() }
		{
			;
		};
		OrderedMap(const OrderedMap<K, T>& that)
			: data_{ that.data_ }, lookup_{ new Lookup() }
		{
			*lookup_ = *that.lookup_;
		}
		OrderedMap& operator=(const OrderedMap&) = delete;
		~OrderedMap() = default;

		using iterator = typename std::vector<Value>::iterator;
		iterator begin() { return data_.begin(); }
		iterator end() { return data_.end(); }

		bool insert(const Value& key_val) {
			if (has(key_val.first))
				return false;

			data_.push_back(key_val);
			(*lookup_)[key_val.first] = data_.size() - 1;
			return true;
		};
		bool emplace(const K& key, const T& val) {
			if (has(key))
				return false;

			data_.emplace_back(key, val);
			(*lookup_)[key] = data_.size() - 1;
			return true;
		};

		// throws an std::out_of_range if key not exists
		T& at(const K& key) {
			return data_[lookup_->at(key)].second;
		};
		const T& at(const K& key) const {
			return data_[lookup_->at(key)].second;
		};

		bool empty() const {
			return data_.empty();
		}

		bool has(const K& k) const {
			return lookup_->find(k) != lookup_->end();
		};

		void clear() {
			data_.clear();  lookup_->clear();
		};

		std::vector<K> get_list_of_keys() const {
			std::vector<K> keys;
			typename std::vector<Value>::const_iterator it;
			for (it = data_.begin(); it < data_.end(); it++)
				keys.push_back(it->first);
			return keys;
		};
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

