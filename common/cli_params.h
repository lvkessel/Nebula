#ifndef __CLI_PARAMS_H_
#define __CLI_PARAMS_H_

/**
 * \brief Very simple command-line parameter parser.
 *
 * It can currently accept two types of arguments. Arugments of the form
 * `--key=value` (an arbitrary number of leading `-` is accepted) are called
 * "flags". Arguments that are not of that form are "positional arguments".
 */

#include <string>
#include <vector>
#include <map>
#include <sstream>

class cli_params
{
public:
	/**
	 * \brief Constructor.
	 * Parameters are the same as what you get in the C++ `main()` function.
	 *
	 * This function parses the arguments and stores them in private variables.
	 *
	 * \param argc Number of command-line arguments.
	 * \param argv The command-line arguments themselves.
	 */
	cli_params(int argc, char** argv)
	{
		program_name = argv[0];
		for (int i = 1; i < argc; ++i)
		{
			std::string item = argv[i];
			if (item[0] == '-')
			{
				size_t key_start = item.find_first_not_of('-');
				size_t split = item.find("=", key_start);
				std::string key = item.substr(key_start, split-key_start);
				std::string value = item.substr(split+1);

				if (value.find('=') != std::string::npos)
					throw std::runtime_error("Two '=' found in " + item);

				auto map_it = flags.find(key);
				if (map_it != flags.end())
					throw std::runtime_error("Key '" + key + "' appears twice in options list");

				flags.insert(map_it, { key, value });
			}
			else
			{
				positional.push_back(item);
			}
		}
	}

	/**
	 * \brief Get the program name, i.e. argv[0].
	 */
	std::string const & get_program_name() const
	{
		return program_name;
	}

	/**
	 * \brief Get the positional arguments (in the correct order).
	 *
	 * The "positional" arguments are arguments that don't follow the --key=value format.
	 */
	std::vector<std::string> const & get_positional() const
	{
		return positional;
	}

	/**
	 * \brief Get an optional flag.
	 *
	 * If the flag designated by `key` was found, fill `value` with the correct
	 * value; otherwise, `value` is left untouched.
	 *
	 * \param[in]  key   Flag name
	 * \param[out] value Variable to store the parameter in
	 */
	template<typename T>
	void get_optional_flag(std::string const & key, T& value) const
	{
		auto map_it = flags.find(key);
		if (map_it != flags.end())
		{
			std::stringstream(map_it->second) >> value;
		}
	}

private:
	std::string program_name;
	std::map<std::string, std::string> flags;
	std::vector<std::string> positional;
};

#endif // __CLI_PARAMS_H_
