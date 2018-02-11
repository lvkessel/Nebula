#ifndef __CLI_PARAMS_H_
#define __CLI_PARAMS_H_

/*
 * Very simple command-line parameter parser.
 */

#include <string>
#include <vector>
#include <map>
#include <sstream>

class cli_params
{
public:
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

	std::string const & get_program_name() const
	{
		return program_name;
	}

	std::vector<std::string> const & get_positional() const
	{
		return positional;
	}

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
