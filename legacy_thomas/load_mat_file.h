#ifndef __LOAD_MAT_FILE_H_
#define __LOAD_MAT_FILE_H_

#include <vector>
#include <fstream>
#include <string>
#include "material.hh"
#include "archive.hh"

material_legacy_thomas load_mat_file(std::string const & filename)
{
	material_legacy_thomas mat;
	
	std::ifstream ifs(filename, std::ifstream::binary);
	if(!ifs.is_open())
	{
		throw std::ios_base::failure("failed to open '"+filename+"' for reading");
	}

	archive::istream ia(ifs);
	ia >> mat;
	ifs.close();
	return mat;
}

#endif // __LOAD_MAT_FILE_H_
