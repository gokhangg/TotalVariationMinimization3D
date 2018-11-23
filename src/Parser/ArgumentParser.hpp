/*
 * Project: Commandline parser
 * Author: Gokhan Gunay, ghngunay@gmail.com
 * Copyright: (C) 2018 by Gokahn Gunay
 * License: GNU GPL v3 (see License.txt)
 */

#ifndef parser_hpp
#define parser_hpp

#include <cstring>
#include <string>


using namespace std;

string get_cstr_arg(string str, const char * arg);

class cparser {
private:
	string input_args;
	string keys;
	unsigned int saved_key_num;
public:

	cparser(int argc,char ** argv);
	cparser();
	~cparser();
	void input(int argc, char ** arg);
	void save_keys(const char* in_arg, const char* key);
	char* return_val(char* arg);
	string operator[](const char* key);
	unsigned int get_saved_key_num(void);
};

/*
* @brief Constructor.
*/
cparser::cparser(int argc, char ** argv)
{
	input(argc, argv);
	saved_key_num = 0;
}

/*
* @brief Default constructor.
*/
cparser::cparser()
{

}

/*
* @brief Destructor.
*/
cparser::~cparser()
{
	
}

/*
* @brief Importing input arguments and count to the argument pool.
*/
void cparser::input(int argc,char ** arg)
{
	for (int cnt=0;cnt<argc;cnt++)
	{
		input_args.append(" ");
		input_args.append(arg[cnt]);
	}
}

/*
* @brief Save a key and corresponding possible input argument.
*/
void cparser::save_keys(const char* key, const  char* in_arg)
{
	if ((key == NULL) || (in_arg == NULL))//If one of the arguments are invalid return without any action.
		return;
	keys.append(key);//append the key after the input argument
	keys.append(" ");//add one space
	keys.append(in_arg);//save input argument
	keys.append(" ");//add one space
	saved_key_num++;//Save the number of saved keys
}

/*
* @brief Retrieving input argument corresponding to the given key.
*/
string cparser::operator[](const char* key)
{
	return get_cstr_arg(input_args, get_cstr_arg(keys, key).c_str());
}

/*
* @brief Getting number of the saved keys.
*/
unsigned int cparser::get_saved_key_num(void)
{
	return saved_key_num;
}

/**
* @brief Retrieving input argument corresponding to the given key.
*/
string get_cstr_arg(string str, const char * arg)
{
	size_t pos=str.find(arg);
	if (pos==-1)
		return static_cast<string>("");//if not found return null
	pos=str.find(" ", pos+1);
	return str.substr(pos+1,str.find(" ", pos+1)-pos-1);
}


#endif
