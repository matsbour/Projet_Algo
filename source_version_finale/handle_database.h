#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <map>
#include <cstring>
#include <vector>
#include <iomanip>
#include <math.h>
#include <endian.h> // for bswap
using namespace std;

class Handle_Database
{
private :
	//Attribut
	string database_path_saved ;
	//Attribut sequence
	char* database_prot_sequence ; //Contiendra toute la database psq
	
	//Attribut header
	char* database_prot_header; //Contiendra la database phr
	map<char,int> hex2int_map ;
	
	//Attribut index
	u_int32_t version;
	u_int32_t database_type ;
	u_int32_t  title_length ;
	char* title ;
	u_int32_t  timestamp_length ;
	char* timestamp;
	u_int32_t  numbers_of_sequence ;
	uint64_t numbers_of_residues;
	u_int32_t  prot_max_length ;
	vector<int>*  header_offset_vector;
	vector<int>*  sequence_offset_vector ;
	
	//Fonction
	char* read_file(const string filepath);
	void generate_prot_index(const string filepath);
	unsigned int number_of_character(unsigned int* position_in_header,string first_byte);
	
public :
	Handle_Database(const string database_path) ;
	~Handle_Database();
	void update_next_protein_header(); //Permet de mettre a jour prot_active_header
	char* fetch_prot_sequence_residu(const unsigned int index, const unsigned int offset);//Renvoie le residu demande par index/offset
	string fetch_prot_header(const unsigned int index);
	
	//Get et set
	u_int32_t get_version();
	u_int32_t get_database_type();
	u_int32_t  get_title_length();
	char* get_title();
	u_int32_t  get_timestamp_length();
	char* get_timestamp();
	u_int32_t  get_numbers_of_sequence();
	uint64_t get_numbers_of_residues();
	u_int32_t  get_prot_max_length();
	const unsigned int get_database_size();
	const unsigned int get_size_sequence_prot(const unsigned int index);
	

};

