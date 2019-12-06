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
	vector<char>* prot_active; //Contiendra la protein lu en cours
	map<int,char> prot_dictionnary ; //Dictionnaire pour lire le fichier binaire cree apr BLAST
	
	//Attribut header
	char* database_prot_header; //Contiendra la database phr
	vector<char>* prot_header_active;
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
	char* read_file(const string filepath);//Function utilitaire
	void generate_prot_index(const string filepath);
	unsigned int number_of_character(unsigned int* position_in_header,string first_byte);
	
public :
	Handle_Database(const string database_path) ;
	~Handle_Database();
	void update_next_protein_header(); //Permet de mettre a jour prot_active_header
	void fetch_prot_sequence(const unsigned int index);//Met dans prot_active la sequence de la prot designe par l index
	void fetch_prot_header(const unsigned int index);
	const vector<char>* get_prot_active();
	const vector<char>* get_prot_header_active();
	const unsigned int get_database_size();
	const unsigned int get_size_sequence_prot(const unsigned int index);

};
