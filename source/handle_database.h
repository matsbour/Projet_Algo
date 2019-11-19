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
	//Attribut equence
	char* prot_sequence ; //Contiendra toute la database psq
	int index_prot_sequence; //Permet de se souvenir a quelle proteine nous sommes
	vector<char>* prot_active; //Contiendra la protein lu en cours
	map<int,char> prot_dictionnary ; //Dictionnaire pour lire le fichier binaire cree apr BLAST
	
	//Attribut header
	char* prot_header; //Contiendra la database phr
	int index_prot_header;
	vector<char>* prot_header_active;
	unsigned int number_of_character(unsigned int* position_in_header,string first_byte);
	
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
	u_int32_t*  header_offset;
	u_int32_t*  sequence_offset ;
	
	//Fonction
	char* read_file(const string filepath);//Function utilitaire
	void generate_prot_index(const string filepath);
	
public :
	Handle_Database(const string database_path) ;
	~Handle_Database();
	void update_next_protein_sequence(); //Permet de mettre a jour prot_active
	void update_next_protein_header(); //Permet de mettre a jour prot_active
	const vector<char>* get_prot_active();
	const vector<char>* get_prot_header_active();

};

