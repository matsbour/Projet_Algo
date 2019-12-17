#include "handle_database.h"

Handle_Database::Handle_Database(const string database_path)
{
	this->database_path_saved = database_path; //sauvegarde le chemin d'accès
	this->database_prot_sequence = NULL;
	this->database_prot_header = NULL;
	this->hex2int_map = { //valeurs hexadécimales et décimales correspondantes
		{'0',0}, {'1',1}, {'2',2}, {'3',3}, {'4',4},
		{'5',5}, {'6',6}, {'7',7}, {'8',8}, {'9',9},
		{'A',10}, {'B',11}, {'C',12}, {'D',14}, {'E',15}, {'F',16},
		{'a',10}, {'b',11}, {'c',12}, {'d',14}, {'e',15}, {'f',16}
	};
	this->header_offset_vector = new vector<int>;
	this->sequence_offset_vector = new vector<int>;
	this->generate_prot_index(database_path+".pin"); //index file
	this->database_prot_sequence = this->read_file(database_path+".psq"); //proteine sequence file
	this->database_prot_header = this->read_file(database_path+".phr"); //header file
}

Handle_Database::~Handle_Database()
{
	delete this->database_prot_sequence;
	delete this->database_prot_header;
	delete this->title;
	delete this->timestamp;
	delete this->header_offset_vector;
	delete this->sequence_offset_vector;
}

char* Handle_Database::read_file(const string filepath)
{
	ifstream file(filepath, ios::in | ios::binary);
	char* char_container = NULL;
	if(file.is_open())
	{
		file.seekg(0,file.end);
		int length = file.tellg(); //Nombre de characters dans le fichier a lire
		file.seekg(0,file.beg);
		char_container = new char [length] ; //servira a stocke tous les char du fichier
		file.read(char_container,length);
		file.close();
	}
	else{
	cout<<"Cannot open file "<< filepath <<endl;
	exit(1);}
	
	return char_container;
}

const unsigned int Handle_Database::get_database_size(){return (this->sequence_offset_vector->size()-1);}
const unsigned int Handle_Database::get_size_sequence_prot(const unsigned int index)
{
	if(index <(this->sequence_offset_vector->size())-1)
	{
		return (unsigned int)(this->sequence_offset_vector->at(index+1) - this->sequence_offset_vector->at(index))-1 ;
	}
	else
	{
		cout << "Out of bound : " << index << " size limit : " << (this->sequence_offset_vector->size()-1) << endl ;
		exit(1);
	}
	return 0;
	}

//Remplis tous les attributs necessaie depuis le fichier .pin
void Handle_Database::generate_prot_index(string filepath)
{
	ifstream file(filepath,std::ifstream::binary);
	if(file.is_open())
	{
		file.read((char*)&version,sizeof(uint32_t));
		version = __builtin_bswap32(version);
	
		file.read((char*)&database_type,sizeof(uint32_t));
		database_type = __builtin_bswap32(database_type);
	
		file.read((char*)&title_length,sizeof(uint32_t));
		title_length = __builtin_bswap32(title_length);
		
		title = new char[(int)title_length+1];
		file.read(title,(int)title_length);
		title[(int)title_length] = '\0';
		
		file.read((char*)&timestamp_length,sizeof(uint32_t));
		timestamp_length = __builtin_bswap32(timestamp_length);
		
		timestamp = new char[(int)timestamp_length+1];
		file.read(timestamp,(int)timestamp_length);
		timestamp[(int)timestamp_length] = '\0';
		
		file.read((char*)&numbers_of_sequence,sizeof(uint32_t));
		numbers_of_sequence = __builtin_bswap32(numbers_of_sequence);
		
		file.read((char*)&numbers_of_residues,sizeof(uint64_t));
		
		file.read((char*)&prot_max_length,sizeof(uint32_t));
		prot_max_length = __builtin_bswap32(prot_max_length);
		 
		u_int32_t* header_offset = new u_int32_t[(int)numbers_of_sequence+1] ; // Tableau d offset representant l ecart entre chaque header
		for(int i=0; i<(int)numbers_of_sequence+1;++i)
		{
			file.read((char*)&(header_offset[i]),sizeof(uint32_t));
			header_offset_vector->push_back((int)__builtin_bswap32(header_offset[i]));
		}
	
		u_int32_t* sequence_offset = new u_int32_t[(int)numbers_of_sequence+1] ; // Tableau d offset representant l ecart entre chaque sequence
		for(int i=0; i<(int)numbers_of_sequence+1;++i)
		{
			file.read((char*)&(sequence_offset[i]),sizeof(uint32_t));
			sequence_offset_vector->push_back((int)__builtin_bswap32(sequence_offset[i]));
		}
		
		delete header_offset;
		delete sequence_offset;
	}
	else{cout << "Cannot read: " << filepath<<endl;}
	file.close();
}

char* Handle_Database::fetch_prot_sequence_residu(const unsigned int index, const unsigned int offset)
{ 
	if(index > this->sequence_offset_vector->size()) // on verifie si on a un numero de prot trop grand
	{
		cout<<"Index is out of bound for sequence_offset" << endl;
		exit(1);
	}
	
	return &( this->database_prot_sequence[(int)(this->sequence_offset_vector->at(index))+offset]);
}

string Handle_Database::fetch_prot_header(const unsigned int index)
{
	if(index >= this->header_offset_vector->size()) // on verifie si on a un numero de prot trop grand
	{
		cout<<"Index is out of bound for header_offset" << endl;
		exit(1);
	}
	
	string value_return = " ";
	unsigned int length = (int)this->header_offset_vector->at(index+1) - (int)this->header_offset_vector->at(index);
	unsigned int size_of_visible_string = 0 ;
	unsigned int int_test = 0;
	unsigned int position_start_string = 0 ;
	unsigned int position_start = (unsigned int)this->header_offset_vector->at(index) ;
	std::stringstream stream;
	for(unsigned int i=position_start; i< length+position_start;++i )
	{
		int_test = (unsigned int) (u_int8_t)database_prot_header[i];
		if(int_test == 26) // Car 1A = 26, on ignore les premier byte designant les types
		{
			int_test = (unsigned int) (u_int8_t)database_prot_header[i+1];
			stream << std::hex << int_test;
			size_of_visible_string = this->number_of_character(&i,stream.str()) ;
			position_start_string = i+2;// si on dit que le premier hex ne peut pas depasser 8 
			stream.str(string()); //vide le stream
			break;
		}
	}

	for(unsigned int j=position_start_string; j<position_start_string+size_of_visible_string; ++j)
	{
		value_return += ((char) database_prot_header[j]);
	}
	
	return value_return;
}

//Fonction permettant de determiner la longueur de la chaine de character visible
unsigned int Handle_Database::number_of_character(unsigned int* position_in_header,string first_byte)
{
	//Dans le cas ou le premier hex est > 8 on doit lire la taille sur les bytes suivant
	if( this->hex2int_map[first_byte[0]] >= 8 )
	{	
		unsigned int number_of_byte = this->hex2int_map[first_byte[1]] ;
		unsigned int size_saved_int = 0;
		unsigned int character_test;
		string size_saved_string ;
		std::stringstream stream;
		for(unsigned int i=2; i<number_of_byte+2;++i)
		{
			character_test = (unsigned int) (u_int8_t)database_prot_header[i+(*position_in_header)];
			stream << std::hex << character_test;
		}
		
		size_saved_string = stream.str();
		unsigned int temp_size = size_saved_string.size();
		for(unsigned int j=0;j<temp_size; ++j)
		{
			size_saved_int += (hex2int_map[size_saved_string[j]])*pow(16,(temp_size-1)-j); // calcul hexadecimal
		}
		*position_in_header += number_of_byte; // met a jour la nouvelle position limite du debut de la string
		
		return size_saved_int;
	}
	//Cas ou les hex donne directement la taille de la string
	else
	{
		return this->hex2int_map[first_byte[0]]*16 + this->hex2int_map[first_byte[1]] ;
	}
	return 0 ;
}



u_int32_t Handle_Database::get_version(){return version;}
u_int32_t Handle_Database::get_database_type(){return database_type;}
u_int32_t  Handle_Database::get_title_length(){return title_length;}
char* Handle_Database::get_title(){return title;}
u_int32_t  Handle_Database::get_timestamp_length(){return timestamp_length;}
char* Handle_Database::get_timestamp(){return timestamp;}
u_int32_t  Handle_Database::get_numbers_of_sequence(){return numbers_of_sequence;}
uint64_t Handle_Database::get_numbers_of_residues(){return this->numbers_of_residues;}
u_int32_t  Handle_Database::get_prot_max_length(){return prot_max_length;}
