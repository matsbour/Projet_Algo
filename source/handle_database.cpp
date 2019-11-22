#include "handle_database.h"

Handle_Database::Handle_Database(const string database_path)
{
	this->database_path_saved = database_path;
	this->index_prot_sequence = 0; //1 Car le premier character est '-'
	this->index_prot_header = 0; 
	this->prot_active = new vector<char>; //Cree sur le tas
	this->prot_header_active = new vector<char>;
	this->prot_dictionnary ={
		{0,'-'}, {1,'A'}, {2,'B'},{3,'C'},{4,'D'},
		{5,'E'}, {6,'F'}, {7,'G'},{8,'H'},{9,'I'},
		{27,'J'}, {10,'K'}, {11,'L'},{12,'M'}, {13,'N'},
		{26,'O'}, {14,'P'}, {15,'Q'},{16,'R'}, {17,'S'},
		{18,'T'}, {24,'U'}, {19,'V'},{20,'W'}, {21,'X'},
		{22,'Y'}, {23,'Z'}, {25,'*'}
	};
	this->hex2int_map = { 
		{'0',0}, {'1',1}, {'2',2}, {'3',3}, {'4',4},
		{'5',5}, {'6',6}, {'7',7}, {'8',8}, {'9',9},
		{'A',10}, {'B',11}, {'C',12}, {'D',14}, {'E',15}, {'F',16},
		{'a',10}, {'b',11}, {'c',12}, {'d',14}, {'e',15}, {'f',16}
	};
	this->header_offset_vector = new vector<int>;
	this->sequence_offset_vector = new vector<int>;
	this->generate_prot_index(database_path+".pin");
	this->sequence_stream.open(database_path+".psq", std::ifstream::binary);
	this->header_stream.open(database_path+".phr", std::ifstream::binary);
}

Handle_Database::~Handle_Database()
{
	delete this->prot_active;
	delete this->prot_sequence;
	delete this->prot_header;
	delete this->prot_header_active;
	delete this->title;
	delete this->timestamp;
	delete this->header_offset_vector;
	delete this->sequence_offset_vector;
	this->sequence_stream.close();
	this->header_stream.close();
}

char* Handle_Database::read_file(const string filepath)
{
	ifstream file(filepath, std::ifstream::binary);
	char* char_container = NULL; // Servira a cree la zone memoire a retourne
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

const vector<char>* Handle_Database::get_prot_active(){return this->prot_active;}
const vector<char>* Handle_Database::get_prot_header_active(){return this->prot_header_active;}
const unsigned int Handle_Database::get_database_size(){return this->sequence_offset_vector->size();}
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
		cout<<"Version : " << version<<endl;
	
		file.read((char*)&database_type,sizeof(uint32_t));
		database_type = __builtin_bswap32(database_type);
		cout<<"Database Type : " << (int)database_type<<endl;
	
		file.read((char*)&title_length,sizeof(uint32_t));
		title_length = __builtin_bswap32(title_length);
		cout<<"Length title: " << (int)title_length<<endl;
		
		title = new char[(int)title_length+1];
		file.read(title,(int)title_length);
		title[(int)title_length] = '\0';
		cout << "Title : " ;
		for(unsigned int i=0;i<strlen(title);++i){cout<<title[i];}
		cout << endl ;
		
		file.read((char*)&timestamp_length,sizeof(uint32_t));
		timestamp_length = __builtin_bswap32(timestamp_length);
		cout<<"Timestamp length: " << (int)timestamp_length<<endl;
		
		timestamp = new char[(int)timestamp_length+1];
		file.read(timestamp,(int)timestamp_length);
		timestamp[(int)timestamp_length] = '\0';
		cout << "Database creation  : " ;
		for(unsigned int i=0;i<strlen(timestamp);++i){cout<<timestamp[i];}
		cout << endl ;
		
		file.read((char*)&numbers_of_sequence,sizeof(uint32_t));
		numbers_of_sequence = __builtin_bswap32(numbers_of_sequence);
		cout<<"Numbers of sequence: " << (int)numbers_of_sequence<<endl;
		
		file.read((char*)&numbers_of_residues,sizeof(uint64_t));
		cout<<"Numbers of residues: " << (int)numbers_of_residues<<endl;
		
		file.read((char*)&prot_max_length,sizeof(uint32_t));
		prot_max_length = __builtin_bswap32(prot_max_length);
		cout<<"Length max: " << (int)prot_max_length<<endl;
		 
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

void Handle_Database::fetch_prot_sequence(const unsigned int index)
{ 
	if(index >= this->sequence_offset_vector->size()) // on verifie si on a un numero de prot trop grand
	{
		cout<<"Index is out of bound for sequence_offset" << endl;
		exit(1);
	}
	
	unsigned int length = (int)this->sequence_offset_vector->at(index+1) - (int)this->sequence_offset_vector->at(index);
	--length; //le dernier byte est un byte null
	string filepath = (this->database_path_saved+".psq");
	char* char_container = NULL; // Servira a cree la zone memoire a retourne
	if(sequence_stream.is_open())
	{
		sequence_stream.seekg((int)this->sequence_offset_vector->at(index),sequence_stream.beg);
		char_container = new char [length] ; //servira a stocke tous les char du fichier
		sequence_stream.read(char_container,length);
	}
	else{
	cout<<"Cannot open file "<< this->database_path_saved << ".psq" <<endl;
	exit(1);}
	
	this->prot_active->clear() ;
	for(unsigned int i=0; i<length;++i)
	{
		prot_active->push_back(prot_dictionnary[(int)char_container[i]]);
	}
	this->index_prot_sequence = index;
	delete char_container;
	
}

void Handle_Database::fetch_prot_header(const unsigned int index)
{
	if(index >= this->header_offset_vector->size()) // on verifie si on a un numero de prot trop grand
	{
		cout<<"Index is out of bound for sequence_offset" << endl;
		exit(1);
	}
	
	unsigned int length = (int)this->header_offset_vector->at(index+1) - (int)this->header_offset_vector->at(index);
	cout << "Length :" << length << endl;
	string filepath = (this->database_path_saved+".phr");
	char* char_container = new char [length]; // servira a stocke tous les char du fichier
	if(header_stream.is_open())
	{
		header_stream.seekg((int)this->header_offset_vector->at(index),header_stream.beg);
		header_stream.read(char_container,length/sizeof(char));  //?
	}
	else{
	cout<<"Cannot open file "<< filepath <<endl;
	exit(1);}
	
	this->prot_header_active->clear() ;
	unsigned int size_of_visible_string = 0 ;
	unsigned int int_test = 0;
	unsigned int position_start_string = 0 ;
	std::stringstream stream;
	for(unsigned int i=0; i< length;++i )
	{
		int_test = (unsigned int) (u_int8_t)char_container[i];
		if(int_test == 26) // Car 1A = 26, on ignore les premier byte designant les types
		{
			int_test = (unsigned int) (u_int8_t)char_container[i+1];
			stream << std::hex << int_test;
			size_of_visible_string = this->number_of_character(&i,stream.str(),char_container) ;
			position_start_string = i+2;// si on dit que le premier hex ne peut pas depasser 8 
			stream.str(string()); //vide le stream
			break;
		}
	}
	
	for(unsigned int j=position_start_string; j<=position_start_string+size_of_visible_string; ++j)
	{
		prot_header_active->push_back((char) char_container[j]);
	}
	cout << endl;
	delete char_container ;
	this->index_prot_header=index;
}

//Fonction permettant de determiner la longueur de la chaine de character visible
unsigned int Handle_Database::number_of_character(unsigned int* position_in_header,string first_byte, char* prot_header)
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
			character_test = (unsigned int) (u_int8_t)prot_header[i+(*position_in_header)];
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

