#include "smith_waterman.h"

Smith_Waterman::Smith_Waterman(const string filepath,string* query_protein)
{
	this->build_blossum_matrix(filepath);
	this->query_protein = query_protein ;
}

Smith_Waterman::~Smith_Waterman()
{
	delete this->blossum_matrix;
}

void Smith_Waterman::build_blossum_matrix(const string filepath)
{
	this->blossum_matrix = new map<char,map<char,int>>;
	string container, order_of_residu;
	ifstream file(filepath, std::ifstream::binary);
	if(file.is_open())
	{
		while(getline(file,container))// permet d eviter les premieres ligne
		{
			if(container[0]!= '#')
			{
				order_of_residu = container;
				break;
			} 
		}
		
		size_t pos_space = 0 ;
		while(true)//Efface les espaces order_of_residu
		{
			pos_space = order_of_residu.find(" ");
			if(pos_space == string::npos){break;}//Il n y plus d espace dans la chaine
			order_of_residu.erase(pos_space,1); //Efface un espace
			
		}
		
		char start_line = ' ';
		map<char,int> line_map ;
		size_t compteur_residu = 0; // permet de savoir auxquelles nous sommes
		while(getline(file,container) and !file.eof())//ligne par ligne
		{
			compteur_residu=0;
			start_line = container[0];
			for(unsigned int i=1; i<container.size();++i)
			{
				if(container[i]!=' ')
				{
					if(container[i]=='-') // si nb negatif
					{
						//Attention convertion char to int mais valeur pas code ascii : (int)char - (int) '0'	
						line_map.insert(pair<char,int>(order_of_residu[compteur_residu],-1*((int)container[i+1] - (int)'0')));		
						++i;
						++compteur_residu; // Passe au residu suivant 
					}
					else//nb a regarder positif
					{
						line_map.insert(pair<char,int>(order_of_residu[compteur_residu],(int)container[i] - (int)'0'));
						++compteur_residu; // Passe au residu suivant 
					}
				}
			}
			this->blossum_matrix->insert(pair<char,map<char,int>>(start_line,line_map));// Ajoute une map a un character dans la matric
			line_map.clear();
		} 
	}
	else{
	cout<<"Cannot open file "<< filepath <<endl;
	exit(1);}
	
	
	file.close();
}

unsigned int Smith_Waterman::score_protein(Handle_Database* database)
{
	/* Plusieurs etapes :- Etablir la matrice de score et retenir le max
	 *					 - Normaliser le score brut obtenu et le retourner */
	 
	 cout << "Database size :" << database->get_database_size() << endl;

	return 0;
}

