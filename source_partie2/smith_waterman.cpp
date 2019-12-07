#include "smith_waterman.h"

//Constante :
const size_t NUMBER_OF_MAX_SAVED = 10 ;

Smith_Waterman::Smith_Waterman(const string filepath,string* query_protein)
{
	this->build_blossum_matrix(filepath);
	this->query_protein = query_protein ;
	this->gap_extension = 1;
	this->gap_opener = 11;
}

Smith_Waterman::~Smith_Waterman()
{
	delete this->blossum_matrix;
}


/*void Smith_Waterman::build_blossum_matrix(const string filepath) // Version Vector[]
{
	
	map<char,int> blos_dictionnary;
	blos_dictionnary= { {'A',0}, {'R',1}, {'N',2}, {'D',3}, {'C',4}, {'Q',5}, {'E',6}, {'G',7}, {'H',8}, {'I',9}, {'L',10}, {'K',11}, {'M',12}, {'F',13}, {'P',14}, {'S',15}, {'T',16}, {'W',17}, {'Y',18}, {'V',19}, {'B',20}, {'Z',21}, {'X',22},{'*',23}};
	vector<vector<int>>* matri= new vector<vector<int>>();
	matri->resize(24);
	for (int i = 0; i<24; i++)
		matri->at(i).resize(24);

int matrice(string Res1, string Res2, const map<char,int>* blos_dictionnary, const vector<vector<int>>* matri){
	i =0;
	res=0;

	while(Res1[i] != ' ' && Res2[i]!= ' '){
		res+=matri->at((*blos_dictionnary)[Res1[i]])[(*blos_dictionnary)[Res2[i]]);
		i++;
	}
	cout <<"Nombre de caractères: " << i << " Résultat: " << res << endl;
	return res;
}
	
	
	
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
}*/

void Smith_Waterman::build_blossum_matrix(const string filepath) //Version Map[]
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
	 *					 - Normaliser le score brut obtenu et le sauvegarder
	 * 					   Sbit = (λ S - ln K)/ ln 2 avec λ = 0.267 et ln(k) = -3.34
	 * */
	 
	 cout << endl << "Debut score_protein : " << endl ;
	 
	 //Cree une matrice de score avec colonnes database et ligne protein query
	 unsigned int size_prot_database;
	 const size_t size_prot_query = this->query_protein->size();
	 vector<vector<int>> score_matrix ; // Contiendra toute les valeurs calcule
	 vector<int> null_vector ;
	 vector<int> line_constructed ;
	 
	 unsigned int index_max_column[size_prot_query+1];//Contiendra l index de la val a utilise pour els gap top
	 unsigned int index_max_line =0 ; // Contiendra la val max de la ligne a utiliser pour les gap left
	 for(unsigned int i=0;i<=size_prot_query; ++i)
	 {
		 null_vector.push_back(0);
		 index_max_column[i]=0;
	 }
	 
	 unsigned int max_abs = 0; //Contiendra le max abs
	 unsigned int index_max_score_line=0; //Contiendra l index de la ligne ou est situe le max relatif au gap
	 unsigned int index_max_score_column=0 ; //Contiendra l index de la colonne ou est situe le max relatif au gap
	 unsigned int max_saved[NUMBER_OF_MAX_SAVED]; //Contiendra dans l ordre decroissant les meilleurs scores normalise
	 fill(max_saved, max_saved+NUMBER_OF_MAX_SAVED,0); //remplis de 0
	 unsigned int index_max_saved[NUMBER_OF_MAX_SAVED]; //Contiendra l index des proteines avec un bon score normalise
	 fill(index_max_saved, index_max_saved+NUMBER_OF_MAX_SAVED,0);
	 int score_left_gap = 0;
	 int score_up_gap = 0 ;
	 
	 map<char,int> *map_database_prot_tested = NULL; // Permet de retenir la map corrpondant au residu test
	 char residu_query;
	 char residu_database ;
	 unsigned int score_saved;
	 for(unsigned int index=115000; index<120000; ++index) //Essais sur les 115.000 128.000 soit 17.000
	 {
		 //Initialisation des variables non constante entre chaque test de protein
		size_prot_database = database->get_size_sequence_prot(index);
		database->fetch_prot_sequence(index); // met a jour la prot _active
		score_matrix.push_back(null_vector); //Premiere ligne de 0
		max_abs = 0 ;
		
		for(unsigned int i=1; i<=size_prot_database; ++i)
		{
			residu_database = database->get_prot_active()->at(i-1);
			map_database_prot_tested = &(this->blossum_matrix->at(residu_database)) ;
			line_constructed.push_back(0); //premiere colonne de la matrix =  zero
			index_max_line = 0 ;
			
			for(unsigned int j=1; j<=size_prot_query; ++j)
			{
				residu_query = this->query_protein->at(j-1);
				
				score_up_gap = score_matrix[index_max_column[j]][j] - this->gap_opener - (i-index_max_column[j])*this->gap_extension;
				score_left_gap = line_constructed[index_max_line] - this->gap_opener - (j-index_max_line)*this->gap_extension;
				
				score_saved = this->max_over_zero(score_left_gap, score_up_gap, 
							  map_database_prot_tested->at(residu_query)+score_matrix[i-1][j-1]) ;
				line_constructed.push_back( score_saved);
				
				if(max_abs <= score_saved) // mise a jour du maximum
				{
					max_abs = score_saved ;
					index_max_score_line = i;
					index_max_score_column = j;
				}
											
				//Check et changement d index pour calculer les gap correctement
				if(score_up_gap < (line_constructed[j]-this->gap_opener)){index_max_column[j]=i;}
				if(score_left_gap < (line_constructed[j]-this->gap_opener)){ index_max_line=j;}
			}
			score_matrix.push_back(line_constructed);
			line_constructed.clear() ;
			map_database_prot_tested = NULL;
		}
		
		locate_replace_max( index, floor((0.267*max_abs +3.34)/(log(2))), max_saved, index_max_saved); //Sbit = (λ S - ln K)/ ln 2 avec λ = 0.267 et ln(k) = -3.34
		
		for(unsigned int i=0;i<=size_prot_query; ++i)
		{
			index_max_column[i] = 0;
		}
		score_matrix.clear() ; //Important pour passer de prot en prot
	 }
	 
	 for(size_t i=0; i<NUMBER_OF_MAX_SAVED; ++i)
	 {
		 cout << "Score Obtenue : " << max_saved[i] << " pour la prot : " << index_max_saved[i] << endl ;
	 }
	 
	return 0;
}

int Smith_Waterman::max_over_zero(int up, int left, int diag) const
{
	int return_value = left ;
	
	if(up>return_value){return_value=up;}
	if(diag>return_value){return_value=diag;}
	if(return_value<0){return_value=0;}
	
	return return_value;
}

void Smith_Waterman::locate_replace_max(const unsigned int index,const unsigned int value, unsigned int max_table[], unsigned int index_max_table[]) 
{
	if(value > max_table[NUMBER_OF_MAX_SAVED-1])//Si la valeur a test est superieur a la plus petite du tableau
	{
		//Possibilite de recherche dichotomique plus mais pas forcement interessante si vraiment petit tableau
		max_table[NUMBER_OF_MAX_SAVED-1] = value ;
		index_max_table[NUMBER_OF_MAX_SAVED-1] = index ;
		int pos_found = NUMBER_OF_MAX_SAVED-2 ;
		while( (value > max_table[pos_found]) and (pos_found>=0) )
		{
			max_table[pos_found+1] = max_table[pos_found] ;
			index_max_table[pos_found+1] = index_max_table[pos_found] ;
			max_table[pos_found] = value ;
			index_max_table[pos_found] = index ;
			--pos_found ; //Si pas trouve on monte dans le tableau
			
		}
	}
}

void Smith_Waterman::set_gap_opener(const unsigned int new_gap_opener){this->gap_opener=new_gap_opener;}
void Smith_Waterman::set_gap_extension(const unsigned int new_gap_extension){this->gap_extension = new_gap_extension;}
