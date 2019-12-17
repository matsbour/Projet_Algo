#include "smith_waterman.h"

//Constante :
const size_t NUMBER_OF_MAX_SAVED = 10 ; //nombre de maximum sauvegardés, valeur constante 

Smith_Waterman::Smith_Waterman(const string filepath,myProtein* query_protein_ini, int gap_opener_penalty, int gap_extension_penalty )
{
	this->prot_dictionnary ={ //dictionnaire du format file sequence
		{'-',0}, {'A',1}, {'B',2},{'C',3},{'D',4},
		{'E',5}, {'F',6}, {'G',7},{'H',8},{'I',9},
		{'J',27}, {'K',10}, {'L',11},{'M',12}, {'N',13},
		{'O',26}, {'P',14}, {'Q',15},{'R',16}, {'S',17},
		{'T',18}, {'U',24}, {'V',19},{'W',20}, {'X',21},
		{'Y',22}, {'Z',23}, {'*',25}
	};
		
	this->build_blossum_matrix(filepath); //construit la matrice BLOSUM
	
	query_protein = new vector<int> ; //construit le vecteur de la query protein en cherchant la valeur correspondant à la lettre dans le dictionnaire
	for(size_t i=0; i<query_protein_ini->getSequence()->size();++i)
	{
		this->query_protein->push_back( prot_dictionnary[(query_protein_ini->getSequence()->at(i))]);
	}
	
	this->query_protein_header = query_protein_ini->getHeader();
	
	this->gap_opener = gap_opener_penalty; //pénalité dans le cas où on ouvre un gap
	this->gap_extension = gap_extension_penalty; //pénalité dans le cas où on propage un gap
	
}

Smith_Waterman::~Smith_Waterman()
{
	delete this->blossum_matrix;
	delete this->query_protein;
}


void Smith_Waterman::build_blossum_matrix(const string filepath) 
{
	/*
	* @desc 
	* @param 
	**//
	
	int flag = 5000; // Sert a remplir la matrice de base
	this->blossum_matrix = new vector<vector<int>>();
	blossum_matrix->resize(28); //28 si uniquement des matrices blosum et il y a 28 éléments dans le prot_dictionnary
	for (int i = 0; i<28; i++)
		blossum_matrix->at(i).resize(28);

	string container;
	vector<int> order_of_residu;
	ifstream file(filepath, std::ifstream::binary);
	if(file.is_open())
	{
		while(getline(file,container))
		{
			if(container[0]!= '#') //ne prend pas en compte les 6 premières lignes du fichier BLOSUM62
			{
				for(size_t i=0; i<container.size();++i)
				{
					if(container[i]!=' '){order_of_residu.push_back( this->prot_dictionnary[container[i]]);} //enlève les espaces
				}
				break;
			} 
		}
		
		
		int start_line = 0;
		vector<int> line_vect;
		for(size_t i=0; i<28; ++i){line_vect.push_back(flag);} //on remplit le vect line
		for(size_t i=0; i<28; ++i){blossum_matrix->at(i) = line_vect;} //on remplit le vect line
		
		
		size_t compteur_residu = 0; // permet de savoir auxquelles nous sommes
		while(getline(file,container) and !file.eof())//ligne par ligne
		{
			compteur_residu=0;
			start_line = prot_dictionnary[container[0]];
			for(unsigned int i=1; i<container.size();++i)
			{
				if(container[i]!=' ')
				{
					if(container[i]=='-') // si nb negatif
					{
						//Attention convertion char to int mais valeur pas code ascii : (int)char - (int) '0'	
						line_vect[order_of_residu[compteur_residu]] = ((-1)*((int)container[i+1] - (int)'0'));		
						++i;
						++compteur_residu; // Passe au residu suivant 
					}
					else//nb a regarder positif
					{
						if(container[i+1]!=' ') // A cause de 11 pour W deux characters
						{
							line_vect[order_of_residu[compteur_residu]] = (((int)container[i]-(int)'0')*10) + ((int)(container[i+1])- (int)'0' );
							++i;
						}
						else{line_vect[order_of_residu[compteur_residu]] = ((int)container[i] - (int)'0');}
						++compteur_residu; // Passe au residu suivant 
					}
				}
			}
			this->blossum_matrix->at(start_line) =  line_vect;
			for(size_t i=0; i<28; ++i){line_vect.push_back(flag);} //on remplit le vect line
		} 
	}
	else{
	cout<<"Cannot open file : ["<< filepath << "] in order to create the blossum matrix" <<endl;
	exit(1);}
	
	int default_value_different = 0 ;//Value qui sert a remplir les cas X:*
	for(size_t i=0;i<28;++i)
	{
		default_value_different = blossum_matrix->at(prot_dictionnary['*']).at(i); // valeur de X:*
		if(default_value_different != flag){break;}
	}
	int default_value_same = blossum_matrix->at(prot_dictionnary['*']).at(prot_dictionnary['*']); //valeur de *:*
	for(size_t i=0; i< 28; ++i)
	{
		for(size_t j=0;j<28;++j)
		{
			if(blossum_matrix->at(i).at(j)==flag)
			{
				if(i==j){blossum_matrix->at(i).at(j)=default_value_same;}
				else{blossum_matrix->at(i).at(j)=default_value_different;}
			}
		}
	}
	
	file.close();
}

unsigned int Smith_Waterman::score_protein(Handle_Database* database)
{
	/*
	* @desc Calcule le score de la comparaison de 2 protéines
	* @param Handle_Database* : pointeur vers la database des protéines
	*
	*
	 *Plusieurs etapes :- Etablir la matrice de score et retenir le max
	 *		    - Normaliser le score brut obtenu et le sauvegarder
	 *                    Sbit = (λ S - ln K)/ ln 2 avec λ = 0.267 et ln(k) = -3.34
	 * */
	 
	 this->display_information(database) ;
	 if(database->get_database_size() == 0)
	 {
		 cout << "Database vide" << endl ;
		 exit(1);
	 }
	 
	 //Crée une matrice de score avec colonnes = protéine de la database et lignes = protein query qu'on veut comparer 
	 unsigned int size_prot_database;
	 const size_t size_prot_query = this->query_protein->size();
	 vector<int> vect_l1 ;
	 vector<int> vect_l2 ;
	 vector<int>* vect_saved; //Pointeur vers le vect qui stocke les valeurs de la ligne du dessus a celle calcule
	 vector<int> null_vector ;
	 vector<int> line_constructed ;
	 unsigned int max_saved[NUMBER_OF_MAX_SAVED]; //Contiendra dans l ordre decroissant les meilleurs scores normalise
	 unsigned int index_max_saved[NUMBER_OF_MAX_SAVED]; //Contiendra l index des proteines avec un bon score normalise
	 fill(index_max_saved, index_max_saved+NUMBER_OF_MAX_SAVED,0);
	 fill(max_saved, max_saved+NUMBER_OF_MAX_SAVED,0); //remplis de 0
	 unsigned int index_max_column[size_prot_query+1];//Contiendra l index de la val a utilise pour les gap top
	 unsigned int max_score_column[size_prot_query+1];
	 unsigned int index_max_line =0 ; // Contiendra la val max de la ligne a utiliser pour les gap left
	 for(unsigned int i=0;i<=size_prot_query; ++i)
	 {
		 null_vector.push_back(0);
		 index_max_column[i]=0;
		 max_score_column[i]=0;
		 vect_l1.push_back(0);
		 vect_l2.push_back(0);
	 }
	 
	 unsigned int max_abs = 0; //Contiendra le max abs
	 unsigned int index_max_score_line=0; //Contiendra l index de la ligne ou est situe le max relatif au gap
	 unsigned int index_max_score_column=0 ; //Contiendra l index de la colonne ou est situe le max relatif au gap
	 int score_left_gap = 0;
	 int score_up_gap = 0 ;
	 
	 vector<int> *vect_database_prot_tested = NULL; // Permet de retenir le vect correspondant au residu test
	 int* residu_query;
	 char* residu_database ;
	 int score_saved;
	 for(unsigned int index=0; index<7000; ++index) //Essais sur i prot de la database (size : database->get_database_size())
	 {
		 
		//Initialisation des variables non constante entre chaque test de protein
		size_prot_database = database->get_size_sequence_prot(index);
		max_abs = 0 ;
		
		for(unsigned int i=1; i<=size_prot_database; ++i)
		{
			residu_database = database->fetch_prot_sequence_residu(index,i-1);
			if(i%2==0){vect_saved = &vect_l1 ;}
			else{vect_saved = &vect_l2 ;}
			vect_database_prot_tested = &(this->blossum_matrix->at((int)(*residu_database))); 
			line_constructed.push_back(0); //premiere colonne de la matrix =  zero
			index_max_line = 0 ;
			
			for(unsigned int j=1; j<=size_prot_query; ++j)
			{
				residu_query = &(this->query_protein->at(j-1));
				
				score_up_gap = max_score_column[j] - this->gap_opener - (i-index_max_column[j])*(this->gap_extension);
				score_left_gap = line_constructed[index_max_line] - this->gap_opener - (j-index_max_line)*(this->gap_extension);
				

				score_saved = this->max_over_zero(score_left_gap, score_up_gap, 
						 vect_database_prot_tested->at(*residu_query)+vect_saved->at(j-1)) ; 
				
				line_constructed.push_back(score_saved);
				
				if(max_abs < score_saved) // mise a jour du maximum
				{
					max_abs = score_saved ;
					index_max_score_line = i;
					index_max_score_column = j;
				}
											
				//Check et changement d index pour calculer les gap correctement
				if(score_up_gap < (score_saved-this->gap_opener))
				{	index_max_column[j] = i ;
					max_score_column[j]=score_saved;}
				if(score_left_gap < (score_saved-this->gap_opener)){ index_max_line=j;}
			}
			
			if(i%2==0){vect_l2 = line_constructed;}
			else{vect_l1 = line_constructed ;}
			line_constructed.clear() ;
			vect_database_prot_tested = NULL;
		}
		
		locate_replace_max( index, floor((0.267*max_abs +3.34)/(log(2))), max_saved, index_max_saved); //Sbit = (λ S - ln K)/ ln 2 avec λ = 0.267 et ln(k) = -3.34
		
		for(unsigned int i=0;i<=size_prot_query; ++i)
		{
			index_max_column[i] = 0;
			max_score_column[i] = 0;
		}
	    //Important pour passer de prot en prot
		vect_l1 = null_vector ;
		vect_l2 = null_vector ;
	 }
	 
	 
	 this->display_max(max_saved, index_max_saved , database);

	 
	return 0;
}

int Smith_Waterman::max_over_zero(int up, int left, int diag) const
{
	/**
	* @desc compare 3 valeurs
	* @param int : un entier en haut, un à gauche et un en diagonale
	* @return return_value : la valeur la plus grande 
	**/
	
	int return_value = left ; //la valeur la plus grande est celle de gauche par défaut 
	
	if(up>return_value){return_value=up;}
	if(diag>return_value){return_value=diag;}
	if(return_value<0){return_value=0;} //si la valeur est négative, on la fixe à 0
	
	return return_value;
}

void Smith_Waterman::locate_replace_max(const unsigned int index,const unsigned int value, unsigned int max_table[], unsigned int index_max_table[]) 
{
	
	/**
	* @desc  
	* @param int : une valeur à comparer et son index, une valeur max d'un tableau et son index 
	**/
	
	if(value > max_table[NUMBER_OF_MAX_SAVED-1])//Vérifie si la valeur à tester est supérieure à la plus petite valeur du tableau
	{
		//Possibilité de recherche dichotomique plus mais pas forcément intéressante si vraiment petit tableau
		max_table[NUMBER_OF_MAX_SAVED-1] = value ;
		index_max_table[NUMBER_OF_MAX_SAVED-1] = index ;
		int pos_found = NUMBER_OF_MAX_SAVED-2 ;
		while((value > max_table[pos_found]) and (pos_found>=0)) //parcourt le tableau tant que la valeur à tester est plus grande
		{
			max_table[pos_found+1] = max_table[pos_found] ;
			index_max_table[pos_found+1] = index_max_table[pos_found] ;
			max_table[pos_found] = value ; //Si une valeur du tableau est plus petite que value, value la remplace
			index_max_table[pos_found] = index ;
			--pos_found ; //Si une valeur du tableau plus grande que value n'a pas été trouvée, on monte dans le tableau
			
		}
	}
}

//Setters 
void Smith_Waterman::set_gap_opener(const unsigned int new_gap_opener){this->gap_opener=new_gap_opener;}
void Smith_Waterman::set_gap_extension(const unsigned int new_gap_extension){this->gap_extension = new_gap_extension;}



const void Smith_Waterman::display_information(Handle_Database* database)
{
	/**
	* @desc affiche des informations sur la database et le protéine query qu'on veut comparer
	* @param Handle_Database* : pointeur vers la database
	**/
	
		cout << "Database title: " ; 
		for(unsigned int i=0;i<strlen(database->get_title());++i){cout<<database->get_title()[i];}
		cout << endl ;
		
		//affiche le nombre de résidus pour les nombre de protéines de la database
		cout << "Database size: " << (int)(database->get_numbers_of_residues()) << " residues in "
			 << database->get_database_size()-1 << " sequences" << endl ; 
		
		//affiche le nombre de résidus de la plus longue protéine de la database
		cout<<"Longest db seq: " << (int)(database->get_prot_max_length()) << " residues"<<endl ;
		
		//affiche l'heure de création de la database
		cout << "Database creation  : " ;
		for(unsigned int i=0;i<strlen(database->get_timestamp());++i){cout<<database->get_timestamp()[i];}
		cout << endl ;
		
		//affiche le nom de la protéine query
		cout <<"Query name :  ";
		for(size_t i=0;i<this->query_protein_header->size();++i)
		{ printf("%c",this->query_protein_header->at(i));}
		cout << endl;
		
		//affiche la taille de la protéine query
		cout << "Query length : " << this->query_protein->size() << " residues" << endl << endl;
		
	
}

const void Smith_Waterman::display_max(unsigned int* max_saved, unsigned int* index_max_saved, Handle_Database* database)
{
	
	/**
	* @desc affiche les protéines ayant les meilleurs scores
	* @param int : pointeur vers les scores maximums sauvegardés ainsi que leur index, Handle_database* : pointeur vers la database
	**/
	
	cout << "Sequences producing significant alignments:" << endl;
	
	string name_display = " ";
	size_t limit_size = 60 ; 
	for(size_t i=0; i<NUMBER_OF_MAX_SAVED; ++i)
	{
		name_display = database->fetch_prot_header(index_max_saved[i]) ;
		if(name_display.size() > limit_size){name_display = name_display.substr(0,limit_size) +"..." ;} //affiche que les 60 premiers caractères du nom de la protéine i
		cout<< "Score :" <<  max_saved[i]  ; //affiche le score de la protéine i 
		cout << " index " << index_max_saved[i] << ":" << name_display << endl; //affiche l'index de la protéine i 
	}
}

