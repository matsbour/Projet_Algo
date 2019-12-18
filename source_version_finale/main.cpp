#include "smith_waterman.h"
#include <cstring>

//Valeurs constantes par défaut pour l'ouverture d'un gap et l'extension d'un gap
const int DEFAULT_GAP_OPENER = 11 ;
const int DEFAULT_GAP_EXTENSION = 1 ;

/*Main du projet final, prend obligatoirement 2 paramètres :
	- Chemin vers la database fasta
	- Chemin vers le fichier fasta de la query
	
De plus, il peut il y avoir 3 autres paramètres optionnels :
	- Chemin vers une matrice Blosum
	- unsigned int gap_opener
	- unsigned int gap_extension
!!Si un nombre négatif est donné pour l'un des gaps, il sera ignoré et si une chaine '1' est donnée pour Blosum il sera ignoré!!
*/

int main(int argc, char* argv[])//Main du projet finale, prend trois parametres
{
	/*Prise en compte des paramètres sous forme -option Paramètre*/
	char* daBase=NULL;
	char* protein_in=NULL;
	char* blosum=NULL;
	int gap_open=DEFAULT_GAP_OPENER, gap_ext=DEFAULT_GAP_EXTENSION;
	int obligatory[2]={0,0};
	int optionnals[3]={0,0,0};
	int arg_num=0;
	while(argv[arg_num] != NULL){
		if (argv[arg_num][0]=='-' && argv[arg_num+1] != NULL && argv[arg_num+1][0] != '-'){
			if(argv[arg_num][1]=='d'){
				daBase=argv[arg_num+1];
				obligatory[0]=1;
			}
			else if(argv[arg_num][1] == 'i' && argv[arg_num][2] =='n'){
				protein_in=argv[arg_num+1];
				obligatory[1]=1;
			}
			else if(argv[arg_num][1] == 'b'){
				blosum=argv[arg_num+1];
				optionnals[0]=1;
			}
			else if(argv[arg_num][1]=='g' && argv[arg_num][2]=='o' && argv[arg_num][3]=='p'){
				gap_open=atoi(argv[arg_num+1]);
				optionnals[1]=1;
			}
			else if(argv[arg_num][1]=='g' && argv[arg_num][2]=='e' && argv[arg_num][3]=='x'){
				gap_ext=atoi(argv[arg_num+1]);
				optionnals[2]=1;
			}
		}
		arg_num++;
	}	
	if(!obligatory[0] || ! obligatory[1]){
		cout << "Erreur nombre de paramètres invalide" << endl;
		cout << "./projet -d DataBase.fasta -in Protein.fasta [-b Blosum] [-gop Gap_opener] [-gex Gap_extension]" << endl;
		exit(1);
	}
	/*------------------------------------------------*/
	Handle_Database* database = new Handle_Database(daBase);
	myProtein* prot_query =  new myProtein(protein_in);
	Smith_Waterman* scoring_algorithm = NULL;

	//On considère les 8 cas où on change ou non les valeurs par défaut de le chemin vers une matrice Blosum, le gap opener et le gap extension 
	
	switch(argc)
	{
	case(3): //Cas où on a les 3 valeurs d'arguments par défaut
		scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query,DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);
	break;
	case(4):
	{
		if(!optionnals[0]){scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query,DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
		
		//chemin différent de Blosum 62
		else{scoring_algorithm = new Smith_Waterman(blosum,prot_query,DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
	}
	break;
	case(5):
	{
		if(!optionnals[0])
		{
			if(gap_open<0){scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query,DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
			
			//chemin différent de Blosum 62 et différente valeur de gap opener
			else{scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query,gap_open, DEFAULT_GAP_EXTENSION);}
		}
		else
		{
			if(gap_open<0){scoring_algorithm = new Smith_Waterman(blosum,prot_query,DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
			else{scoring_algorithm = new Smith_Waterman(blosum,prot_query,gap_open, DEFAULT_GAP_EXTENSION);}
		}

	}
	break;
	case(6):
	{
		if(!optionnals[0])
		{
			if(gap_open<0)
			{
				if(gap_ext<0){scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query,DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
				
				//différente valeur de gap extension
				else{scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query,DEFAULT_GAP_OPENER, gap_ext);}
			}
			else
			{
				if(gap_ext<0){scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query,gap_open, DEFAULT_GAP_EXTENSION);}
				
				//différentes valeurs de gap opener et gap extension
				else{scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query,gap_open, gap_ext);}
			}
		}
		else
		{
			if(gap_open<0)
			{
				if(gap_ext<0){scoring_algorithm = new Smith_Waterman(blosum,prot_query,DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
				
				//différent chemin de blosum 62 et différent gap extension
				else{scoring_algorithm = new Smith_Waterman(blosum,prot_query,DEFAULT_GAP_OPENER, gap_ext);}
			}
			else
			{
				if(gap_ext<0){scoring_algorithm = new Smith_Waterman(blosum,prot_query,gap_open, DEFAULT_GAP_EXTENSION);}
				
				//les 3 paramètres sont différents de ceux par défaut 
				else{scoring_algorithm = new Smith_Waterman(blosum,prot_query,gap_open, gap_ext);}
			}
		}
	}
	break;
	default: 
	{
		cout << "Erreur nombre de parametre invalide : " << argc << endl ;
		exit(1);
	}
	}
	
	
	scoring_algorithm->score_protein(database); //calcule le score de la protéine
	
	delete scoring_algorithm; //Libère la mémoire 
	delete database;
	delete prot_query;
	
return 0;
}

