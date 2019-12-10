#include "myProtein.h"
#include "smith_waterman.h"
#include <cstring>

//Valeurs constantes par défaut pour l'ouverture d'un gap et l'extension d'un gap
const int DEFAULT_GAP_OPENER = 11 ;
const int DEFAULT_GAP_EXTENSION = 1 ;

/*Attention a la fin de la databse verif si aucun prbl*/

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
	if(argc < 3) //Verification si le nb de parametre est valide
	{
		cout << "Erreur nombre de parametre invalide : " << argc << endl ;
		exit(1);
	}
	
	Handle_Database* database = new Handle_Database(argv[1]);
	myProtein* prot_query =  new myProtein(argv[2]);
	Smith_Waterman* scoring_algorithm = NULL;
	
	//On considère les 8 cas où on change ou non les valeurs par défaut de le chemin vers une matrice Blosum, le gap opener et le gap extension 
	
	switch(argc) 
	{
	case(3): //Cas où on a les 3 valeurs d'arguments par défaut
		scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query->getSequence(),DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);
	break;
	case(4): 
	{
		if(atoi(argv[3]) == 1){scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query->getSequence(),DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
		
		//chemin différent de Blosum 62
		else{scoring_algorithm = new Smith_Waterman(argv[3],prot_query->getSequence(),DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
	}
	break;
	case(5):
	{
		if(atoi(argv[3]) == 1) 
		{	
			if(atoi(argv[4])<0){scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query->getSequence(),DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
			else{scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query->getSequence(),DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
		}
		else
		{
			if(atoi(argv[4])<0){scoring_algorithm = new Smith_Waterman(argv[3],prot_query->getSequence(),DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
			
			//chemin différent de Blosum 62 et différente valeur de gap opener
			else{scoring_algorithm = new Smith_Waterman(argv[3],prot_query->getSequence(),atoi(argv[4]), DEFAULT_GAP_EXTENSION);}
		}

	}
	break;
	case(6):
	{
		if(atoi(argv[3]) == 1)
		{
			if(atoi(argv[4])<0)
			{
				if(atoi(argv[5])<0){scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query->getSequence(),DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
				
				//différente valeur de gap extension
				else{scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query->getSequence(),DEFAULT_GAP_OPENER, atoi(argv[5]));}
			}
			else
			{
				//différente valeur de gap opener
				if(atoi(argv[5])<0){scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query->getSequence(),atoi(argv[4]), DEFAULT_GAP_EXTENSION);}
				
				//différentes valeurs de gap opener et gap extension
				else{scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query->getSequence(),atoi(argv[4]), atoi(argv[5]));}
			}
		}
		else
		{
			if(atoi(argv[4])<0)
			{
				if(atoi(argv[5])<0){scoring_algorithm = new Smith_Waterman(argv[3],prot_query->getSequence(),DEFAULT_GAP_OPENER, DEFAULT_GAP_EXTENSION);}
				
				//différent chemin de blosum 62 et différent gap extension
				else{scoring_algorithm = new Smith_Waterman(argv[3],prot_query->getSequence(),DEFAULT_GAP_OPENER, atoi(argv[5]));}
			}
			else
			{
				if(atoi(argv[5])<0){scoring_algorithm = new Smith_Waterman(argv[3],prot_query->getSequence(),atoi(argv[4]), DEFAULT_GAP_EXTENSION);}
				
				//les 3 paramètres sont différents de ceux par défaut 
				else{scoring_algorithm = new Smith_Waterman(argv[3],prot_query->getSequence(),atoi(argv[4]), atoi(argv[5]));}
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

