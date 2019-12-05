#include "myProtein.h"
#include "smith_waterman.h"
#include <cstring>

//Main prend un argument le chemin vers la database, vers la query et vers la matrice blossum

/*
 * 
 * Attention a la fin !!!!
 * 
 * 
 */

void display_vector(const vector<char>* subject)
{
	for(unsigned int i=0; i<subject->size();++i)
	{
		if(i%80==0){cout << endl;}
		cout << subject->at(i) ;
	}
	
	cout << endl;
}

int main(int argc, char* argv[])//Main du projet finale, prend trois parametres
{
	Handle_Database* database = new Handle_Database(argv[1]);
	myProtein* prot_query =  new myProtein(argv[2]);
	
	//Blossum 62 de bae mais faire en sorte de pouvoir toute les prendre
	Smith_Waterman* scoring_algorithm = new Smith_Waterman("BLOSUM62",prot_query->getSequence());
	
	scoring_algorithm->score_protein(database);
	
	delete scoring_algorithm;
	delete database;
	delete prot_query;
	
return 0;
}

