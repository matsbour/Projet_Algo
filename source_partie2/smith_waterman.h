#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include "handle_database.h"
using namespace std;

class Smith_Waterman
{
private:
	 //Permet de retrouver les valeurs necessaire a la construction de la matrice score 
	 //Accede a une valeur par (blossum_matrix->at('X')).at('Y')
	map<char,map<char,int>>* blossum_matrix;
	string* query_protein; // Contient la protein a chercher dans la database
	
	//Fct privee
	void build_blossum_matrix(const string filepath);//Fct appele dans le constructeur, remplis la matrix blossum
	
public:

	Smith_Waterman(const string filepath,string* query_protein); // Constucteur
	~Smith_Waterman();//Destructeur
	
	// Permet de retourner le score d une proteine de la database comparer a celle de la query
	unsigned int score_protein(Handle_Database* database);

};

