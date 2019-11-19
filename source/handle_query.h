#include <iostream>
#include <string>
using namespace std;

/*Classe contenant toute les fonctions necessaires a la gestion de
 * la base de donnée et de la query au format FASTA qui sera fournis.
 */
 
class Handle_query
{
private:
	 string database_filepath ;
	 string query_description ; //Contiendra les informations de l entete du fichier
	 string query_prot_chain ; //Contiendra la chaine a comparer
	 
public:
	 Handle_query(string query_file_path, string database_filepath); //Constructeur passage du fichier a ouvrir en parametre
	 
	 //Set et get
	 string get_prot_chain() const;
	 string get_query_description() const;
	 
};
 

