#include "handle_database.h" // include permet aussi d'amener toutes les fonctionnalités telles que vector,map...
#include "myProtein.h"
#include <thread>
#include <mutex>
using namespace std;

class Smith_Waterman
{
private:
	 //Permet de retrouver les valeurs necessaires à la construction de la matrice score 
	vector<vector<int>>* blossum_matrix;  // Dépend du prot_dictionnary pour la traduction de char to int
	vector<int>* query_protein; // Contient la protéine à chercher dans la database qui sera traduite char->int avec le prot_dico de handle_database
	string* query_protein_header ;
	int gap_opener ; // 11 par défaut, modifié si besoin par un set
	int gap_extension ; //1 par défaut, modifié si besoin par un set
	Handle_Database* database ; // Contient le pointeur vers la database
	map<char,int> prot_dictionnary ; //Dictionnaire pour lire le int associé à chaque char
	
	//Fonction privée
	void build_blossum_matrix(const string filepath);//Fct appelée dans le constructeur et remplit la matrix blosum
	void locate_replace_max(const unsigned int index,const unsigned int value, unsigned int max_table[], unsigned int index_max_table[]) ;
	void display_information() const;
	void display_max() const;
	int max_over_zero(int up, int left, int diag) const;
	
public:

	Smith_Waterman(const string filepath,myProtein* query_protein_ini, int gap_opener_penalty, int gap_extension_penalty,Handle_Database* database); // Constucteur
	~Smith_Waterman();//Destructeur
	
	// Permet de retourner le score d'une protéine de la database comparée à celle de la query
	void score_protein(int identifier);
	void setup_score_protein();
	
	//Setters
	void set_gap_opener(const unsigned int new_gap_opener);
	void set_gap_extension(const unsigned int new_gap_extension);

};

