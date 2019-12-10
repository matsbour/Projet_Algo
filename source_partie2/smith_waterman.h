#include "handle_database.h" // include permet aussi d amener toute les fonctionnalite tel vector,map...
using namespace std;

class Smith_Waterman
{
private:
	 //Permet de retrouver les valeurs necessaire a la construction de la matrice score 
	 //Accede a une valeur par (blossum_matrix->at('X')).at('Y')
	//map<int,map<int,int>>* blossum_matrix;  // Depend du prot_dictionnary pour la trad char to int
	vector<vector<int>>* blossum_matrix;  // Depend du prot_dictionnary pour la trad char to int
	vector<int>* query_protein; // Contient la protein a chercher dans la database qui sera traduite char->int avec le prot_dico de handle_database
	int gap_opener ; // 11 par defaut modifie si besoin par un set
	int gap_extension ; //1 par defaut modifier si besoin par un set
	
	map<char,int> prot_dictionnary ; //Dictionnaire pour lire le int associer a chaque char
	
	//Fct privee
	void build_blossum_matrix(const string filepath);//Fct appele dans le constructeur, remplis la matrix blossum
	void locate_replace_max(const unsigned int index,const unsigned int value, unsigned int max_table[], unsigned int index_max_table[]) ;
	
public:

	Smith_Waterman(const string filepath,string* query_protein_ini, int gap_opener_penalty, int gap_extension_penalty ); // Constucteur
	~Smith_Waterman();//Destructeur
	
	// Permet de retourner le score d une proteine de la database comparer a celle de la query
	unsigned int score_protein(Handle_Database* database);
	int max_over_zero(int up, int left, int diag) const;
	
	//Get et set
	void set_gap_opener(const unsigned int new_gap_opener);
	void set_gap_extension(const unsigned int new_gap_extension);

};

