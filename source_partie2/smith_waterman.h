#include "handle_database.h" // include permet aussi d amener toute les fonctionnalite tel vector,map...
using namespace std;

class Smith_Waterman
{
private:
	 //Permet de retrouver les valeurs necessaire a la construction de la matrice score 
	 //Accede a une valeur par (blossum_matrix->at('X')).at('Y')
	map<char,map<char,int>>* blossum_matrix; 
	const string* query_protein; // Contient la protein a chercher dans la database
	unsigned int gap_opener ; // 11 par defaut modifie si besoin par un set
	unsigned int gap_extension ; //1 par defaut modifier si besoin par un set
	
	//Fct privee
	void build_blossum_matrix(const string filepath);//Fct appele dans le constructeur, remplis la matrix blossum
	void locate_replace_max(const unsigned int index,const unsigned int value, unsigned int max_table[], unsigned int index_max_table[]) ;
	
public:

	Smith_Waterman(const string filepath,string* query_protein); // Constucteur
	~Smith_Waterman();//Destructeur
	
	// Permet de retourner le score d une proteine de la database comparer a celle de la query
	unsigned int score_protein(Handle_Database* database);
	int max_over_zero(int up, int left, int diag) const;
	
	//Get et set
	void set_gap_opener(const unsigned int new_gap_opener);
	void set_gap_extension(const unsigned int new_gap_extension);

};

