#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class myProtein {
	
	public:
	
	myProtein(const string proteinfile);
	
	
	 // Getters
	 
	string* getSequence();
	string* getHeader();
	const int getSize() const;
	
	
	private:

	
	// Variables d'instance
	string header;
	string sequence;
	int size;
   

	// Déclaration de la fonction output_header() qui extrait le header du fichier avec la data base 
	string output_header(const string proteinfile);
	
	// Déclaration de la fonction output_sequence() qui extrait la séquence du fichier avec la data base 
	string output_sequence(const string proteinfile);
	
		
};
