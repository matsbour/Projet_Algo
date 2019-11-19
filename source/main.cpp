#include "handle_database.h"
#include <cstring>

//Main prend un argument le chemin vers la database.psq et .phr

void display_header(Handle_Database* database) // Permet d afficher le header de la prot suivante
{
	database->update_next_protein_header();
	for(int i=0; i<database->get_prot_header_active()->size();++i)
	{
		if(i%80==0){cout << endl;}
		cout << database->get_prot_header_active()->at(i) ;
	}
	cout << endl;
}

int main(int argc, char* argv[])
{
	Handle_Database* database = new Handle_Database(argv[1]);
	database->update_next_protein_sequence();
	
	for(int i=0; i<10; ++i)
	{
		display_header(database);
	}
	
	delete database ;
	return 0 ;
}


