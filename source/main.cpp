#include "handle_database.h"
#include <cstring>
#include "myProtein.h"

//Main prend un argument le chemin vers la database

void display_vector(const vector<char>* subject)
{
	for(unsigned int i=0; i<subject->size();++i)
	{
		if(i%80==0){cout << endl;}
		cout << subject->at(i) ;
	}
	cout << endl;
}

int main(int argc, char* argv[])
{
	Handle_Database* database = new Handle_Database(argv[1]);
	
	cout << endl << "From fetch : "<< endl ;
	database->fetch_prot_header(150);
	database->fetch_prot_sequence(150);
	display_vector(database->get_prot_header_active());
	display_vector(database->get_prot_active());
	
	
	delete database ;
	return 0 ;
}


int main() 
{ 
    // Creating object of the class 
    myProtein object ; 
  
    // Extracting the header 
   cout << object.header << endl;
   cout << object.sequence << endl;

  
    return 0; 
}

